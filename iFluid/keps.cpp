/*******************************************************************
 * 	keps.cpp	
 *******************************************************************/
/*
 //TODO: ??? I assume this is related to incomplete parallelization?
 //
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include "keps.h"

static void printField2d(double*,const char*,int*,int*,int*);
static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};
static int find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);


//--------------------------------------------------------------------------
// 		KE_CARTESIAN
//--------------------------------------------------------------------------

KE_CARTESIAN::KE_CARTESIAN(Front &front)
    : front(&front)
{}

KE_CARTESIAN::~KE_CARTESIAN()
{}

bool KE_CARTESIAN::activated = false;
void KE_CARTESIAN::activateKE() {activated = true;}
void KE_CARTESIAN::deactivateKE() {activated = false;}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup edges:		
// 2) setup cell_center
//---------------------------------------------------------------
void KE_CARTESIAN::initMesh(void)
{
	int i,j,k,index;
	double crds[MAXD];
	int icoords[MAXD];
	int num_cells;
	int cell_index;
	IF_RECTANGLE       rectangle;

	// init vertices,edges & cell_center
	if (debugging("trace")) printf("Entering initMesh()\n");
	FT_MakeGridIntfc(front);
	setDomain();

	num_cells = 1;
	for (i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	}
	cell_center.insert(cell_center.end(),num_cells,rectangle);
	
	// setup vertices
	// left to right, down to up
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
            index = d_index1d(i,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    }
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
            index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
	    	crds[2] = top_L[2] + top_h[2]*k;
            index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}

	setComponent();
	FT_FreeGridIntfc(front);
	if (debugging("trace")) printf("Leaving initMesh()\n");
}

void KE_CARTESIAN::setComponent(void)
{
    static STATE *state = NULL;
    double *coords;
    int *icoords;
    int size = (int)cell_center.size();
    HYPER_SURF_ELEMENT *hse;
    HYPER_SURF *hs;
    double t[MAXD],point[MAXD];
    int n;

    if (state == NULL)
        FT_ScalarMemoryAlloc((POINTER*)&state, sizeof(STATE));

    for (int i = 0; i < size; i++)
    {
        coords = cell_center[i].coords;
        if (cell_center[i].comp != -1 &&
            cell_center[i].comp != top_comp[i])
        {
            if (FT_FindNearestIntfcPointInRange(front,top_comp[i],coords,
                            INCLUDE_BOUNDARIES,point,t,&hse,&hs,2))
            {
                if (!FrontNearestIntfcState(front,coords,top_comp[i],
                            (POINTER)state))
                {
                    (void) printf("In setComponent()\n");
                    (void) printf("FrontNearestIntfcState() failed\n");
                    (void) printf("old_comp = %d new_comp = %d\n",
                                    cell_center[i].comp,top_comp[i]);
                    clean_up(ERROR);
                }
            }
            else
            {
                double temp_nb = 0.0;
                int ic[MAXD],index;
                icoords = cell_center[i].icoords;
                n = 0;

                for (int ii = 0; ii < dim; ++ii)
                {
                    for (int jj = 0; jj < dim; ++jj)
                        ic[jj] = icoords[jj];

                    ic[ii] = (icoords[ii] == 0) ? 0 : icoords[ii] - 1;
                    index = d_index(ic,top_gmax,dim);
                    
                    if (cell_center[index].comp != -1 &&
                        cell_center[index].comp == top_comp[index])
                    {
                        n++;
                    }
                    
                    ic[ii] = (icoords[ii] == top_gmax[ii]) ? top_gmax[ii]
                                    : icoords[ii] + 1;
                    
                    index = d_index(ic,top_gmax,dim);
                    
                    if (cell_center[index].comp != -1 &&
                        cell_center[index].comp == top_comp[index])
                    {
                        //temp_nb += field->temperature[index];
                        n++;
                    }
                }
                //field->temperature[i] = temp_nb/n;
            }
    
        }
        cell_center[i].comp = top_comp[i];
    }
}	/* end setComponent */

void KE_CARTESIAN::setInitialCondition(void)
{
	int i;
	double coords[MAXD];
	INTERFACE *intfc = front->interf;
	POINT *p;
    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	int c;
	short unsigned int seed[3] = {2,72,7172};

	FT_MakeGridIntfc(front);
    setDomain();

    //TODO: find inlet velocity by checking rect boundaries in order
    //      to set U0 (U_infty) which is used to limit TKE, epsilon, and mu_t.

    double nu = eqn_params->mu/eqn_params->rho;

    double k_inlet = 1.5*sqr(eqn_params->U0*eqn_params->I);
    double eps_inlet = eqn_params->Cmu*sqr(k_inlet)/nu/eqn_params->ViscRatioFar;
    //TODO: Use near or far visc ratio for inlet boundary?

    //save inlet boundary conditions for use in solver
    eqn_params->k_inlet = k_inlet;
    eqn_params->eps_inlet = eps_inlet;

    double k0 = sqr(eqn_params->mu0/eqn_params->l0/eqn_params->rho);
    double eps0 = eqn_params->Cmu*pow(k0,1.5)/eqn_params->l0;
        //double k0 = 1.5*sqr(eqn_params->U0*eqn_params->I);
        //double eps0 = eqn_params->Cmu*sqr(k0)/nu/eqn_params->ViscRatioFar;

    //save the initial conditions for activation time
    eqn_params->k0 = k0;
    eqn_params->eps0 = eps0;

	for (i = 0; i < cell_center.size(); i++)
	{
	    field->k[i] = 0.0;
	    field->eps[i] = 0.0;
        field->mu_t[i] = 0.0;
        
        /*
	    c = top_comp[i];
	    getRectangleCenter(i,coords);
	    field->k[i] = k0;
	    field->eps[i] = eps0;
	    if (keps_model == REALIZABLE)
            field->Cmu[i] = eqn_params->Cmu;
	    field->mu_t[i] = eqn_params->mu0;
        */
	}
	printf("k0 = %e, eps0 = %e\n",k0,eps0);
	printf("mu0 = %e\n",eqn_params->mu0);
}	/* end setInitialCondition */

void KE_CARTESIAN::applyInitialConditions()
{
    double k0 = eqn_params->k0;
    double eps0 = eqn_params->eps0;
    double mu0 = eqn_params->mu0;
    double Cmu = eqn_params->Cmu;

    double rho = eqn_params->rho;
    double nu = eqn_params->mu/rho;

	for (int i = 0; i < cell_center.size(); ++i)
	{
	    field->k[i] = 0.0;
	    field->eps[i] = 0.0;
            //field->k[i] = k0;
            //field->eps[i] = eps0;
	    
        //double nu_t = std::max(Cmu*sqr(k0)/eps0,0.0001*nu);
        //field->mu_t[i] = nu_t*rho;
	        //field->mu_t[i] = mu0;
        field->mu_t[i] = 0.0;
        
        if (keps_model == REALIZABLE)
            field->Cmu[i] = Cmu;
	}
}

/*compute lift force*/
//TODO: verify computation
//not in use..
void KE_CARTESIAN::computeLiftDrag(Front* front)
{
	double force[MAXD];
	double torque;
	static boolean first = YES;
	boolean skip = NO;
	char fname[200];
	FILE* vfile;
	CURVE **c;
	INTERFACE *intfc = front->interf;
	int dim = front->rect_grid->dim;
	double dt = front->dt;
	double tol[MAXD];
	double ref_p[MAXD];
	int i;

	for (i = 0; i < dim; i++)
	    tol[i]  = 2*front->rect_grid->h[i];
	
	for (c = intfc->curves; c && *c; ++c)
	{
	    for (i = 0; i < dim; i++)
	    {
            ref_p[i] = Coords((*c)->first->start)[i];
	        if (ref_p[i] < front->rect_grid->L[i] + 2.0*tol[i] ||
                ref_p[i] > front->rect_grid->U[i] - 2.0*tol[i])
            {
                skip = YES;
            }
	    }

	    if (skip) continue;

	    if (wave_type(*c) == NEUMANN_BOUNDARY || wave_type(*c) == ELASTIC_BOUNDARY)
            ifluid_compute_force_and_torque(front,Hyper_surf(*c),dt,force,&torque);
	}
	sprintf(fname,"%s/force",OutName(front));
	if (first)
	{
	    first = NO;
	    vfile = fopen(fname,"w");
	    fclose(vfile);		
	}
	vfile = fopen(fname,"a");
	fprintf(vfile,"%f %f %f\n",front->time,force[0],force[1]);
	fclose(vfile);
}

void KE_CARTESIAN::setIndexMap(COMPONENT sub_comp)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];
	int count;

	if (debugging("trace")) printf("Entering setIndexMap()\n");
	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 1:
		count = (imax - imin + 1);
	    	FT_VectorMemoryAlloc((POINTER*)&i_to_I,top_gmax[0]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_i,count,1,INT);
	    	break;
	    case 2:
		count = (imax - imin + 1)*(jmax - jmin + 1);
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,
					top_gmax[1]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ij,count,2,INT);
	    	break;
	    case 3:
		count = (imax - imin + 1)*(jmax - jmin + 1)*(kmax - kmin + 1);
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,
					top_gmax[1]+1,top_gmax[2]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ijk,count,3,INT);
	    	break;
	    }
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = (lbuf[i] != 0) ? lbuf[i] : 1;
	    uubuf[i] = (ubuf[i] != 0) ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    i_to_I[i] = index + ilower;
	    	    index++;
		}
		else
		    i_to_I[i] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)i_to_I);
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ij_to_I[i][j] = index + ilower;
		    I_to_ij[index][0] = i;
                    I_to_ij[index][1] = j;
	    	    index++;
		}
		else
		    ij_to_I[i][j] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ijk_to_I[i][j][k] = index + ilower;
		    I_to_ijk[index][0] = i;
                    I_to_ijk[index][1] = j;
                    I_to_ijk[index][2] = k;
	    	    index++;
		}
		else
		    ijk_to_I[i][j][k] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
	if (debugging("trace")) printf("Leaving setIndexMap()\n");

}	/* end setIndexMap */

void KE_CARTESIAN::computeAdvection()
{
	int i;
	COMPONENT sub_comp[2];

	sub_comp[0] = SOLID_COMP;
	sub_comp[1] = LIQUID_COMP2;
	
	for (i = 0; i < 2; ++i)
    {
        if(sub_comp[i] == SOLID_COMP) continue;

	   setGlobalIndex(sub_comp[i]);
	   computeAdvectionK(sub_comp[i]);
	   
       //NOTE: ONLY RNG model works
	   switch (keps_model)
	   {
		case STANDARD:
	       	    computeAdvectionE_STD(sub_comp[i]);
		    break;
		case RNG:
	       	    computeAdvectionE_RNG(sub_comp[i]);
		    break;
		case REALIZABLE:
	       	    computeAdvectionE_REAL(sub_comp[i]);
		    break;
		default:
		    printf("Unknown k-eps model\n");
		    clean_up(ERROR);
	   }
	}
}

//TODO: incomplete
//not in use anywhere
void KE_CARTESIAN::findBdryPoint()
{
	HYPER_SURF *hs;
        STATE *intfc_state;
        INTERFACE* grid_intfc = front->grid_intfc;
        int icoords[MAXD];
	double crx_coords[MAXD];
	int index,i,j,l,m,nc;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	boolean fr_crx_grid_seg;
	COMPONENT comp;
	static CRXING *crxs[10];

	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
		icoords[0] = i; icoords[1] = j;
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		if (comp != LIQUID_COMP2)
                    continue;
		for (l = 0; l < dim; l++)
	 	for (m = 0; m < 2; m++)
		{
		    nc = GridSegCrossing(crxs,icoords,dir[l][m],grid_intfc);
		}
	}
}

void KE_CARTESIAN::chienComputeK(COMPONENT sub_comp)
{
    int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
    int gmin[MAXD],ipn[MAXD];
    double crx_coords[MAXD];
	double nor[MAXD];
    double K0,K_nb,D,lambda,coeff,coeff_nb,rhs;
    COMPONENT comp;
    
    PETSc solver;
    double *x;
    int num_iter = 0;
    double residual = 0.0;
    
    boolean fr_crx_grid_seg;
    const GRID_DIRECTION dir[3][2] =
            {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    
    double *K = field->k;
	double *Pk = field->Pk;
	double *mu_t = field->mu_t;
    double *K_prev = field->k_prev;
	double *gamma = field->gamma;

	double *dist_array = field->dist;
    double **nor_array = field->nor;
    double **f_surf = field->f_surf;

	double Cmu = eqn_params->Cmu;
	double delta_k = eqn_params->delta_k;
	double rho = eqn_params->rho;
	double nu = eqn_params->mu/eqn_params->rho;
	double y_p = eqn_params->y_p;
	
    double v[MAXD],v_wall[MAXD],crds_wall[MAXD],k_wall;
    double eta;
	double Ut,Ut_old;
	int    niter;

	/*For boundary state*/
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *intfc_state;
	INTERFACE* grid_intfc = front->grid_intfc;
	double coords[MAXD];
	boolean if_adj_pt;
	double point[MAXD],t[MAXD];

    start_clock("chienComputeK");
    if (debugging("trace")) printf("Entering chienComputeK()\n");

    for (i = 0; i < dim; ++i) gmin[i] = 0;

    setIndexMap(sub_comp);
    if (debugging("trace"))
    {
        int domain_size = 1;
        printf("ilower = %d  iupper = %d\n",ilower,iupper);
        for (i = 0; i < dim; ++i)
        printf("domain_size = %d\n",domain_size);
    }

    start_clock("set_coefficients");

    switch(dim)
    {
    case 2:

        solver.Create(ilower, iupper-1, 5, 5);

        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            ic = d_index2d(i,j,top_gmax);
            I = ij_to_I[i][j];
            
            //Save K for use in wall functions for epsilon
            K_prev[ic] = K[ic];
            //TODO: function for copying K and computing gamma?
            //      Or could put in computeMuTurb()??

            comp = top_comp[ic];
            if (comp != sub_comp) continue;

            K0 = K[ic];
            rhs = K0 + m_dt*Pk[ic];

            //Update the linearization parameter gamma
            if (fabs(mu_t[ic]) > MACH_EPS)
            {
                gamma[ic] = std::max(Cmu*K0*rho/mu_t[ic],0.0);
            }
            else
            {
                gamma[ic] = 0.0;
            }

            double alpha = gamma[ic] + 2.0*nu/dist_array[ic];

            coeff = 1.0 + m_dt*alpha;
                //coeff = 1.0 + m_dt*gamma[ic];


            if (isinf(coeff) || isnan(coeff) || isinf(rhs) || isnan(rhs))
            {
                printf("In chienComputeK(): ");
                printf("icoords[%d %d], index = %d\n",i,j,ic);
                printf("coeff=%f, K=%e, E=%e, mu_t=%e, Pk=%e\n",
                coeff,K0,field->eps[ic],mu_t[ic],Pk[ic]);
                clean_up(ERROR);
            }

            for (l = 0; l < dim; ++l) v[l] = 0.0;

            if (field->vel != NULL)
            {
                for (l = 0; l < dim; ++l)
                    v[l] = field->vel[l][ic];
            }
            
            D = nu + mu_t[ic]/eqn_params->delta_k/rho;

            for (l = 0; l < dim; ++l)
            {
                lambda = D*m_dt/sqr(top_h[l]);
                eta = v[l]*m_dt/(top_h[l]); //upwind difference
                double eta_p = std::max(eta, 0.0);
                double eta_m = std::min(eta, 0.0);
                
                coeff += eta_p - eta_m;
                coeff += 2.0*lambda;

                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index2d(ipn[0],ipn[1],top_gmax);
                    I_nb = ij_to_I[ipn[0]][ipn[1]];
            
                    //This version returns provides the hypersurf element hse needed for wall stress
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing2(front,
                            icoords,dir[l][m],comp,(POINTER*)&intfc_state,
                            &hs,&hse,crx_coords);
                    
                    if (!fr_crx_grid_seg) 
                    {
                        coeff_nb = -lambda;
                        coeff_nb += (m == 0) ? -eta_p : eta_m;
                        solver.Add_A(I,I_nb,coeff_nb);
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                             wave_type(hs) == ELASTIC_BOUNDARY)
                    {

                        K_nb = 0.0;
                            //rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb);

                        //Compute wall tangential shear stress (force/length for 2d)
                        double tau_wall[MAXD] = {0.0};
                        computeTangentialStress(icoords,tau_wall);
                        //TODO: update y+ here??
                        
                        int num_bonds;
                        BOND* bonds[10];
                        BondAndNeighbors(hse,hs,&num_bonds,bonds,3);//1 bond on each side for total of 3 bonds
                        //Try other orders 4,5 etc.
                        
                            //BOND* nearest_bond = Bond_of_hse(hse);
                            //double len = bond_length(nearest_bond);

                        double fwall[MAXD] = {0.0};
                        for (int ib = 0; ib < num_bonds; ++ib)
                        {
                            double len = bond_length(bonds[ib]);
                            for (int kk = 0; kk < dim; ++kk)
                                fwall[kk] += tau_wall[kk]*len;
                        }
                        
                        double* force_wall = intfc_state->shear_force;

                        for (int kk = 0; kk < dim; ++kk)
                        {
                            force_wall[kk] = fwall[kk];
                            f_surf[kk][ic] += -1.0*force_wall[kk]/rho;
                        }

                        //TODO: For elastic boundaries force_wall should be coupled
                        //      to the interface/canopy somehow... compute acceleration with it?
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: convective outlet boundary similar
                        //      to iF_flowThroughBoundaryState() needed?
                        if (boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            K_nb = K0;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                                //rhs += lambda*K_nb - (pow(-1,m+1)*eta)*K_nb;
                        }
                        else
                        {
                            //INLET
                            K_nb = eqn_params->k_inlet;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb);
                        }
                    }
                    else
                    {
                        printf("Unknown boundary condition! \n");
                        clean_up(ERROR);
                    }
                }
            }

            solver.Add_A(I,I,coeff);
            solver.Add_b(I,rhs);
        }
        break;

    case 3:

        solver.Create(ilower, iupper-1, 7, 7);
    
        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            I = ijk_to_I[i][j][k];
            ic = d_index3d(i,j,k,top_gmax);

            //Save K for use in wall functions for epsilon
            K_prev[ic] = K[ic];
            //TODO: function for copying K and computing gamma?
            //      Or could put in computeMuTurb()??

            comp = top_comp[ic];
            if (comp != sub_comp) continue;

            K0 = K[ic];
            Cmu = eqn_params->Cmu;
            rhs = K0 + m_dt*Pk[ic];

            //std and rng, not realizable
            if (fabs(mu_t[ic]) > MACH_EPS)
            {
                gamma[ic] = std::max(Cmu*K0*rho/mu_t[ic],0.0);
            }
            else
            {
                gamma[ic] = 0.0;
            }
                
            double alpha = gamma[ic] + 2.0*nu/dist_array[ic];

            coeff = 1.0 + m_dt*alpha;
                //coeff = 1.0 + m_dt*gamma[ic];

            
            if (isinf(coeff) || isnan(coeff))
            {
                printf("In computeAdvectionK(): ");
                printf("coeff=%f, K=%e, E=%e, mu_t=%e, Pk=%e\n",
                coeff,K0,field->eps[ic],mu_t[ic],Pk[ic]);
                clean_up(ERROR);
            }

            for (l = 0; l < dim; ++l) v[l] = 0.0;

            if (field->vel != NULL)
            {
                for (l = 0; l < dim; ++l)
                    v[l] = field->vel[l][ic];
            }

            D = nu + mu_t[ic]/eqn_params->delta_k/rho;
            
            for (l = 0; l < dim; ++l)
            {
                lambda = D*m_dt/sqr(top_h[l]);
                eta = v[l]*m_dt/(top_h[l]); //upwind difference
                double eta_p = std::max(eta, 0.0);
                double eta_m = std::min(eta, 0.0);

                coeff += eta_p - eta_m;
                coeff += 2.0*lambda;

                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                    I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
            
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing2(front,
                            icoords,dir[l][m],comp,(POINTER*)&intfc_state,
                            &hs,&hse,crx_coords);
                    
                    if (!fr_crx_grid_seg) 
                    {
                        coeff_nb = -lambda;
                        coeff_nb += (m == 0) ? -eta_p : eta_m;
                        solver.Add_A(I,I_nb,coeff_nb);
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                             wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        K_nb = 0.0;
                            //rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 

                        //Compute wall tangential shear stress (force/area)
                        double tau_wall[MAXD] = {0.0};
                        computeTangentialStress(icoords,tau_wall);
                        //TODO: update y+ here??

                        int num_tris;
                        TRI* tris[30];//TODO: is 30 enough tris?
                        TriAndFirstRing(hse,hs,&num_tris,tris);
                        //TODO: Or use PointAndFirstRingTris();

                        double fwall[MAXD] = {0.0};
                        for (int it = 0; it < num_tris; ++it)
                        {
                            double area = tri_area(tris[it]);
                            for (int kk = 0; kk < dim; ++kk)
                                fwall[kk] += tau_wall[kk]*area;
                        }
                        
                        double* force_wall = intfc_state->shear_force;

                        for (int kk = 0; kk < dim; ++kk)
                        {
                            force_wall[kk] = fwall[kk];
                            f_surf[kk][ic] += -1.0*force_wall[kk]/rho;
                        }

                        //TODO: For elastic boundaries force_wall should be coupled
                        //      to the interface/canopy somehow... compute acceleration with it?
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: convective outlet boundary similar
                        //      to iF_flowThroughBoundaryState() needed?
                        if (boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            K_nb = K0;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                        }
                        else
                        {
                            //INLET
                            K_nb = eqn_params->k_inlet;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                        }
                    }
                    else
                    {
                        printf("Unknown boundary condition %d! \n",wave_type(hs));
                        clean_up(ERROR);
                    }
                }
            }

            solver.Add_A(I,I,coeff);
            solver.Add_b(I,rhs);
        }
    
        break;
    }

    stop_clock("set_coefficients");
    
    solver.SetMaxIter(10000);
    solver.SetTol(1e-8);
    
    start_clock("petsc_solve");
    solver.Solve();
    stop_clock("petsc_solve");

	solver.GetNumIterations(&num_iter);
    solver.GetResidualNorm(&residual);

    if (debugging("PETSc"))
    {
        (void) printf("KE_CARTESIAN::chienComputeK: "
                    "num_iter = %d, residual = %g \n",
                    num_iter, residual);
    }

	if (residual > 1)
	{
	    printf("KE_CARTESIAN::chienComputeK(): \
                Solution diverges!\n");
	    LOC(); clean_up(ERROR);
	}

    //TODO: not parallelized? size = iupper - ilower ??
    FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
    solver.Get_x(x);

    start_clock("scatter_data");
    switch (dim)
    {
    case 1:
        for (i = imin; i <= imax; i++)
        {
            I = i_to_I[i];
            ic = d_index1d(i,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    case 2:
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            ic = d_index2d(i,j,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    case 3:
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            ic = d_index3d(i,j,k,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    }

    scatMeshArray();
    
    switch (dim)
    {
    case 1:
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index1d(i,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                K[ic] = array[ic];
        }
        break;
    case 2:
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index2d(i,j,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                K[ic] = array[ic];
        }
        break;
    case 3:
        for (k = 0; k <= top_gmax[2]; ++k)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index3d(i,j,k,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                K[ic] = array[ic];
        }
        break;
    }
    stop_clock("scatter_data");
    FT_FreeThese(1,x);

    FT_ParallelExchGridVectorArrayBuffer(f_surf,front);

    if (debugging("trace")) printf("Leaving chienComputeK()\n");
    stop_clock("chienComputeK");
}       /* end chienComputeK */

void KE_CARTESIAN::computeAdvectionK(COMPONENT sub_comp)
{
    int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
    int gmin[MAXD],ipn[MAXD];
    double crx_coords[MAXD];
	double nor[MAXD];
    double K0,K_nb,D,lambda,coeff,coeff_nb,rhs;
    COMPONENT comp;
    
    PETSc solver;
    double *x;
    int num_iter = 0;
    double residual = 0.0;
    
    boolean fr_crx_grid_seg;
    const GRID_DIRECTION dir[3][2] =
            {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    
    double *K = field->k;
	double *Pk = field->Pk;
	double *mu_t = field->mu_t;
    double *K_prev = field->k_prev;
	double *gamma = field->gamma;

	double *dist_array = field->dist;
    double **nor_array = field->nor;
    double **f_surf = field->f_surf;

	double Cmu = eqn_params->Cmu;
	double delta_k = eqn_params->delta_k;
	double rho = eqn_params->rho;
	double nu = eqn_params->mu/eqn_params->rho;
	double y_p = eqn_params->y_p;
	
    double v[MAXD],v_wall[MAXD],crds_wall[MAXD],k_wall;
    double eta;
	double Ut,Ut_old;
	int    niter;

	//For boundary state
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *intfc_state;
	INTERFACE* grid_intfc = front->grid_intfc;
	double coords[MAXD];
	boolean if_adj_pt;
	double point[MAXD],t[MAXD];

    start_clock("computeAdvectionK");
    if (debugging("trace")) printf("Entering computeAdvectionK()\n");

    for (i = 0; i < dim; ++i) gmin[i] = 0;

    setIndexMap(sub_comp);
    if (debugging("trace"))
    {
        int domain_size = 1;
        printf("ilower = %d  iupper = %d\n",ilower,iupper);
        for (i = 0; i < dim; ++i)
        printf("domain_size = %d\n",domain_size);
    }

    start_clock("set_coefficients");

    switch(dim)
    {
    case 2:

        solver.Create(ilower, iupper-1, 5, 5);

        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            ic = d_index2d(i,j,top_gmax);
            I = ij_to_I[i][j];
            
            //Save K for use in wall functions for epsilon
            K_prev[ic] = K[ic];
            //TODO: function for copying K and computing gamma?
            //      Or could put in computeMuTurb()??

            comp = top_comp[ic];
            if (comp != sub_comp) continue;

            K0 = K[ic];
            rhs = K0 + m_dt*Pk[ic];

            //Update the linearization parameter gamma
            if (fabs(mu_t[ic]) > MACH_EPS)
            {
                gamma[ic] = std::max(Cmu*K0*rho/mu_t[ic],0.0);
            }
            else
            {
                gamma[ic] = 0.0;
            }
                
            coeff = 1.0 + m_dt*gamma[ic];

            //TODO: REALIZABLE AND STD disabled
            //
            //if (keps_model == REALIZABLE)
            //    coeff = 1.0 + m_dt*std::max(field->eps[ic],0.0);
            //else
            //    coeff = 1.0 + m_dt*std::max(Cmu*K0*rho/mu_t[ic],0.0);
            //

            if (isinf(coeff) || isnan(coeff) || isinf(rhs) || isnan(rhs))
            {
                printf("In computeAdvectionK(): ");
                printf("icoords[%d %d], index = %d\n",i,j,ic);
                printf("coeff=%f, K=%e, E=%e, mu_t=%e, Pk=%e\n",
                coeff,K0,field->eps[ic],mu_t[ic],Pk[ic]);
                clean_up(ERROR);
            }

            for (l = 0; l < dim; ++l) v[l] = 0.0;

            if (field->vel != NULL)
            {
                for (l = 0; l < dim; ++l)
                    v[l] = field->vel[l][ic];
            }
            
            D = nu + mu_t[ic]/eqn_params->delta_k/rho;

            for (l = 0; l < dim; ++l)
            {
                lambda = D*m_dt/sqr(top_h[l]);
                eta = v[l]*m_dt/(top_h[l]); //upwind difference
                double eta_p = std::max(eta, 0.0);
                double eta_m = std::min(eta, 0.0);
                
                coeff += eta_p - eta_m;
                coeff += 2.0*lambda;

                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index2d(ipn[0],ipn[1],top_gmax);
                    I_nb = ij_to_I[ipn[0]][ipn[1]];
            
                    //This version returns provides the hypersurf element hse needed for wall stress
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing2(front,
                            icoords,dir[l][m],comp,(POINTER*)&intfc_state,
                            &hs,&hse,crx_coords);
                    
                    ////fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    //        grid_intfc,icoords,dir[l][m],comp,
                    //        (POINTER*)&intfc_state,&hs,crx_coords);


                    if (!fr_crx_grid_seg) 
                    {
                        coeff_nb = -lambda;
                        coeff_nb += (m == 0) ? -eta_p : eta_m;
                            //coeff_nb = -lambda + (pow(-1,m+1)*eta);
                        solver.Add_A(I,I_nb,coeff_nb);
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                             wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        //NOTE: Would this require resolution of the viscous sublayer?
                            //setTKEatWall(icoords,l,m,comp,hs,intfc_state,K,&K_nb);

                        double v_slip[MAXD] = {0.0};
                        setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_slip);

                        //use wall function
                        double u_t =
                            std::max(pow(Cmu,0.25)*sqrt(std::max(K[ic],0.0)),
                                    Mag2d(v_slip)/y_p);

                        //
                        //use wall function
                        //double u_t = std::max(
                        //          pow(eqn_params->Cmu, 0.25)*sqrt(
                        //              std::max(field->k[ic], 0.0)),
                        //          Mag2d(intfc_state->vel)/eqn_params->y_p);
                        //

                        //TODO: Note that intfc_state->vel is not
                        //      vel tangential to the wall.
                        //
                        //      project out the normal component,
                        //      required by the free-slip condition u dot n = 0.
                        
                        K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                        rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb);
                            //rhs += lambda*K_nb - (pow(-1,m+1)*eta)*K_nb;

                        //Compute wall tangential shear stress (force/length for 2d)
                        double tau_wall[MAXD] = {0.0};
                        for (int kk = 0; kk < dim; ++kk)
                        {
                            tau_wall[kk] = u_t/y_p*v_slip[kk];
                        }
                        
                            //tau_wall[0] = u_t/y_p*v_slip[0];
                            //tau_wall[1] = u_t/y_p*v_slip[1];
                                //tau_wall[0] = -1.0*u_t/y_p*v_slip[0];
                                //tau_wall[1] = -1.0*u_t/y_p*v_slip[1];
                        
                        //TODO: Need more than one interface element?
                        //      see BondAndNeighbors() usage in surfaceTension()
                        
                            //BOND* nearest_bond = Bond_of_hse(hse);
                            //double len = bond_length(nearest_bond);
                        
                        int num_bonds;
                        BOND* bonds[10];
                        BondAndNeighbors(hse,hs,&num_bonds,bonds,3);//1 bond on each side for total of 3 bonds
                        //Try other orders 2, 3

                        double fwall[MAXD] = {0.0};
                        for (int ib = 0; ib < num_bonds; ++ib)
                        {
                            double len = bond_length(bonds[ib]);
                            for (int kk = 0; kk < dim; ++kk)
                                fwall[kk] += tau_wall[kk]*len;
                        }
                        
                        double* force_wall = intfc_state->shear_force;
                        //force_wall[0] = tau_wall[0]*len;
                        //force_wall[1] = tau_wall[1]*len;
                        //f_surf[0][ic] += -1.0*force_wall[0]/rho;
                        //f_surf[1][ic] += -1.0*force_wall[1]/rho;

                        for (int kk = 0; kk < dim; ++kk)
                        {
                            force_wall[kk] = fwall[kk];
                            f_surf[kk][ic] += -1.0*force_wall[kk]/rho;
                        }

                        //TODO: For elastic boundaries force_wall should be coupled
                        //      to the interface/canopy somehow... compute acceleration with it?
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: convective outlet boundary similar
                        //      to iF_flowThroughBoundaryState() needed?
                        if (boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            K_nb = K0;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                                //rhs += lambda*K_nb - (pow(-1,m+1)*eta)*K_nb;
                        }
                        else
                        {
                            //INLET
                            
                                ///K_nb = eqn_params->Cbc
                                // * (sqr(intfc_state->vel[0])
                                // +  sqr(intfc_state->vel[1]));
                            
                            K_nb = eqn_params->k_inlet;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb);
                                //rhs += lambda*K_nb - (pow(-1,m+1)*eta)*K_nb;
                        }
                    }
                    else
                    {
                        printf("Unknown boundary condition! \n");
                        clean_up(ERROR);
                    }
                }
            }

            solver.Add_A(I,I,coeff);
            solver.Add_b(I,rhs);
        }
        break;

    case 3:

        solver.Create(ilower, iupper-1, 7, 7);
    
        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            I = ijk_to_I[i][j][k];
            ic = d_index3d(i,j,k,top_gmax);

            //Save K for use in wall functions for epsilon
            K_prev[ic] = K[ic];
            //TODO: function for copying K and computing gamma?
            //      Or could put in computeMuTurb()??

            comp = top_comp[ic];
            if (comp != sub_comp) continue;

            K0 = K[ic];
            Cmu = eqn_params->Cmu;
            rhs = K0 + m_dt*Pk[ic];

            //std and rng, not realizable
            if (fabs(mu_t[ic]) > MACH_EPS)
            {
                gamma[ic] = std::max(Cmu*K0*rho/mu_t[ic],0.0);
            }
            else
            {
                gamma[ic] = 0.0;
            }
                
            coeff = 1.0 + m_dt*gamma[ic];

            
            if (isinf(coeff) || isnan(coeff))
            {
                printf("In computeAdvectionK(): ");
                printf("coeff=%f, K=%e, E=%e, mu_t=%e, Pk=%e\n",
                coeff,K0,field->eps[ic],mu_t[ic],Pk[ic]);
                clean_up(ERROR);
            }

            for (l = 0; l < dim; ++l) v[l] = 0.0;

            if (field->vel != NULL)
            {
                for (l = 0; l < dim; ++l)
                    v[l] = field->vel[l][ic];
            }

            D = nu + mu_t[ic]/eqn_params->delta_k/rho;
            
            for (l = 0; l < dim; ++l)
            {
                lambda = D*m_dt/sqr(top_h[l]);
                eta = v[l]*m_dt/(top_h[l]); //upwind difference
                double eta_p = std::max(eta, 0.0);
                double eta_m = std::min(eta, 0.0);

                coeff += eta_p - eta_m;
                coeff += 2.0*lambda;

                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                    I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
            
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing2(front,
                            icoords,dir[l][m],comp,(POINTER*)&intfc_state,
                            &hs,&hse,crx_coords);
                    
                    //
                    //fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    //        grid_intfc,icoords,dir[l][m],comp,
                    //        (POINTER*)&intfc_state,&hs,crx_coords);
                    //
                    
                    if (!fr_crx_grid_seg) 
                    {
                        coeff_nb = -lambda;
                        coeff_nb += (m == 0) ? -eta_p : eta_m;
                        solver.Add_A(I,I_nb,coeff_nb);
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                             wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        //setTKEatWall(icoords,l,m,comp,hs,intfc_state,K,&K_nb);

                        double v_slip[MAXD] = {0.0};
                        setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_slip);

                        //use wall function
                        double u_t =
                            std::max(pow(Cmu,0.25)*sqrt(std::max(K[ic],0.0)),
                                    Mag3d(v_slip)/y_p);

                        //
                        //use wall function
                        //double u_t = std::max(
                        //        pow(eqn_params->Cmu, 0.25)*sqrt(
                        //            std::max(field->k[ic], 0.0)),
                        //        Mag3d(intfc_state->vel)/eqn_params->y_p);
                        //
                        //TODO: Note that intfc_state->vel is not
                        //      vel tangential to the wall
                        
                        K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                        rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 

                        //Compute wall tangential shear stress (force/area)
                        double tau_wall[MAXD] = {0.0};
                        for (int kk = 0; kk < dim; ++kk)
                        {
                            tau_wall[kk] = u_t/y_p*v_slip[kk];
                        }

                        int num_tris;
                        TRI* tris[30];//TODO: is 30 enough tris?
                        TriAndFirstRing(hse,hs,&num_tris,tris);
                        //TODO: Or use PointAndFirstRingTris();

                        double fwall[MAXD] = {0.0};
                        for (int it = 0; it < num_tris; ++it)
                        {
                            double area = tri_area(tris[it]);
                            for (int kk = 0; kk < dim; ++kk)
                                fwall[kk] += tau_wall[kk]*area;
                        }
                        
                        double* force_wall = intfc_state->shear_force;

                        for (int kk = 0; kk < dim; ++kk)
                        {
                            force_wall[kk] = fwall[kk];
                            f_surf[kk][ic] += -1.0*force_wall[kk]/rho;
                        }

                        //TODO: For elastic boundaries force_wall should be coupled
                        //      to the interface/canopy somehow... compute acceleration with it?
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: convective outlet boundary similar
                        //      to iF_flowThroughBoundaryState() needed?
                        if (boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            K_nb = K0;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                        }
                        else
                        {
                            //INLET

                            ////K_nb = eqn_params->Cbc
                            // * (sqr(intfc_state->vel[0])
                            // +  sqr(intfc_state->vel[1])
                            // +  sqr(intfc_state->vel[2]));
                            
                            K_nb = eqn_params->k_inlet;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                        }
                    }
                    else
                    {
                        printf("Unknown boundary condition %d! \n",wave_type(hs));
                        clean_up(ERROR);
                    }
                }
            }

            solver.Add_A(I,I,coeff);
            solver.Add_b(I,rhs);
        }
    
        break;
    }

    stop_clock("set_coefficients");
    
    solver.SetMaxIter(10000);
    solver.SetTol(1e-8);
    
    start_clock("petsc_solve");
    solver.Solve();
    stop_clock("petsc_solve");

	solver.GetNumIterations(&num_iter);
    solver.GetResidualNorm(&residual);

    if (debugging("PETSc"))
    {
        (void) printf("KE_CARTESIAN::computeAdvectionK: "
                    "num_iter = %d, residual = %g \n",
                    num_iter, residual);
    }

	if (residual > 1)
	{
	    printf("KE_CARTESIAN::computeAdvectionK(): \
                Solution diverges!\n");
	    LOC(); clean_up(ERROR);
	}

    //TODO: not parallelized? size = iupper - ilower ??
    FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
    solver.Get_x(x);

    start_clock("scatter_data");
    switch (dim)
    {
    case 1:
        for (i = imin; i <= imax; i++)
        {
            I = i_to_I[i];
            ic = d_index1d(i,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    case 2:
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            ic = d_index2d(i,j,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    case 3:
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            ic = d_index3d(i,j,k,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    }

    scatMeshArray();
    
    switch (dim)
    {
    case 1:
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index1d(i,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                K[ic] = array[ic];
        }
        break;
    case 2:
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index2d(i,j,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                K[ic] = array[ic];
        }
        break;
    case 3:
        for (k = 0; k <= top_gmax[2]; ++k)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index3d(i,j,k,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                K[ic] = array[ic];
        }
        break;
    }
    stop_clock("scatter_data");
    FT_FreeThese(1,x);

    FT_ParallelExchGridVectorArrayBuffer(f_surf,front);

    if (debugging("trace")) printf("Leaving computeAdvectionK()\n");
    stop_clock("computeAdvectionK");
}       /* end computeAdvectionK */

double KE_CARTESIAN::computePointFieldC2_RNG(int* icoords)
{
	int index = d_index(icoords,top_gmax,dim);

	double mu_t = field->mu_t[index];
	double K0 = field->k[index];
	double E0 = field->eps[index];
    
    double Cmu = eqn_params->Cmu;
    double C2 = eqn_params->C2;

	double S = computePointFieldStrain(icoords);
	
	double r = 0.0;
    if (mu_t == 0.0 || E0 == 0.0)
	    r = 0.0;
	else
	    r = S*K0/E0;
	
    double C2_star = C2 + (Cmu*r*r*r*(1.0 - r/4.38))/(1.0 + 0.012*r*r*r);

    return C2_star;
}

double KE_CARTESIAN::computePointFieldC1_REAL(int* icoords,double S)
{
	int index = d_index(icoords,top_gmax,dim);
	double K0 = field->k[index];
    double E0 = field->eps[index];

	if (E0 < 10000*MACH_EPS)
	{
	    return 1.0;
	}
	else
	{
	    double r = S*K0/E0;
	    r = std::max(r,0.0);
	    return std::max(0.43,r/(r + 5.0));
	}
}

void KE_CARTESIAN::computeAdvectionE_RNG(COMPONENT sub_comp)
{
	return computeAdvectionE_STD(sub_comp);
}

void KE_CARTESIAN::computeAdvectionE_REAL(COMPONENT sub_comp)
{
	return  computeAdvectionE_STD(sub_comp);
}       /* end computeAdvectionE */

void KE_CARTESIAN::computeAdvectionE_STD(COMPONENT sub_comp)
{
    int i,j,k,l,ll,m,ic,icn,I,I_nb,icoords[MAXD];
    int gmin[MAXD],ipn[MAXD];
    double crx_coords[MAXD];
    double K0,E0,E_nb,S,D,lambda,coeff,coeff_nb,rhs,C2;
    COMPONENT comp;
    PETSc solver;
    double *x;
    int num_iter = 0;
    double residual = 0.0;
    
    boolean fr_crx_grid_seg;
    const GRID_DIRECTION dir[3][2] =
            {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    
    double *E = field->eps;
	double *Pk = field->Pk;
	double *mu_t = field->mu_t;
    double *K_prev = field->k_prev;
	double *gamma = field->gamma;

    double Cmu = eqn_params->Cmu;
	double delta_eps = eqn_params->delta_eps;
	double rho = eqn_params->rho;
    double v[MAXD],v_wall[MAXD];
    double eta;
	double Ut; /*friction velocity*/
	double nu = eqn_params->mu/eqn_params->rho;
	
    /*For distance y*/
	double y;
	double center[MAXD], t[MAXD], point[MAXD],crds_wall[MAXD];
	double dist;
	
    /*For boundary state*/
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *intfc_state;
	INTERFACE* grid_intfc = front->grid_intfc;
	double coords[MAXD];
	double nor[MAXD];
	boolean if_adj_pt = NO;

    start_clock("computeAdvectionE");
    if (debugging("trace")) printf("Entering computeAdvectionE()\n");

    for (i = 0; i < dim; ++i) gmin[i] = 0;

    setIndexMap(sub_comp);
    if (debugging("trace"))
    {
        int domain_size = 1;
        printf("ilower = %d  iupper = %d\n",ilower,iupper);
        for (i = 0; i < dim; ++i)
            domain_size *= (imax-imin+1);
        printf("domain_size = %d\n",domain_size);
    }

    start_clock("set_coefficients");
	double dbg_max_rhs = 0, dbg_max_aii = 0;

    switch (dim)
    {
    case 2:

        solver.Create(ilower, iupper-1, 5, 5);

        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            ic = d_index2d(i,j,top_gmax);
            comp = top_comp[ic];
            I = ij_to_I[i][j];

            if (comp != sub_comp) continue;
            
            E0 = E[ic];
            K0 = K_prev[ic]; //K0 = field->k[ic];

            //computed in computeAdvectionK()
            double Gamma = gamma[ic];

            if (keps_model == REALIZABLE)
            {
                S = computePointFieldStrain(icoords);
                rhs = E0 + m_dt*std::max(computePointFieldC1_REAL(icoords,S)*S*E0,0.0);
            }
            else
            {
                rhs = E0 + m_dt*std::max(Pk[ic]*eqn_params->C1*Gamma,0.0); 
            }


            if (keps_model == RNG)
            {
                C2 = computePointFieldC2_RNG(icoords);
                coeff = 1.0 + std::max(C2*Gamma,0.0)*m_dt;
            }
            else if (keps_model == REALIZABLE)
            {
                C2 = eqn_params->C2;
                coeff = 1.0 + std::max(C2*E0/(K0+sqrt(std::max(nu*E0,0.0))),0.0)*m_dt;
            }
            else
            { 
                C2 = eqn_params->C2;
                coeff = 1.0 + std::max(C2*Gamma,0.0)*m_dt;
            }

            
            for (l = 0; l < dim; ++l) v[l] = 0.0;

            if (field->vel != NULL)
            {
                for (l = 0; l < dim; ++l)
                    v[l] = field->vel[l][ic];
            }

            /*set values at points adjacent to wall interface*/
            /*skip these points after settings*/

            //TODO: Disabled for now -- figure out what they were trying to do...    
            
            if_adj_pt = NO;
            
            /*
            for (l = 0; l < dim; ++l)
            for (m = 0; m < 2; ++m)
            {
                fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                    grid_intfc,icoords,dir[l][m],comp,
                                    (POINTER*)&intfc_state,&hs,crx_coords);
                
                if (fr_crx_grid_seg &&
                        (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                         wave_type(hs) == ELASTIC_BOUNDARY))
                {
                    if_adj_pt = YES;
                    for (ll = 0; ll < dim; ll++)
                        crds_wall[ll] = crx_coords[ll];
                }
            }
            */
            
            //disabled
            if (if_adj_pt == YES)
            {
                /*found adjacent point, use wall function*/
                getRectangleCenter(ic,center);
                dist = distance_between_positions(center,crds_wall,dim);
                
                /*set a lower bound for dist, since y+ > 11.067*/
                if (field->k[ic] > 0.0)
                {
                    dist = std::max(nu*eqn_params->y_p/(pow(eqn_params->Cmu,
                                    0.25)*pow(field->k[ic],0.25)),dist);
                    rhs = pow(eqn_params->Cmu,0.75)*pow(field->k[ic],1.5)/(0.41*dist);
                }
                else
                {
                    rhs = 0.0;
                }

                coeff = 1.0;
            }
            else
            {
                D = nu + mu_t[ic]/eqn_params->delta_eps/rho;

                for (l = 0; l < dim; ++l)
                {
                    lambda = D*m_dt/sqr(top_h[l]);
                    eta = v[l]*m_dt/(top_h[l]); //upwind difference
                    double eta_p = std::max(eta, 0.0);
                    double eta_m = std::min(eta, 0.0);
                    
                    coeff += eta_p - eta_m;
                    coeff += 2.0*lambda;

                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
                        I_nb = ij_to_I[ipn[0]][ipn[1]];
            
                        fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icoords,dir[l][m],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords);
                    
                        if (!fr_crx_grid_seg) 
                        {
                            coeff_nb = -lambda;
                            coeff_nb += (m == 0) ? -eta_p : eta_m;
                                //coeff_nb = -lambda + pow(-1,m+1)*eta;
                            solver.Add_A(I,I_nb,coeff_nb);
                        }
                        else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                                 wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                                 wave_type(hs) == ELASTIC_BOUNDARY)
                        {
                            double v_slip[MAXD] = {0.0};
                            setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_slip);

                            //use wall function 
                            double u_t =
                                std::max(pow(Cmu,0.25)*sqrt(std::max(K_prev[ic],0.0)),
                                        Mag2d(v_slip)/eqn_params->y_p);

                            //With K_prev so that friction velocity computation
                            //matches computeAdvectionK() friction velocity.
                            
                            /*
                            //use wall function
                            double u_t = std::max(pow(eqn_params->Cmu, 0.25)
                                *sqrt(std::max(field->k[ic], 0.0)), 
                                Mag2d(intfc_state->vel)/eqn_params->y_p);
                            */

                            //TODO: Note that intfc_state->vel is not
                            //      vel tangential to the wall
                            
                            E_nb = pow(u_t, 4.0)/(0.41*eqn_params->y_p*nu);
                            rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                        
                            //TODO: compute wall tangential stress (force/area)
                            //
                            //tau_wall[0] = -u_t/eqn_params->y_p*u_tangential[0];
                            //tau_wall[1] = -u_t/eqn_params->y_p*u_tangential[1];

                            //TODO: store tau_wall and apply in momentum eqns with f_surf[][] array??
                            //      would need to find area of nearest wall triangle, or bond length
                            //      in 2d case.
                        }
                        else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                        {
                            //TODO: convective outlet boundary similar
                            //      to iF_flowThroughBoundaryState() needed?
                            if (boundary_state_function_name(hs) &&
                                    strcmp(boundary_state_function_name(hs),
                                        "flowThroughBoundaryState") == 0)
                            {
                                //OUTLET
                                E_nb = E0;
                                rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                                    //rhs += lambda*E_nb + eta*pow(-1,m)*E_nb;
                            }
                            else
                            {
                                //INLET
                                
                                /*E_nb = eqn_params->Cmu
                                    *pow(eqn_params->Cbc,1.5)
                                    *pow(sqr(intfc_state->vel[0])
                                         + sqr(intfc_state->vel[1]),3.0)
                                    /eqn_params->l0;*/
                                
                                E_nb = eqn_params->eps_inlet;
                                rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                                    //rhs += lambda*E_nb + eta*pow(-1,m)*E_nb;
                            }
                        }
                        else
                        {
                            printf("Unknown boundary condition! \n");
                            clean_up(ERROR);
                        }
 
                    }  /*m*/

                }  /*l*/
     
            }  /*if_adj_pt*/

            if (isnan(coeff) || isinf(coeff) || isnan(rhs) || isinf(rhs))
            {
                printf("Warning: at [%d %d], coeff = %f, rhs = %f\n",i,j,coeff,rhs);
                LOC(); clean_up(ERROR);
            }

            solver.Add_A(I,I,coeff);
            solver.Add_b(I,rhs);
    
            dbg_max_aii = std::max(coeff, dbg_max_aii);
            dbg_max_rhs = std::max(rhs, dbg_max_rhs);
        }
        break;

    case 3:

        solver.Create(ilower, iupper-1, 7, 7);
        
        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            ic = d_index3d(i,j,k,top_gmax);
            comp = top_comp[ic];
            I = ijk_to_I[i][j][k];
            
            if (comp != sub_comp) continue;

            E0 = E[ic];
            K0 = K_prev[ic]; //K0 = field->k[ic];
    
            //computed in computeAdvectionK()
            double Gamma = gamma[ic];

            if (keps_model == REALIZABLE)
            {
                S = computePointFieldStrain(icoords);
                rhs = E0 + m_dt*computePointFieldC1_REAL(icoords,S)*S*E0;
            }
            else
            {
                rhs = E0 + m_dt*std::max(Pk[ic]*eqn_params->C1*Gamma,0.0); 
            }

            if (keps_model == RNG)
            {
                C2 = computePointFieldC2_RNG(icoords);
                coeff = 1.0 + std::max(C2*Gamma,0.0)*m_dt;
            }
            else if (keps_model == REALIZABLE)
            {
                C2 = eqn_params->C2;
                coeff = 1.0 + std::max(C2*E0/(K0+sqrt(std::max(nu*E0,0.0))),0.0)*m_dt;
            }
            else
            { 
                C2 = eqn_params->C2;
                coeff = 1.0 + std::max(C2*Gamma,0.0)*m_dt;
            }
            
            
            for (l = 0; l < dim; ++l) v[l] = 0.0;
            
            if (field->vel != NULL)
            {
                for (l = 0; l < dim; ++l)
                    v[l] = field->vel[l][ic];
            }
        
            D = nu + mu_t[ic]/eqn_params->delta_eps/rho;

            for (l = 0; l < dim; ++l)
            {
                lambda = D*m_dt/sqr(top_h[l]);
                eta = v[l]*m_dt/(top_h[l]); //upwind difference
                double eta_p = std::max(eta, 0.0);
                double eta_m = std::min(eta, 0.0);
                
                coeff += eta_p - eta_m;
                coeff += 2.0*lambda;

                for (m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                    I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
            
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                            grid_intfc,icoords,dir[l][m],comp,
                            (POINTER*)&intfc_state,&hs,crx_coords);
            
                    if (!fr_crx_grid_seg) 
                    {
                        coeff_nb = -lambda;
                        coeff_nb += (m == 0) ? -eta_p : eta_m;
                        solver.Add_A(I,I_nb,coeff_nb);
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                             wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        double v_slip[MAXD] = {0.0};
                        setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_slip);

                        //use wall function 
                        double u_t =
                            std::max(pow(Cmu,0.25)*sqrt(std::max(K_prev[ic],0.0)),
                                    Mag3d(v_slip)/eqn_params->y_p);
                        
                        //With K_prev so that friction velocity computation
                        //matches computeAdvectionK() friction velocity.

                        /*
                        double u_t = std::max(pow(eqn_params->Cmu, 0.25)
                            *sqrt(std::max(field->k[ic], 0.0)), 
                            Mag3d(intfc_state->vel)/eqn_params->y_p);
                        */

                        //TODO: Note that intfc_state->vel is not
                        //      vel tangential to the wall
                        
                        E_nb = pow(u_t, 4.0)/(0.41*eqn_params->y_p*nu);
                        rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: convective outlet boundary similar
                        //      to iF_flowThroughBoundaryState() needed?
                        if (boundary_state_function_name(hs) &&
                        strcmp(boundary_state_function_name(hs),
                                        "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            E_nb = E0;
                            rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                        }
                        else
                        {
                            //INLET
                            
                            /*E_nb = eqn_params->Cmu
                             *pow(eqn_params->Cbc,1.5)
                             *pow(Mag3d(intfc_state->vel), 3.0)
                             /eqn_params->l0;*/

                            E_nb = eqn_params->eps_inlet;
                            rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                        }
                    }
                    else
                    {
                        printf("Unknown boundary condition! \n");
                        clean_up(ERROR);
                    }
                }  /*m*/

            }  /*l*/
     
            solver.Add_A(I,I,coeff);
            solver.Add_b(I,rhs);
    
            dbg_max_aii = std::max(coeff, dbg_max_aii);
            dbg_max_rhs = std::max(rhs, dbg_max_rhs);
        }
        
        break;
    
    }//end switch (dim)
        
    stop_clock("set_coefficients");
    
    solver.SetMaxIter(10000);
    solver.SetTol(1e-8);
    
    start_clock("petsc_solve");
    solver.Solve();
    stop_clock("petsc_solve");
	
    solver.GetNumIterations(&num_iter);
    solver.GetResidualNorm(&residual);
    
    if (debugging("PETSc"))
    {
        (void) printf("KE_CARTESIAN::computeAdvectionE: "
                    "num_iter = %d, residual = %g \n",
                    num_iter, residual);
        printf("max rhs = %f\n", dbg_max_rhs);
        printf("max aii = %f\n", dbg_max_aii);
    }

	if (residual > 1)
	{
	    printf("KE_CARTESIAN::computeAdvectionE(): \
                Solution diverges!\n");
	    LOC(); clean_up(ERROR);
	}

    //TODO: not parallelized? size = iupper - ilower ??
    FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
    solver.Get_x(x);

    start_clock("scatter_data");
    switch (dim)
    {
    case 1:
        for (i = imin; i <= imax; i++)
        {
            I = i_to_I[i];
            ic = d_index1d(i,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    case 2:
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            ic = d_index2d(i,j,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    case 3:
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            ic = d_index3d(i,j,k,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                array[ic] = x[I-ilower];
            else
                array[ic] = 0.0;
        }
        break;
    }
    
    scatMeshArray();
    
    switch (dim)
    {
    case 1:
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index1d(i,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                E[ic] = array[ic];
        }
        break;
    case 2:
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index2d(i,j,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                E[ic] = array[ic];
        }
        break;
    case 3:
        for (k = 0; k <= top_gmax[2]; ++k)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index3d(i,j,k,top_gmax);
            comp = cell_center[ic].comp;
            if (comp == sub_comp)
                E[ic] = array[ic];
        }
        break;
    }
    stop_clock("scatter_data");
    FT_FreeThese(1,x);

    if (debugging("trace")) printf("Leaving computeAdvectionE()\n");
    stop_clock("computeAdvectionE");
}       /* end computeAdvectionE */

void KE_CARTESIAN::solve(double dt)
{
	if (debugging("trace")) printf("Entering solve()\n");
	m_dt = dt;
	start_clock("solve");

    if (!activated)
    {
	    if (front->time <= eqn_params->t0) return;

        applyInitialConditions();
        activateKE();

        printf("\n\nTurbulence Model Activated\n\n");

        //TODO: needs to be called every time when interface is allowed to move
            //computeDistances();
    }

    setDomain();
    if (debugging("trace")) printf("Passing setDomain()\n");

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");

    computeSource();
    if (debugging("trace")) printf("Passing computeSource()\n");

	computeAdvection();
	if (debugging("trace")) printf("Passing computeAdvection()\n");

	computeMuTurb();
	if (debugging("trace")) printf("Passing computeMuTurb()\n");

        //setAdvectionDt();
        //if (debugging("trace")) printf("Passing setAdvectionDt()\n");

	stop_clock("solve");
	if (debugging("trace")) printf("Leaving solve()\n");
}

static void printField2d(double *var,
		       const char* varname, 
		       int* lmin, 
		       int* lmax,
		       int* top_gmax)
{
	int index;
	FILE* outfile;
	outfile = fopen(varname,"w");

	for (int j = lmin[1]; j <= lmax[1]; ++j)
	{
	    for (int i = lmin[0]; i <= lmax[0]; ++i)
	    {
	        index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%e ",var[index]);
	    }
	    fprintf(outfile,"\n");
	}
	fclose(outfile);	
}

static void printField3d(double *var,
		       const char* varname, 
		       int* lmin, 
		       int* lmax,
		       int* top_gmax)
{
	int index;
	FILE* outfile;
	outfile = fopen(varname,"w");
	
    int i = (lmin[0] + lmax[0])/2;
	for (int k = lmin[2]; k <= lmax[2]; ++k)
	{
	    for (int j = lmin[1]; j <= lmax[1]; ++j)
	    {
	        index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%e ",var[index]);
	    }
	    fprintf(outfile,"\n");
	}
	fclose(outfile);	
}

//S = sqrt(2*Sij*Sij)
//Pk*rho/mu_t = 0.5*(Uij + Uji)^2 = 2*Sij*Sij
//Sij = 0.5*(Uij + Uji) and Uij = dUi/dxj
double KE_CARTESIAN::computePointFieldStrain(int* icoords)
{
	int index = d_index(icoords,top_gmax,dim);

    if (fabs(field->mu_t[index]) > MACH_EPS)
        return sqrt(field->Pk[index]*eqn_params->rho/field->mu_t[index]);
    else
        return 0.0;
}

double KE_CARTESIAN::computePointFieldCmu(int* icoords)
{
	int i,j,k,l,m,index;
	char fname[200];
	static int count = 0;
	int lmin[MAXD], lmax[MAXD];
	INTERFACE *grid_intfc = front->grid_intfc;
	
	boolean fr_crx_grid_seg;
    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

    STATE* intfc_state;
	COMPONENT comp;
	double crx_coords[MAXD];
	boolean Adj_Pt = NO;
	double Ut, y, dist;
	HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
	double center[MAXD], t[MAXD], point[MAXD],nor[MAXD];
	double rho = eqn_params->rho;
	double nu = eqn_params->mu/rho;
	double S[MAXD][MAXD],R[MAXD][MAXD],Cmu,J,U,W,phi,A0,As;
	double **vel = field->vel;
	double *K = field->k;
	double *E = field->eps;
	double d_h[2],vel_nb[2],v_tmp[MAXD];
	int index_nb,nb;

	index = d_index(icoords,top_gmax,dim);
  	comp = top_comp[index];
	getRectangleCenter(index,center);

 	if (!ifluid_comp(comp)) return 0.0;
	
    /*compute module of the strain rate tensor*/
	for (l = 0; l < dim; l++)
	for (m = 0; m < dim; m++)
	{
        //l derivatives in m direction
	    for (nb = 0; nb < 2; nb++)
	    {
            fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    grid_intfc,icoords,dir[m][nb],comp,
                    (POINTER*)&intfc_state,&hs,crx_coords);
            
            d_h[nb] = top_h[m]; 
		
            if (!fr_crx_grid_seg)
            {
                index_nb = next_index_in_dir(icoords,dir[m][nb],top_gmax,dim);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                     wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                     wave_type(hs) == ELASTIC_BOUNDARY)
            {
                setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,field->vel,v_tmp);
                vel_nb[nb] = v_tmp[l];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                if (boundary_state_function_name(hs) &&
                        strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                {
                    //OUTLET
                    vel_nb[nb] = vel[l][index];
                }
                else
                {
                    //INLET
                    vel_nb[nb] = intfc_state->vel[l];
                }
            }
        }

        S[l][m] = 0.5*(vel_nb[1]- vel_nb[0])/(d_h[1]+d_h[0]);
	    R[l][m] = S[l][m];

        //m derivatives in l direction
	    for (nb = 0; nb < 2; nb++)
        {
            fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    grid_intfc,icoords,dir[l][nb],comp,
                    (POINTER*)&intfc_state,&hs,crx_coords);
            
            d_h[nb] = top_h[l];
		
            if (!fr_crx_grid_seg)
            {
                index_nb = next_index_in_dir(icoords,dir[l][nb],top_gmax,dim);
                vel_nb[nb] = vel[m][index_nb];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                     wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                     wave_type(hs) == ELASTIC_BOUNDARY)
            {
                setSlipBoundary(icoords,l,nb,comp,hs,intfc_state,field->vel,v_tmp);
                vel_nb[nb] = v_tmp[m];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                if (boundary_state_function_name(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                {
                    //OUTLET
                    vel_nb[nb] = vel[m][index];
                }
                else
                {
                    //INLET
                    vel_nb[nb] = intfc_state->vel[m];
                }
            }
	    }

        S[l][m] += 0.5*(vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
	    R[l][m] -= 0.5*(vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);

	    if (isnan(S[l][m]))
        {
            printf("wave_type = %d, vel_nb = [%f %f], d_h = [%f %f]\n",
                    wave_type(hs),vel_nb[0],vel_nb[1],d_h[0],d_h[1]);
        }
	}

	for (l = 0; l < dim; l++)
	for (m = 0; m < dim; m++)
	{
	    S[l][l] -= 1.0/(double)dim*S[m][m];   
	}

	J = 0.0;
	for (l = 0; l < dim; l++)
	for (m = 0; m < dim; m++)
	{
	    J += S[m][l]*S[m][l];
	}

	J = sqrt(J);
	
    if (J != 0.0)
	{
	    W = 1.0;
	    for (l = 0; l < dim; l++)
	    for (m = l+1; m < dim; m++)
	    {
	        W *= S[l][m]/J;   
	    }
	}
	else
    {
	    W = 0.0;
    }

	W = W > 0 ? std::min(W,1.0/sqrt(6.0)) : std::max(W,-1.0/sqrt(6.0));
	phi = sqrt(6.0)*W;
	phi = 1.0/3.0*acos(phi);
	A0 = 4.04; As = sqrt(6.0)*cos(phi);

	U = 0.0;
	for (l = 0; l < dim; l++)
	for (m = 0; m < dim; m++)
	{
	    U += S[l][m]*S[l][m]+R[l][m]*R[l][m];
	}

    U = sqrt(U);

	if (E[index] == 0)
	    Cmu = 0.0;
	else
	    Cmu = 1.0/(A0+As*U*std::max(K[index]/E[index],0.0));
	Cmu = std::min(Cmu,eqn_params->Cmu);
	
	if (isnan(Cmu) || isinf(Cmu))
	{
	    printf("Warning: Cmu = %f, E = %f, K = %f, A0 = %f, As = %f, U = %f\n",
		    Cmu,E[index],K[index],A0,As,U);
	    printf("phi = %f, W = %f\n",phi,W);
	    LOC(); clean_up(ERROR);
	}

	return Cmu;
}

void KE_CARTESIAN::computeMuTurb()
{
	int i,j,k,l,m,ll,index;
	char fname[200];
	static int count = 0;
	int lmin[MAXD], lmax[MAXD];
	INTERFACE *grid_intfc = front->grid_intfc;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
        POINTER intfc_state;
	COMPONENT comp;
	int icoords[MAXD];
	double crx_coords[MAXD];
	boolean Adj_Pt = NO;
	boolean fr_crx_grid_seg;
	double Ut, y, dist,vn;
	HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double center[MAXD], t[MAXD], point[MAXD],nor[MAXD],v_wall[MAXD];
	double rho = eqn_params->rho;
	double nu = eqn_params->mu/rho;
	double Cmu = eqn_params->Cmu;

	lmin[0] = imin; lmin[1] = jmin; lmin[2] = kmin;
	lmax[0] = imax; lmax[1] = jmax; lmax[2] = kmax;

	switch(dim)
	{
	    case 2:
		for (i = imin; i <= imax; i++)
		for (j = jmin; j <= jmax; j++)
		{
            icoords[0] = i;
            icoords[1] = j;
            index = d_index2d(i,j,top_gmax);
            comp = top_comp[index];
        
            if (keps_model == REALIZABLE)
		    {
			    Cmu = computePointFieldCmu(icoords);
			    field->Cmu[index] = Cmu;
		    }
		    else
            {
                Cmu = eqn_params->Cmu;
            }

		    if (field->eps[index] != 0.0)
            {
                field->mu_t[index] = Cmu*sqr(field->k[index])/field->eps[index]*eqn_params->rho;
            }
		    else
            {
                field->mu_t[index] = 0.0;
                //field->mu_t[index] = 0.0001*eqn_params->mu;
            }

            if (isnan(field->mu_t[index]) || isinf(field->mu_t[index]))
		    {
                printf("Warning: mu_t=%f,Cmu=%f, k=%f, eps=%f\n",
                        field->mu_t[index],Cmu,field->k[index],field->eps[index]);
                LOC(); clean_up(EXIT_FAILURE);
                //field->mu_t[index] = 0.0001*eqn_params->mu;
		    }

            //field->mu_t[index] = std::max(field->mu_t[index],0.0001*eqn_params->mu);
		}
		break;

        case 3:
		for (k = kmin; k <= kmax; k++)
		for (j = jmin; j <= jmax; j++)
		for (i = imin; i <= imax; i++)
		{
            icoords[0] = i;
            icoords[1] = j;
		    icoords[2] = k;
            
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];

		    if (keps_model == REALIZABLE)
		    {
                Cmu = computePointFieldCmu(icoords);
                field->Cmu[index] = Cmu;
		    }
		    else
            {
                Cmu = eqn_params->Cmu;
            }

		    if (field->eps[index] != 0.0)
            {
                field->mu_t[index] = Cmu*sqr(field->k[index])/field->eps[index]*eqn_params->rho;
            }
            else
            {
                field->mu_t[index] = 0.0;
                //field->mu_t[index] = std::max(field->mu_t[index],0.0001*eqn_params->mu);
            }

            if (isnan(field->mu_t[index]) || isinf(field->mu_t[index]))
		    {
                printf("Warning: mu_t=%f,Cmu=%f, k=%f, eps=%f\n",
                        field->mu_t[index],Cmu,field->k[index],field->eps[index]);
                LOC(); clean_up(EXIT_FAILURE);
            }

            //TODO: Limit far field eddy viscosity with ViscRatioFar.
            //      Need to determine if "far" enough to impose freestream
            //      boundary condition.
            //
            /*status = FT_FindNearestIntfcPointInRange(front,comp,coords,
                    NO_BOUNDARIES,point,t,&hse,&hs,range);*/
            //if (!status) .... mu_t[index] = ViscRatioFar*mu;
		}
        break;
	}

	FT_ParallelExchGridArrayBuffer(field->k,front,NULL);
	FT_ParallelExchGridArrayBuffer(field->eps,front,NULL);
	FT_ParallelExchGridArrayBuffer(field->mu_t,front,NULL);
	
    if (keps_model == REALIZABLE)
	    FT_ParallelExchGridArrayBuffer(field->Cmu,front,NULL);

	if (dim == 2)
	{
	    sprintf(fname,"%s/K_field",OutName(front));
	    printField2d(field->k,fname,lmin,lmax,top_gmax);
	    sprintf(fname,"%s/E_field",OutName(front));
	    printField2d(field->eps,fname,lmin,lmax,top_gmax);
	}
	else if (dim == 3)
	{
            sprintf(fname,"%s/K_field",OutName(front));
            printField3d(field->k,fname,lmin,lmax,top_gmax);
            sprintf(fname,"%s/E_field",OutName(front));
            printField3d(field->eps,fname,lmin,lmax,top_gmax);
	    if (keps_model == REALIZABLE)
	    {
		sprintf(fname,"%s/Cmu",OutName(front));
		printField3d(field->Cmu,fname,lmin,lmax,top_gmax);
	    }
	}

}

//TODO: update y+ here?
void KE_CARTESIAN::computeDistances()
{
    double* dist_array = field->dist;
    double** nor_array = field->nor;

    int index;
    COMPONENT comp;
    double coords[MAXD];
    double intfc_coords[MAXD];
    double nor[MAXD];
    bool status;

    std::string fname = OutName(front);
    fname +=  "/dist-nor.txt";
    FILE* file = fopen(fname.c_str(),"w");

	switch (dim)
	{
	    case 2:
		
            for (int i = imin; i <= imax; ++i)
            for (int j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                comp = cell_center[index].comp;

                if (!ifluid_comp(comp))
                {
                    dist_array[index] = -1;
                    for (int l = 0; l < dim; ++l)
                        nor_array[l][index] = 0;
                    fprintf(file,"index = %d coords = %f %f \n\t \
                            dist = %g    nor = %f %f \n", index,coords[0],coords[1],
                            dist_array[index],nor_array[0][index],nor_array[1][index]);
                    continue;
                }

                getRectangleCenter(index,coords);
                status = getNearestInterfacePoint(comp,coords,intfc_coords,nor);
                
                if (status)
                    dist_array[index] = distance_between_positions(coords,intfc_coords,dim);
                else
                    dist_array[index] = -2;
                
                for (int l = 0; l < dim; ++l)
                    nor_array[l][index] = nor[l];

                fprintf(file,"index = %d coords = %f %f \n\t \
                        dist = %g    nor = %f %f \n", index,coords[0],coords[1],
                        dist_array[index],nor_array[0][index],nor_array[1][index]);
            }

	    case 3:
		
            for (int i = imin; i <= imax; ++i)
            for (int j = jmin; j <= jmax; ++j)
            for (int k = kmin; k <= kmax; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                
                comp = cell_center[index].comp;
                if (!ifluid_comp(comp))
                {
                    dist_array[index] = -1;
                    for (int l = 0; l < dim; ++l)
                        nor_array[l][index] = 0;
                    fprintf(file,"index = %d coords = %f %f %f \n\t \
                            dist = %g    nor = %f %f %f\n", index,coords[0],coords[1],coords[2],
                            dist_array[index],nor_array[0][index],nor_array[1][index],nor_array[2][index]);
                    continue;
                }

                getRectangleCenter(index,coords);
                status = getNearestInterfacePoint(comp,coords,intfc_coords,nor);

                if (status)
                    dist_array[index] = distance_between_positions(coords,intfc_coords,dim);
                else
                    dist_array[index] = -2;

                for (int l = 0; l < dim; ++l)
                    nor_array[l][index] = nor[l];

                fprintf(file,"index = %d coords = %f %f %f \n\t \
                        dist = %g    nor = %f %f %f\n", index,coords[0],coords[1],coords[2],
                        dist_array[index],nor_array[0][index],nor_array[1][index],nor_array[2][index]);
            }

    }
    fclose(file);

    FT_ParallelExchGridArrayBuffer(dist_array,front,NULL);
    FT_ParallelExchGridVectorArrayBuffer(nor_array,front);
}

bool KE_CARTESIAN::getNearestInterfacePoint(
        COMPONENT comp,     // mesh point component
        double *p,          // mesh point coords
	    double *q,          // intfc point coords
        double* nor)        // intfc normal at point
        //double* kappa)      // intfc curvature at point
{
	INTERFACE *intfc = front->interf;
	int dim = front->rect_grid->dim;
	double t[MAXD],mag_nor;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	POINT *pts[MAXD];
	BOND *bond;
	TRI *tri;
	

    static double **pts_nor, *pts_kappa;
 
    if (pts_kappa == NULL)
    {
        FT_VectorMemoryAlloc((POINTER*)&pts_kappa,MAXD,FLOAT);
        FT_MatrixMemoryAlloc((POINTER*)&pts_nor,MAXD,MAXD,FLOAT);
    }

    //TODO: use nearest_similar_interface_point3d() instead?
	

    boolean status;

    status = nearest_interface_point(p,comp,intfc,NO_BOUNDARIES,NULL,q,t,&hse,&hs);
        
        //int range = 10;
        //status = FT_FindNearestIntfcPointInRange(front,comp,p,NO_BOUNDARIES,q,t,&hse,&hs,range);
    
    if (hse != NULL)
	{
	    switch (dim)
        {
            case 2:

            bond = Bond_of_hse(hse);
            pts[0] = bond->start;
            pts[1] = bond->end;
            for (int i = 0; i < dim; ++i)
            {
                //GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
                GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
            }
            //*kappa = (1.0 - t[0])*pts_kappa[0] + t[0]*pts_kappa[1];
            for (int i = 0; i < dim; ++i)
                nor[i] = (1.0 - t[0])*pts_nor[0][i] + t[0]*pts_nor[1][i];
            mag_nor = mag_vector(nor,dim);
            for (int i = 0; i < dim; ++i)
                nor[i] /= mag_nor;
            break;

            
            case 3:

            tri = Tri_of_hse(hse);	
            for (int i = 0; i < dim; ++i)
            {
                pts[i] = Point_of_tri(tri)[i];
                //GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
                GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
            }
            for (int i = 0; i < dim; ++i)
            {
                //*kappa = t[i]*pts_kappa[i];
                for (int j = 0; j < dim; ++j)
                nor[j] = t[i]*pts_nor[i][j];
            }
            mag_nor = mag_vector(nor,dim);
            for (int i = 0; i < dim; ++i)
                nor[i] /= mag_nor;
        }
	}
	else
	{
	    //*kappa = 0.0;
	    for (int i = 0; i < dim; ++i)
            nor[i] = 0.0;
	}


    if (status) return true;
    return false;
}	/* end getNearestInterfacePoint */


	
void KE_CARTESIAN::setAdvectionDt()
{
	int i;
	double mu_max = -HUGE;
	double D,Dl,Ds;

	for (i = 0; i < comp_size; i++)
		mu_max = std::max(field->mu_t[i],mu_max);
	Dl = mu_max/eqn_params->delta_k   + eqn_params->mu;
	Ds = mu_max/eqn_params->delta_eps + eqn_params->mu;
	D = std::max(Dl,Ds)/eqn_params->rho;
	m_dt = sqr(hmin)/D*Time_step_factor(front);

	front->dt = std::min(m_dt,front->dt);
	if (debugging("trace"))
	{
	    printf("In setAdvectionDt: m_dt = %24.18g min_dt = %f\n",
				m_dt,min_dt);
	}
}	/* end setAdvectionDt */

void KE_CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void KE_CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int KE_CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void KE_CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

void KE_CARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(cell_center[index0].coords[i] +
	    		     cell_center[index1].coords[i]);
	}
}

int KE_CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int KE_CARTESIAN::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 1:
	    index = d_index1d(icoords[0],top_gmax);
	    return top_comp[index];
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	}
}

void KE_CARTESIAN::save(char *filename)
{
	RECT_GRID *rect_grid = front->rect_grid;
	INTERFACE *intfc    = front->interf;
		
	int i, j;
	int xmax = rect_grid->gmax[0];
	int ymax = rect_grid->gmax[1];
	double x, y;
		
#if defined(HAVE_MPI)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(HAVE_MPI) */

	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		clean_up(EXIT_FAILURE);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
	    CURVE		**cur;
	    CURVE		*curve;
	    BOND		*bond;
			
	    for(cur=intfc->curves; cur && *cur; cur++)
	    {
		curve = *cur;
		fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", 
				curve->num_points, 1);
		bond=curve->first;
		fprintf(hfile, "%.4f %.4f \n",bond->start->_coords[0], 
				bond->start->_coords[1]);
		for(bond=curve->first; bond!=NULL; bond=bond->next)
		    fprintf(hfile, "%.4f %.4f \n",bond->end->_coords[0], 
		    		bond->end->_coords[1]);
		}					
	}		
	fclose(hfile);
}

void KE_CARTESIAN::deleteGridIntfc()
{
	FT_FreeGridIntfc(front);
}

void KE_CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
}

void KE_CARTESIAN::setGlobalIndex(COMPONENT sub_comp)
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(int));
	}
	NLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];

        //TODO: see comment at top of file
    for (i = 1; i <= myid; i++)
    {
        ilower += n_dist[i-1];
        iupper += n_dist[i];
    }	
}

void KE_CARTESIAN::printFrontInteriorState(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double *K = field->k;
	double *eps = field->eps;

	sprintf(filename,"%s/state.ts%s-keps",out_name,right_flush(front->step,7));
        
#if defined(HAVE_MPI)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(HAVE_MPI) */
	outfile = fopen(filename,"w");
	
        /* Initialize states at the interface */
        fprintf(outfile,"Interface states:\n");
        int count = 0;
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            count++;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",sl->temperature,
				sr->temperature);
        }

	fprintf(outfile,"\nInterior states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%24.18g %24.18g\n",K[index],eps[index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g %24.18g\n",K[index],eps[index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g %24.18g\n",K[index],eps[index]);
	    }
	    break;
	}
	fclose(outfile);
}

void KE_CARTESIAN::readFrontInteriorState(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double x;
	double *K = field->k;
	double *eps = field->eps;
	double *mu_t = field->mu_t;

	char fname[100];
	sprintf(fname,"%s-keps",restart_name);
#if defined(HAVE_MPI)
        if (pp_numnodes() > 1)
            sprintf(fname,"%s-nd%s",fname,right_flush(pp_mynode(),4));
#endif /* defined(HAVE_MPI) */

	infile = fopen(fname,"r");

        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            sl->temperature = x;
            fscanf(infile,"%lf",&x);
            sr->temperature = x;
        }

	FT_MakeGridIntfc(front);
        setDomain();

        /* Initialize states in the interior regions */

	next_output_line_containing_string(infile,"Interior states:");

	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	    	fscanf(infile,"%lf %lf",&K[index],&eps[index]);
		mu_t[index] = eqn_params->Cmu*sqr(K[index])/eps[index];
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf %lf",&K[index],&eps[index]);
		mu_t[index] = eqn_params->Cmu*sqr(K[index])/eps[index];
	    }
	    break;
	}
	fclose(infile);
}

void KE_CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i;

	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	if (field == NULL)
	    FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(IF_FIELD));
	switch (dim)
	{
	case 1:
            if (first)
            {
		comp_size = top_gmax[0]+1;
                FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->Pk,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->k,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->eps,comp_size,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field = field;
	    break;
	case 2:
	    if (first)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1);
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->Pk,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->k,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->eps,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->k_prev,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->mu_t,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->gamma,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->dist,comp_size,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&field->nor,2,comp_size,sizeof(double));
		if (keps_model == REALIZABLE)
		    FT_VectorMemoryAlloc((POINTER*)&field->Cmu,comp_size,FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field = field;
	    break;
	case 3:
	    if (first)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->Pk,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->k,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->eps,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->k_prev,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->mu_t,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->gamma,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->dist,comp_size,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&field->nor,3,comp_size,sizeof(double));
		if (keps_model == REALIZABLE)
		    FT_VectorMemoryAlloc((POINTER*)&field->Cmu,comp_size,FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field = field;
	    break;
	}
}	/* end setDomain */

void KE_CARTESIAN::initMovieVariables()
{
	switch (dim)
	{
	case 2:
	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,"mu_t",0,
				field->mu_t,getStateTemp,0,0);
	    break;
	case 3:
        /* Added for vtk movie of scalar field */
        FT_AddVtkScalarMovieVariable(front,"mu_t",field->mu_t);
	    break;
	}

        if (debugging("trace"))
            printf("Leaving initMovieVariables()\n");
}	/* end initMovieVariables */

static int find_state_at_crossing(
	Front *front,
	int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
	boolean status;
	INTERFACE *grid_intfc = front->grid_intfc;

	status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
				comp,state,hs,crx_coords);
        
    if (status == NO)
        return NO_PDE_BOUNDARY;

    if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
    {
	    return NO_PDE_BOUNDARY;
    }
	else if (wave_type(*hs) == DIRICHLET_BOUNDARY)
	{
        return DIRICHLET_PDE_BOUNDARY;
	}
	else if (wave_type(*hs) == GROWING_BODY_BOUNDARY)
	{
        return DIRICHLET_PDE_BOUNDARY;
	}
	else if (wave_type(*hs) == NEUMANN_BOUNDARY || 
            wave_type(*hs) == ELASTIC_BOUNDARY)
    {
        return NEUMANN_PDE_BOUNDARY;
    }
}       /* find_state_at_crossing */

//TODO: incomplete
//NOTE: only used for 2d
void KE_CARTESIAN::setTKEatWall(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double* K,
	double* K_nb)
{
	int             i,j,index;
        int             ic[MAXD];
        double          coords[MAXD],coords_ref[MAXD],crx_coords[MAXD];
        double          nor[MAXD],vn,v[MAXD];
        GRID_DIRECTION  ldir[3] = {WEST,SOUTH,LOWER};
        GRID_DIRECTION  rdir[3] = {EAST,NORTH,UPPER};
        GRID_DIRECTION  dir;
        double  vel_ref[MAXD];

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
        {
            coords[i] = top_L[i] + icoords[i]*top_h[i];
            ic[i] = icoords[i];
        }
	dir = (nb == 0) ? ldir[idir] : rdir[idir];
        FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);
	ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
        for (j = 0; j < dim; ++j)
            coords_ref[j] = top_L[j] + ic[j]*top_h[j];

        /* Reflect ghost point through intfc-mirror at crossing */
        coords_ref[idir] = 2.0*crx_coords[idir] - coords_ref[idir];
        vn = 0.0;
        for (j = 0; j < dim; ++j)
        {
            v[j] = coords_ref[j] - crx_coords[j];
            vn += v[j]*nor[j];
        }
        for (j = 0; j < dim; ++j)
            v[j] = 2.0*vn*nor[j] - v[j];
        for (j = 0; j < dim; ++j)
            coords_ref[j] = crx_coords[j] + v[j];

        /* Interpolate the state at the reflected point */
        /*for (j = 0; j < dim; ++j)
            FT_IntrpStateVarAtCoords(front,comp,coords_ref,K,
        	getStateK,K_nb,NULL);*/
	if (rect_in_which(coords_ref,ic,top_grid))
	{
	    index = d_index(ic,top_gmax,dim);
	}
	else
	{
	    printf("ERROR: point [%f %f] is out of domain\n",
		    coords_ref[0],coords_ref[1]);
	    LOC();
	    clean_up(ERROR);
	}
	*K_nb = K[index];
}

//TODO: see iFbasic.cpp implemtation -- should be changed, but saving for now,
//      since kepsilon solver development is on hold for the moment.
void KE_CARTESIAN::setSlipBoundary(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double** vel,
	double* v_tmp)
{
	int             i,j,index;
    int             ic[MAXD];
    double          coords[MAXD],coords_ref[MAXD],crx_coords[MAXD];
    double          nor[MAXD],vn,v[MAXD];
    
    GRID_DIRECTION  ldir[3] = {WEST,SOUTH,LOWER};
    GRID_DIRECTION  rdir[3] = {EAST,NORTH,UPPER};
    GRID_DIRECTION  dir;
    double  vel_intfc[MAXD];

	index = d_index(icoords,top_gmax,dim);
	
    for (i = 0; i < dim; ++i)
    {
        vel_intfc[i] = (*getStateVel[i])(state);
        coords[i] = top_L[i] + icoords[i]*top_h[i];
        ic[i] = icoords[i];
    }
	
    dir = (nb == 0) ? ldir[idir] : rdir[idir];

    FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);
	
    ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;

    //ghost coords
    for (j = 0; j < dim; ++j)
        coords_ref[j] = top_L[j] + ic[j]*top_h[j];

        
    /* Reflect ghost point through intfc-mirror at crossing */
    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ref[idir];
    vn = 0.0;
    for (j = 0; j < dim; ++j)
    {
        v[j] = coords_ref[j] - crx_coords[j];
        vn += v[j]*nor[j];
    }

    for (j = 0; j < dim; ++j)
        v[j] = 2.0*vn*nor[j] - v[j];

    for (j = 0; j < dim; ++j)
        coords_ref[j] = crx_coords[j] + v[j];

    /* Interpolate the state at the reflected point */
    for (j = 0; j < dim; ++j)
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,vel[j],
                getStateVel[j],&v_tmp[j],&vel[j][index]);

    double vel_rel[MAXD] = {0.0};
    vn = 0.0;

    for (j = 0; j < dim; ++j)
    {
        //Relative velocity of reflected point to boundary
        vel_rel[j] = v_tmp[j] - vel_intfc[j];
                //vel_rel[j] = vel_reflect[j] - vel_intfc[j];
        //Normal component of the relative velocity
            //vn += vel_rel[j]*nor[j];
    }

    /*normal component equal to zero while tangential component is permitted*/

    /*
    //double v_slip[MAXD] = {0.0};
    for (j = 0; j < dim; ++j)
	    v_tmp[j] -= vn*nor[j];
	    //v_slip[j] = v_tmp[j] - vn*nor[j];

    fprint_general_vector(stdout,"nor",nor,dim,"\n");
    fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
    */

    //TODO: Use normal vector instead??
    //      May actually  be equivalent to the vector v
    //      being computed.
	
    for (j = 0; j < dim; ++j)
        v[j] = coords_ref[j] - (top_L[j] + ic[j]*top_h[j]);

    for (j = 0; j < dim; ++j)
        v[j] /= mag_vector(v,dim);

    vn = 0.0;
    for (j = 0; j < dim; ++j)
        vn += v[j] * vel_rel[j];//accounts for moving interface	
        //vn += v[j] * v_tmp[j]; 	

    for (j = 0; j < dim; ++j)
	    v_tmp[j] -= vn*v[j]; 

    /*
    fprint_general_vector(stdout,"v",v,dim,"\n");
    fprint_general_vector(stdout,"v_tmp",v_tmp,dim,"\n");
    */
}

void KE_CARTESIAN::computeSource()
{
	int i,j,k,l,m,nb,index;
	int index_nb,icrds[MAXD];
	double J,S;
	double **vel = field->vel;
	double vel_nb[2],d_h[2];
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double rho = eqn_params->rho;
	double *Pk = field->Pk;
	double *mu_t = field->mu_t;

	/*find crx*/
	HYPER_SURF *hs;
	STATE *intfc_state;
	INTERFACE* grid_intfc = front->grid_intfc;
	int fr_crx_grid_seg;
	int comp;
	double crx_coords[MAXD],center[MAXD],v_tmp[MAXD];

	switch (dim)
	{
	    case 2:
		for (i = imin; i <= imax; i++)
		for (j = jmin; j <= jmax; j++)
		{
		    icrds[0] = i;
		    icrds[1] = j;
		    index = d_index2d(i,j,top_gmax);
		    getRectangleCenter(index,center);
		    comp = cell_center[index].comp;

		    if (!ifluid_comp(comp)) continue;

            //TODO: Missing a boundary condition for Pk??

		    J = 0.0;
		    /*compute module of the strain rate tensor*/
		    for (l = 0; l < dim; l++)
		    for (m = 0; m < dim; m++)
		    {
                //l derivatives in m direction
                for (nb = 0; nb < 2; nb++)
                {
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                            grid_intfc,icrds,dir[m][nb],comp,
                            (POINTER*)&intfc_state,&hs,crx_coords); 
    
                    d_h[nb] = top_h[m]; 
    
                    if (!fr_crx_grid_seg)
                    {
                        index_nb = next_index_in_dir(icrds,dir[m][nb],top_gmax,dim);
                        vel_nb[nb] = vel[l][index_nb];
                    }
                    else if (fr_crx_grid_seg &&
                        (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                         wave_type(hs) == ELASTIC_BOUNDARY))
                    {
                        setSlipBoundary(icrds,m,nb,comp,hs,intfc_state,field->vel,v_tmp);
                        vel_nb[nb] = v_tmp[l];
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: can't we use vel computed from iF_flowThroughBoundaryState()?
                        if (boundary_state_function_name(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            vel_nb[nb] = vel[l][index];
                        }
                        else
                        {
                            //INLET
                            vel_nb[nb] = intfc_state->vel[l];
                        }
    
                        d_h[nb] = distance_between_positions(center,crx_coords,dim);
                    }
                }
    
                S = (vel_nb[1]- vel_nb[0])/(d_h[1]+d_h[0]);

                //m derivatives in l direction
                for (nb = 0; nb < 2; nb++)
                {
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                            grid_intfc,icrds,dir[l][nb],comp,
                            (POINTER*)&intfc_state,&hs,crx_coords); 
    
                    d_h[nb] = top_h[l]; 
    
                    if (!fr_crx_grid_seg)
                    {
                        index_nb = next_index_in_dir(icrds,dir[l][nb],top_gmax,dim);
                        vel_nb[nb] = vel[m][index_nb];
                    }
                    else if (fr_crx_grid_seg &&
                        (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                         wave_type(hs) == ELASTIC_BOUNDARY))
                    {
                        setSlipBoundary(icrds,l,nb,comp,hs,intfc_state,field->vel,v_tmp);
                        vel_nb[nb] = v_tmp[m];
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: can't we use vel computed from iF_flowThroughBoundaryState()?
                        if (boundary_state_function_name(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            vel_nb[nb] = vel[m][index];
                        }
                        else
                        {
                            //INLET
                            vel_nb[nb] = intfc_state->vel[m];
                        }
     
                        d_h[nb] = distance_between_positions(center,crx_coords,dim);
                    }
 
                }
    
                S += (vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
                J += (S*S);
		    }

            Pk[index] = 0.5*mu_t[index]*J/rho; 
		}	
		break;
	    
        case 3:
		for (i = imin; i <= imax; i++)
		for (j = jmin; j <= jmax; j++)
		for (k = kmin; k <= kmax; k++)
		{
		    icrds[0] = i;
		    icrds[1] = j;
		    icrds[2] = k;
		    index = d_index3d(i,j,k,top_gmax);
		    getRectangleCenter(index,center);
		    comp = cell_center[index].comp;

		    if (!ifluid_comp(comp)) continue;

		    J = 0.0;
		    /*compute module of the strain rate tensor*/
		    for (l = 0; l < dim; l++)
		    for (m = 0; m < dim; m++)
		    {
                //l derivatives in m direction
                for (nb = 0; nb < 2; nb++)
                {
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                            grid_intfc,icrds,dir[m][nb],comp,
                            (POINTER*)&intfc_state,&hs,crx_coords); 
        
                    d_h[nb] = top_h[m]; 
            
                    if (!fr_crx_grid_seg)
                    {
                        index_nb = next_index_in_dir(icrds,dir[m][nb],top_gmax,dim);
                        vel_nb[nb] = vel[l][index_nb];
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                            wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                            wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        setSlipBoundary(icrds,m,nb,comp,hs,intfc_state,field->vel,v_tmp);
                        vel_nb[nb] = v_tmp[l];
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: can't we use vel computed from iF_flowThroughBoundaryState()?
                        if (boundary_state_function_name(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            vel_nb[nb] = vel[l][index];
                        }
                        else
                        {
                            //INLET
                            vel_nb[nb] = intfc_state->vel[l];
                        }

                        d_h[nb] = distance_between_positions(center,crx_coords,dim);
                    }
                }
    
                S = (vel_nb[1]- vel_nb[0])/(d_h[1]+d_h[0]);

                //m derivatives in l direction
                for (nb = 0; nb < 2; nb++)
                {
                    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                            grid_intfc,icrds,dir[l][nb],comp,
                            (POINTER*)&intfc_state,&hs,crx_coords); 
    
                    d_h[nb] = top_h[l]; 
			    
                    if (!fr_crx_grid_seg)
                    {
                        index_nb = next_index_in_dir(icrds,dir[l][nb],top_gmax,dim);
                        vel_nb[nb] = vel[m][index_nb];
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                            wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                            wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        setSlipBoundary(icrds,l,nb,comp,hs,intfc_state,field->vel,v_tmp);
                        vel_nb[nb] = v_tmp[m];
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //TODO: can't we use vel computed from iF_flowThroughBoundaryState()?
                        if (boundary_state_function_name(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            vel_nb[nb] = vel[m][index];

                        }
                        else
                        {
                            //INLET
                            vel_nb[nb] = intfc_state->vel[m];
                        }

                        d_h[nb] = distance_between_positions(center,crx_coords,dim);
                    }
    
                }
    
                S += (vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
                J += (S*S);

            }

            Pk[index] = 0.5*mu_t[index]*J/rho; 
		}	
		
        break;

        default: 
		printf("In computeSource(), Unknown dim = %d\n",dim);
		clean_up(ERROR);
 
    }

    return;
}
                                    
void KE_CARTESIAN::computeTangentialStress(int* icoords, double* tau_wall)
{
	INTERFACE *grid_intfc = front->grid_intfc;
	
	boolean fr_crx_grid_seg;
    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

    STATE* intfc_state;
	HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;

	double coords[MAXD], crx_coords[MAXD];
    double t[MAXD], point[MAXD],nor[MAXD];
	
    double rho = eqn_params->rho;
	double nu = eqn_params->mu/rho;
	double S[MAXD][MAXD];

	double **vel = field->vel;
	double **nor_array = field->nor;

	double d_h[2],vel_nb[2],v_tmp[MAXD];
	int index_nb;

	int index = d_index(icoords,top_gmax,dim);
  	COMPONENT comp = top_comp[index];
	getRectangleCenter(index,coords);

 	if (!ifluid_comp(comp))
    {
        for (int i = 0; i < dim; ++i)
            tau_wall[i] = 0.0;
        return;
    }

	for (int l = 0; l < dim; l++)
	for (int m = 0; m < dim; m++)
	{
        //l derivatives in m direction
	    for (int nb = 0; nb < 2; nb++)
	    {
            fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    grid_intfc,icoords,dir[m][nb],comp,
                    (POINTER*)&intfc_state,&hs,crx_coords);
            
            d_h[nb] = top_h[m]; 
		
            if (!fr_crx_grid_seg)
            {
                index_nb = next_index_in_dir(icoords,dir[m][nb],top_gmax,dim);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                     wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                     wave_type(hs) == ELASTIC_BOUNDARY)
            {
                setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,field->vel,v_tmp);
                vel_nb[nb] = v_tmp[l];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                if (boundary_state_function_name(hs) &&
                        strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                {
                    //OUTLET
                    vel_nb[nb] = vel[l][index];
                }
                else
                {
                    //INLET
                    vel_nb[nb] = intfc_state->vel[l];
                }
            }
        }

        S[l][m] = (vel_nb[1]- vel_nb[0])/(d_h[1]+d_h[0]);
            //S[l][m] = nu*(vel_nb[1]- vel_nb[0])/(d_h[1]+d_h[0]);

        //m derivatives in l direction
	    for (int nb = 0; nb < 2; nb++)
        {
            fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    grid_intfc,icoords,dir[l][nb],comp,
                    (POINTER*)&intfc_state,&hs,crx_coords);
            
            d_h[nb] = top_h[l];
		
            if (!fr_crx_grid_seg)
            {
                index_nb = next_index_in_dir(icoords,dir[l][nb],top_gmax,dim);
                vel_nb[nb] = vel[m][index_nb];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                     wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                     wave_type(hs) == ELASTIC_BOUNDARY)
            {
                setSlipBoundary(icoords,l,nb,comp,hs,intfc_state,field->vel,v_tmp);
                vel_nb[nb] = v_tmp[m];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                if (boundary_state_function_name(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                {
                    //OUTLET
                    vel_nb[nb] = vel[m][index];
                }
                else
                {
                    //INLET
                    vel_nb[nb] = intfc_state->vel[m];
                }
            }
	    }

        S[l][m] += (vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
            //S[l][m] += nu*(vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);

	    if (isnan(S[l][m]))
        {
            printf("wave_type = %d, vel_nb = [%f %f], d_h = [%f %f]\n",
                    wave_type(hs),vel_nb[0],vel_nb[1],d_h[0],d_h[1]);
        }
	}

    
    for (int l = 0; l < dim; ++l)
        nor[l] = nor_array[l][index];

    double ncomp = 0.0;
    double tmp[MAXD] = {0.0};

    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            tmp[i] += nu*S[i][j]*nor[j];
                //tmp[i] += S[i][j]*nor[j];
        }
        
        ncomp += tmp[i]*nor[i];
    }

    for (int i = 0; i < dim; ++i)
    {
        tau_wall[i] = tmp[i] - ncomp*nor[i];
    }
}

/*read k-epsilon parameters*/
void KE_CARTESIAN::read_params(
	char *inname,
	KE_PARAMS *eqn_params)
{
	char string[100];
	FILE* infile;
	infile = fopen(inname,"r");

	/*default parameters*/
	keps_model = RNG;
	eqn_params->delta_k = 0.7194;
	eqn_params->delta_eps = 0.7194;
	eqn_params->Cmu = 0.0845;
	eqn_params->C1 = 1.42;
	eqn_params->C2 = 1.68;
	eqn_params->Cbc = 0.01; /*typicaly 0.003 ~ 0.01*/
	eqn_params->B = 5.2; /*5.2 for smooth wall*/
	eqn_params->y_p = 30.0;
	eqn_params->t0 = 0.0;

	eqn_params->I = 0.01;
	eqn_params->ViscRatioFar = 0.1;
	eqn_params->ViscRatioNear = 100;
    
	eqn_params->U0 = 1.0;
	eqn_params->rho = 1.0;

	/*end default parameter*/
	CursorAfterStringOpt(infile,"Enter type of k-eps model:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'S' || string[0] == 's')
	    keps_model = STANDARD;
	else if (string[0] == 'R' || string[0] == 'r')
	{
	    if (string[1] == 'N' || string[1] == 'n')
	        keps_model = RNG;
	    else
		keps_model = REALIZABLE;
	}

	CursorAfterStringOpt(infile,"Enter turbulent Prandtl number for k:");
	fscanf(infile,"%lf",&eqn_params->delta_k);
	(void) printf("%f\n",eqn_params->delta_k);

	CursorAfterStringOpt(infile,"Enter turbulent Prandtl number for epsilon:");
	fscanf(infile,"%lf",&eqn_params->delta_eps);
	(void) printf("%f\n",eqn_params->delta_eps);

	CursorAfterStringOpt(infile,"Enter C1:");
	fscanf(infile,"%lf",&eqn_params->C1);
	(void) printf("%f\n",eqn_params->C1);

	CursorAfterStringOpt(infile,"Enter C2:");
	fscanf(infile,"%lf",&eqn_params->C2);
	(void) printf("%f\n",eqn_params->C2);

	CursorAfterStringOpt(infile,"Enter Cmu:");
	fscanf(infile,"%lf",&eqn_params->Cmu);
	(void) printf("%f\n",eqn_params->Cmu);

	CursorAfterStringOpt(infile,"Enter Cbc:");
	fscanf(infile,"%lf",&eqn_params->Cbc);
	(void) printf("%f\n",eqn_params->Cbc);
	
    CursorAfterString(infile,"Enter I:");
	fscanf(infile,"%lf",&eqn_params->I);
	(void) printf("%f\n",eqn_params->I);

    CursorAfterString(infile,"Enter ViscRatioFar:");
	fscanf(infile,"%lf",&eqn_params->ViscRatioFar);
	(void) printf("%f\n",eqn_params->ViscRatioFar);

    CursorAfterString(infile,"Enter ViscRatioNear:");
	fscanf(infile,"%lf",&eqn_params->ViscRatioNear);
	(void) printf("%f\n",eqn_params->ViscRatioNear);

	CursorAfterString(infile,"Enter l0:");
	fscanf(infile,"%lf",&eqn_params->l0);
	(void) printf("%f\n",eqn_params->l0);

	CursorAfterString(infile,"Enter U0:");
	fscanf(infile,"%lf",&eqn_params->U0);
	(void) printf("%f\n",eqn_params->U0);

	CursorAfterStringOpt(infile,"Enter mu0:");
	fscanf(infile,"%lf",&eqn_params->mu0);
	(void) printf("%f\n",eqn_params->mu0);

	CursorAfterStringOpt(infile,"Enter y+:");
	fscanf(infile,"%lf",&eqn_params->y_p);
	(void) printf("%f\n",eqn_params->y_p);

    CursorAfterStringOpt(infile,"Enter B:");
    fscanf(infile,"%lf",&eqn_params->B);
    (void) printf("%f\n",eqn_params->B);

    CursorAfterStringOpt(infile,"Enter time to active turbulence model:");
    fscanf(infile,"%lf",&eqn_params->t0);
    (void) printf("%f\n",eqn_params->t0);

    fclose(infile);
	return;
}

//////////////////////////////////////////////////////
/////           DOMAIN STATISTICS                ////
////////////////////////////////////////////////////

double KE_CARTESIAN::computeEddyViscosityMean()
{
    switch (dim)
    {
        case 2:
            return computeEddyViscosityMean2d();
            break;
        case 3:
            return computeEddyViscosityMean3d();
            break;
        default:
            printf("ERROR: Invalid Dimension!\n");
            LOC(); clean_up(EXIT_FAILURE);
    }
}

double KE_CARTESIAN::computeEddyViscosityMean2d()
{
    int index;
    COMPONENT comp;
    double* mu_t = field->mu_t;

    int count = 0;
    double mut_avg = 0.0;

    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    {
        index = d_index2d(i,j,top_gmax);
        comp = cell_center[index].comp;
        if (!ifluid_comp(comp)) continue;

        mut_avg += mu_t[index];
        count++;
    }

    mut_avg /= (double)count;
    return mut_avg;
}

double KE_CARTESIAN::computeEddyViscosityMean3d()
{
    int index;
    COMPONENT comp;
    double* mu_t = field->mu_t;

    int count = 0;
    double mut_avg = 0.0;

    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    for (int k = kmin; k <= kmax; ++k)
    {
        index = d_index3d(i,j,k,top_gmax);
        comp = cell_center[index].comp;
        if (!ifluid_comp(comp)) continue;

        mut_avg += mu_t[index];
        count++;
    }
    
    mut_avg /= (double)count;
    return mut_avg;
}

double KE_CARTESIAN::computeEddyViscosityVariance()
{
    switch (dim)
    {
        case 2:
            return computeEddyViscosityVariance2d();
            break;
        case 3:
            return computeEddyViscosityVariance3d();
            break;
        default:
            printf("ERROR: Invalid Dimension!\n");
            LOC(); clean_up(EXIT_FAILURE);
    }
}

double KE_CARTESIAN::computeEddyViscosityVariance2d()
{
    int index;
    COMPONENT comp;
    double* mu_t = field->mu_t;

    int count = 0;
    double mut_var = 0.0;

    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    {
        index = d_index2d(i,j,top_gmax);
        comp = cell_center[index].comp;
        if (!ifluid_comp(comp)) continue;

        mut_var += sqr(mu_t[index]);
        count++;
    }

    mut_var /= (double)count;

    double mut_avg = computeEddyViscosityMean2d();
    mut_var -= sqr(mut_avg);
    
    return mut_var;
}

double KE_CARTESIAN::computeEddyViscosityVariance3d()
{
    int index;
    COMPONENT comp;
    double* mu_t = field->mu_t;

    int count = 0;
    double mut_var = 0.0;

    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    for (int k = kmin; k <= kmax; ++k)
    {
        index = d_index3d(i,j,k,top_gmax);
        comp = cell_center[index].comp;
        if (!ifluid_comp(comp)) continue;

        mut_var += sqr(mu_t[index]);
        count++;
    }
    
    mut_var /= (double)count;
    
    double mut_avg = computeEddyViscosityMean3d();
    mut_var -= sqr(mut_avg);
    
    return mut_var;
}

double KE_CARTESIAN::computeEddyViscosityStdDeviation()
{
    double mut_var = computeEddyViscosityVariance();
    return sqrt(mut_var);
}





