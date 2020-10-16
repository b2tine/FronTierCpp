/*******************************************************************
 * 	keps.cpp	
 *******************************************************************/
/*
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include "keps.h"

static void printField2d(double*,const char*,int*,int*,int*);
static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};
static int find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
static int next_index_in_dir(int*,GRID_DIRECTION,int,int*);
static std::string dir2String(GRID_DIRECTION dir);

//----------------------------------------------------------------
//		KE_RECTANGLE
//----------------------------------------------------------------

KE_RECTANGLE::KE_RECTANGLE()
    : index(-1), comp(-1)
{
}

void KE_RECTANGLE::setCoords(
	double *crds,
	int dim)
{
	for (int i = 0; i < dim; ++i)
	    coords[i] = crds[i];
}
//--------------------------------------------------------------------------
// 		KE_CARTESIAN
//--------------------------------------------------------------------------

KE_CARTESIAN::~KE_CARTESIAN()
{
}

KE_CARTESIAN::KE_CARTESIAN(Front &front)
    : front(&front)
{
}

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
	KE_RECTANGLE       rectangle;

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
	int i;

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

        for (i = 0; i < size; i++)
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
                    int ii,jj,ic[MAXD],index;
                    icoords = cell_center[i].icoords;
                    n = 0;
                    for (ii = 0; ii < dim; ++ii)
                    {
                        for (jj = 0; jj < dim; ++jj) ic[jj] = icoords[jj];
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
	double coords[MAXD];
	INTERFACE *intfc = front->interf;
	POINT *p;
    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	short unsigned int seed[3] = {2,72,7172};

	FT_MakeGridIntfc(front);
    setDomain();

    //compute mixing length limits
    lmax = HUGE;
    lmin = -HUGE;
    for (int i = 0; i < dim; ++i)
    {
        double domain_size = top_U[i] - top_L[i];
        if (domain_size < lmax)
            lmax = domain_size;
        if (top_h[i] > lmin)
            lmin = top_h[i];
    }
    
    //lmax /= 3.0;
    //lmin *= 2.0;

    if (eqn_params->l0 < lmin)
        eqn_params->l0 = lmin;


	double k0 = sqr(eqn_params->mu0/eqn_params->rho/eqn_params->l0);
	double eps0 = eqn_params->Cmu*pow(k0,1.5)/eqn_params->l0;
    
    //save the initial conditions for activation time
    eqn_params->k0 = k0;
    eqn_params->eps0 = eps0;

    //TODO: These should already be zeroed -- doing for sanity check
	for (int i = 0; i < cell_center.size(); ++i)
	{
        field->k[i] = 0.0;
	    field->eps[i] = 0.0;
	    field->mu_t[i] = 0.0;
    }

	printf("mu0 = %e\n",eqn_params->mu0);
    printf("lmin = %e, l0 = %e, lmax = %e\n",lmin,eqn_params->l0,lmax);
    printf("k0 = %e, eps0 = %e\n",k0,eps0);
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
        field->k[i] = k0;
	    field->eps[i] = eps0;
	    
        //double nu_t = std::max(Cmu*sqr(k0)/eps0,0.0001*nu);
        //field->mu_t[i] = nu_t*rho;
	        //field->mu_t[i] = mu0;
        field->mu_t[i] = 0.0;
        
        if (keps_model == REALIZABLE)
            field->Cmu[i] = Cmu;
	}
}


//TODO: update this method
//
//NOTE: not called anywhere
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
	        if (ref_p[i] < front->rect_grid->L[i] + 2*tol[i] ||
                ref_p[i] > front->rect_grid->U[i] - 2*tol[i])
            {
                skip = YES;
            }
	    }
	    
        if (skip) continue;
	    
        //TODO: force computation should include effects of shear stress from
        //      turbulence model + wall functions (see to computeDiffusionCN() todos).
        
        //if (wave_type(*c) == NEUMANN_BOUNDARY || wave_type(*c) == ELASTIC_BOUNDARY)
        if (wave_type(*c) == NEUMANN_BOUNDARY)
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
	fprintf(vfile,"%f %f %f %f\n",front->time,force[0],force[1],force[2]);
	fclose(vfile);
}

//TODO: Note the difference between Incompress_Solver_Smooth_Basis::setIndexMap()
//      Can probably get ride of I_to_ijk etc.
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

void KE_CARTESIAN::computeKE()
{
	COMPONENT sub_comp[2];
	sub_comp[0] = SOLID_COMP;
	sub_comp[1] = LIQUID_COMP2;
	
	for (int i = 0; i < 2; ++i)
	{
        if(sub_comp[i] == SOLID_COMP) continue;
        
        //TODO: Check iupper and ilower vals per the comment at top of file
        setGlobalIndex(sub_comp[i]);

        implicitComputeK(sub_comp[i]);
        implicitComputeE(sub_comp[i]);
    }
}

void KE_CARTESIAN::implicitComputeK(COMPONENT sub_comp)
{
    switch (dim)
    {
        case 2:
            implicitComputeK2d(sub_comp);
            break;
        case 3:
            printf("implicitComputeK3d() not implemented!\n"); LOC();
            clean_up(EXIT_FAILURE);
            //implicitComputeK3d(sub_comp);
            break;
        default:
            printf("Invalid Dimension! dim = %d\n",dim); LOC();
            clean_up(EXIT_FAILURE);
    }
}

void KE_CARTESIAN::implicitComputeE(COMPONENT sub_comp)
{
    switch (dim)
    {
        case 2:
            implicitComputeE2d(sub_comp);
            break;
        case 3:
            printf("implicitComputeE3d() not implemented!\n"); LOC();
            clean_up(EXIT_FAILURE);
            //implicitComputeE3d(sub_comp);
            break;
        default:
            printf("Invalid Dimension! dim = %d\n",dim); LOC();
            clean_up(EXIT_FAILURE);
    }
}

void KE_CARTESIAN::computeAdvection()
{
	COMPONENT sub_comp[2];
	sub_comp[0] = SOLID_COMP;
	sub_comp[1] = LIQUID_COMP2;
	
	for (int i = 0; i < 2; ++i)
	{
        if(sub_comp[i] == SOLID_COMP) continue;
        
        setGlobalIndex(sub_comp[i]);
        computeAdvectionK(sub_comp[i]);

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

//TODO: Remove
void KE_CARTESIAN::updateGamma()
{
    switch (dim)
    {
        case 2:
            updateGamma2d();
            break;
        case 3:
            printf("3d code not implemented!\n"); LOC();
            clean_up(EXIT_FAILURE);
            //updateGamma3d();
            break;
        default:
            printf("Invalid Dimension! dim = %d\n",dim); LOC();
            clean_up(EXIT_FAILURE);
    }
}

//TODO: Remove
void KE_CARTESIAN::updateGamma2d()
{
    double *gamma = field->gamma;
    double *K = field->k;
	double *mu_t = field->mu_t;
	double rho = eqn_params->rho;
    double Cmu = eqn_params->Cmu;

    //double icoords[MAXD] = {0};

    for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax; ++i)
    {
        //icoords[0] = i;
        //icoords[1] = j;
        int index = d_index2d(i,j,top_gmax);
        
        COMPONENT comp = top_comp[index];
        if (!ifluid_comp(comp))
        {
            gamma[index] = 0.0;
            continue;
        }

        //TODO: may need a gamma_k and gamma_eps that are updated
        //      slightly differently with values from different time levels
        double nu_t = mu_t[index]/rho;
        if (fabs(nu_t) < MACH_EPS)
            gamma[index] = 0.0;
        else
            gamma[index] = Cmu*K[index]/nu_t;

        gamma[index] = std::max(gamma[index],0.0);
    }

    //TODO: need to scatter?
    /*
    scatMeshArray();

    for (int j = 0; j <= top_gmax[1]; ++j)
    for (int i = 0; i <= top_gmax[0]; ++i)
    {
        ic = d_index2d(i,j,top_gmax);
        comp = cell_center[ic].comp;
        if (comp == sub_comp)
            gamma[ic] = array[ic];
    }
    */
}

//TODO: Remove
void KE_CARTESIAN::updateGamma3d()
{
    double *gamma = field->gamma;
    double *K = field->k;
	double *mu_t = field->mu_t;
	double rho = eqn_params->rho;
    double Cmu = eqn_params->Cmu;

    //double icoords[MAXD] = {0};

    for (int k = kmin; k <= kmax; ++k)
    for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax; ++i)
    {
        //icoords[0] = i;
        //icoords[1] = j;
        //icoords[2] = k;
        int index = d_index3d(i,j,k,top_gmax);
        
        COMPONENT comp = top_comp[index];
        if (!ifluid_comp(comp))
        {
            gamma[index] = 0.0;
            continue;
        }

        //TODO: may need a gamma_k and gamma_eps that are updated
        //      slightly differently with values from different time levels
        double nu_t = mu_t[index]/rho;
        if (fabs(nu_t) < MACH_EPS)
            gamma[index] = 0.0;
        else
            gamma[index] = Cmu*K[index]/nu_t;

        gamma[index] = std::max(gamma[index],0.0);
    }

    //TODO: need to scatter?
    //
    /*
    scatMeshArray();

    for (int k = 0; k <= top_gmax[2]; ++k)
    for (int j = 0; j <= top_gmax[1]; ++j)
    for (int i = 0; i <= top_gmax[0]; ++i)
    {
        ic = d_index3d(i,j,k,top_gmax);
        comp = cell_center[ic].comp;
        if (comp == sub_comp)
            gamma[ic] = array[ic];
    }
    */
}

/*
void KE_CARTESIAN::explicitComputeKE(COMPONENT sub_comp)
{
    switch (dim)
    {
        case 2:
            explicitComputeKE2d(sub_comp);
            break;
        case 3:
            printf("explicitComputeKE3d() not implemented!\n"); LOC();
            clean_up(EXIT_FAILURE);
            //explicitComputeKE3d(sub_comp);
            break;
        default:
            printf("Invalid Dimension! dim = %d\n",dim); LOC();
            clean_up(EXIT_FAILURE);
    }
}
*/

/*
//TODO: Debug/Implement and extend to 3d if successful
void KE_CARTESIAN::explicitComputeKE2d(COMPONENT sub_comp)
{
    explicitComputeK2d(sub_comp);
    explicitComputeE2d(sub_comp);
}
*/

/*
//DOOMED
void KE_CARTESIAN::explicitComputeK2d(COMPONENT sub_comp)
{
	INTERFACE* grid_intfc = front->grid_intfc;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	STATE *intfc_state;

    int icoords[MAXD], icoords_nb[MAXD];
    double crx_coords[MAXD];
	double nor[MAXD];
    COMPONENT comp;
    
    boolean fr_crx_grid_seg;
    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    //double *gamma = field->gamma;
	double *Pk = field->Pk;
    double *K = field->k;
	double *mu_t = field->mu_t;

	double rho = eqn_params->rho;
	double nu = eqn_params->mu/eqn_params->rho;
    double Cmu = eqn_params->Cmu;
    double C1 = eqn_params->C1;
    double C2 = eqn_params->C2;
    double Cbc = eqn_params->Cbc;
	double sigma_k = eqn_params->sigma_k;
	double sigma_eps = eqn_params->sigma_eps;
    double l0 = eqn_params->l0;
	

    start_clock("computeK");
    if (debugging("trace")) printf("Entering explicitComputeK2d()\n");

        //setIndexMap(sub_comp);

    if (debugging("trace"))
    {
        int domain_size = 1;
        printf("ilower = %d  iupper = %d\n",ilower,iupper);
        for (int i = 0; i < dim; ++i)
            printf("domain_size = %d\n",domain_size);
    }

    for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax; ++i)
    {
        icoords[0] = i;
        icoords[1] = j;
        int ic = d_index2d(i,j,top_gmax);
        
            //I = ij_to_I[i][j];
        
        comp = top_comp[ic];
        if (comp != sub_comp) 
        {
            Karray[ic] = 0.0;
            continue;
        }
        

        //gammaK updated using previous values (time n)
        double nu_t = mu_t[ic]/rho;
        double gammaK = std::max(Cmu*K[ic]/nu_t,0.0);

        //K[ic] at time n+1 written into Karray[ic] to prevent data race
        Karray[ic] = (1.0 - gammaK*m_dt)*K[ic] + Pk[ic]*m_dt;

        for (int l = 0; l < dim; ++l)
        {
            double alpha = 0.5*m_dt*field->vel[l][ic]/top_h[l];
            double beta = m_dt*(nu + nu_t/sigma_k)/sqr(top_h[l]);
                //double beta = m_dt*nu_t/sigma_k/sqr(top_h[l]);
        
            Karray[ic] -= 2.0*beta*K[ic];
        
            for (int nb = 0; nb < 2; ++nb)
            {
                for (int n = 0; n < dim; ++n)
                    icoords_nb[n] = icoords[n];

                int sign = pow(-1,nb);

                icoords_nb[l] =
                    (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;

                int icnb = d_index(icoords_nb,top_gmax,dim);

                fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                        grid_intfc,icoords,dir[l][nb],comp,
                        (POINTER*)&intfc_state,&hs,crx_coords);
                
                if (!fr_crx_grid_seg) 
                {
                    //no cross
                    Karray[ic] += sign*alpha*K[icnb] + beta*K[icnb];
                }
                else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    // Inlet/Outlet on rectangular boundaries only
                    if (boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                    {
                        //OUTLET: nor dot grad(k) = 0
                        double K_nb = K[ic];
                        Karray[ic] += sign*alpha*K_nb + beta*K_nb;
                    }
                    else
                    {
                        //INLET: k = Cbc*|u|^2 where Cbc \in [0.003,0.01]
                        double sqrmag_vel = Dot2d(intfc_state->vel,intfc_state->vel);
                        double K_nb = eqn_params->Cbc*sqrmag_vel;
                        Karray[ic] += sign*alpha*K_nb + beta*K_nb;
                    }
                }
                else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                {
                    double v_tan[MAXD];
                    setSlipBoundary(icoords,l,nb,comp,hs,intfc_state,field->vel,v_tan);

                    double mag_vtan = Magd(v_tan,dim);

                    //friction velocity
                    double u_t = std::max(
                            pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                            mag_vtan/eqn_params->y_p);
                    
                    //double u_t = 0.41*mag_vtan/log(eqn_params->E*eqn_params->y_p);
                    double K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                    Karray[ic] += sign*alpha*K_nb + beta*K_nb;
                }
                else
                {
                    printf("UNKNOWN BOUNDARY TYPE!\n"); LOC();
                    clean_up(EXIT_FAILURE);
                }
            
            }//nb

        }//l
    
    }//i,j
            
        
    
	FT_ParallelExchGridArrayBuffer(Karray,front,NULL);

    //Write Karray back into K
    for (int j = 0; j <= top_gmax[1]; ++j)
    for (int i = 0; i <= top_gmax[0]; ++i)
    {
        int ic = d_index2d(i,j,top_gmax);
        comp = cell_center[ic].comp;
        if (comp == sub_comp)
        {
            K[ic] = Karray[ic];
        }
    }
    
    
    if (debugging("trace")) printf("Leaving explicitComputeK2d()\n");
    stop_clock("computeK");
}*/       /* end explicitComputeK2d */
   
/*
//DOOMED
void KE_CARTESIAN::explicitComputeE2d(COMPONENT sub_comp)
{
	INTERFACE* grid_intfc = front->grid_intfc;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	STATE *intfc_state;

    int icoords[MAXD], icoords_nb[MAXD];
    double crx_coords[MAXD];
	double nor[MAXD];
    COMPONENT comp;
    
    boolean fr_crx_grid_seg;
    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    //double *gamma = field->gamma;
	double *Pk = field->Pk;
    double *K = field->k;
    double *E = field->eps;
	double *mu_t = field->mu_t;

	double rho = eqn_params->rho;
	double nu = eqn_params->mu/rho;
    double Cmu = eqn_params->Cmu;
    double C1 = eqn_params->C1;
    double C2 = eqn_params->C2;
    double Cbc = eqn_params->Cbc;
	double sigma_k = eqn_params->sigma_k;
	double sigma_eps = eqn_params->sigma_eps;
    double l0 = eqn_params->l0;
	

    start_clock("computeK");
    if (debugging("trace")) printf("Entering explicitComputeK2d()\n");

        //setIndexMap(sub_comp);

    if (debugging("trace"))
    {
        int domain_size = 1;
        printf("ilower = %d  iupper = %d\n",ilower,iupper);
        for (int i = 0; i < dim; ++i)
            printf("domain_size = %d\n",domain_size);
    }

    for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax; ++i)
    {
        icoords[0] = i;
        icoords[1] = j;
        int ic = d_index2d(i,j,top_gmax);
        
            //I = ij_to_I[i][j];
        
        comp = top_comp[ic];
        if (comp != sub_comp) 
        {
            Karray[ic] = 0.0;
            Earray[ic] = 0.0;
            continue;
        }
        
        //gammaE and Peps updated using updated values (time n+1),
        //except for E which has not been updated yet and nu_t which
        //is only needed to recover the common term of Pk and Peps.
        double nu_t = mu_t[ic]/rho;
        double gammaE = std::max(C2*E[ic]/K[ic],0.0);
        double Peps = std::max(C1*K[ic]*Pk[ic]/nu_t,0.0);

        //E[ic] at time n+1 written into Earray[ic] to prevent data race
        Earray[ic] = (1.0 - gammaE*m_dt)*E[ic] + Peps*m_dt;

        for (int l = 0; l < dim; ++l)
        {
            double alpha = 0.5*m_dt*field->vel[l][ic]/top_h[l];
            double beta = m_dt*(nu + nu_t/sigma_eps)/sqr(top_h[l]);
                //double beta = m_dt*nu_t/sigma_eps/sqr(top_h[l]);
        
            Earray[ic] -= 2.0*beta*E[ic];
        
            for (int nb = 0; nb < 2; ++nb)
            {
                for (int n = 0; n < dim; ++n)
                    icoords_nb[n] = icoords[n];

                int sign = pow(-1,nb);

                icoords_nb[l] =
                    (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;

                int icnb = d_index(icoords_nb,top_gmax,dim);

                fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                        grid_intfc,icoords,dir[l][nb],comp,
                        (POINTER*)&intfc_state,&hs,crx_coords);
                
                if (!fr_crx_grid_seg) 
                {
                    //no cross
                    Earray[ic] += sign*alpha*E[icnb] + beta*E[icnb];
                }
                else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    // Inlet/Outlet on rectangular boundaries only
                    if (boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                    {
                        //OUTLET: nor dot grad(eps) = 0
                        double E_nb = E[ic];
                        Earray[ic] += sign*alpha*E_nb + beta*E_nb;
                    }
                    else
                    {
                        //INLET: k = Cbc*|u|^2 where Cbc \in [0.003,0.01]
                        double sqrmag_vel = Dot2d(intfc_state->vel,intfc_state->vel);
                        double K_nb = eqn_params->Cbc*sqrmag_vel;

                        // eps = Cmu*k^(3/2)/l0
                        double E_nb = Cmu*pow(K_nb,1.5)/l0;
                        Earray[ic] += sign*alpha*E_nb + beta*E_nb;
                    }
                }
                else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                {
                    double v_tan[MAXD];
                    setSlipBoundary(icoords,l,nb,comp,hs,intfc_state,field->vel,v_tan);

                    double mag_vtan = Magd(v_tan,dim);

                    //friction velocity
                    double u_t = std::max(
                            pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                            mag_vtan/eqn_params->y_p);

                    //double u_t = 0.41*mag_vtan/log(eqn_params->E*eqn_params->y_p);
                    double E_nb = pow(u_t,4.0)/(0.41*eqn_params->y_p*nu);
                    Earray[ic] += sign*alpha*E_nb + beta*E_nb;
                }
                else
                {
                    printf("UNKNOWN BOUNDARY TYPE!\n"); LOC();
                    clean_up(EXIT_FAILURE);
                }
            
            }//nb

        }//l
    
    }//i,j
            
        
    
	FT_ParallelExchGridArrayBuffer(Earray,front,NULL);

    //Write Karray and Earray back into K and E
    for (int j = 0; j <= top_gmax[1]; ++j)
    for (int i = 0; i <= top_gmax[0]; ++i)
    {
        int ic = d_index2d(i,j,top_gmax);
        comp = cell_center[ic].comp;
        if (comp == sub_comp)
        {
            E[ic] = Earray[ic];
        }
    }
    
    
    if (debugging("trace")) printf("Leaving explicitComputeE2d()\n");
    stop_clock("computeK");
}*/       /* end explicitComputeE2d */
           

//TODO: Rename and implement Crank Nicholson scheme (2nd order time)
//implicit time explicit space
void KE_CARTESIAN::computeAdvectionK(COMPONENT sub_comp)
{
    int i,j,k,ll,ic,icn,I,I_nb,icoords[MAXD];
    int gmin[MAXD],ipn[MAXD];
    double crx_coords[MAXD];
	double nor[MAXD];
    double K0,K_nb,D,lambda,coeff,coeff_nb,rhs;
    COMPONENT comp;
    PETSc solver;
    double *x;
    int num_iter = 0;
    double rel_residual = 0;
    boolean fr_crx_grid_seg;
    const GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

    double *K = field->k;
	double *Pk = field->Pk;
	double *mu_t = field->mu_t;
	double sigma_k = eqn_params->sigma_k;
	double rho = eqn_params->rho;
	double nu = eqn_params->mu/eqn_params->rho;
	double y_pp,dist,center[MAXD];
    double v[MAXD],v_wall[MAXD],crds_wall[MAXD],k_wall;
    double eta;
	double Ut,Ut_old,vn;
	int    niter;
	
    /*For boundary state*/
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *intfc_state;
	INTERFACE* grid_intfc = front->grid_intfc;
	double coords[MAXD];
	boolean if_adj_pt;
	double point[MAXD],t[MAXD];
	double Cmu;

        start_clock("computeAdvectionK");
        if (debugging("trace")) printf("Entering computeAdvectionK()\n");

        for (i = 0; i < dim; ++i) gmin[i] = 0;

        setIndexMap(sub_comp);
        if (debugging("trace"))
        {
            //TODO: Unfinished feature that comment at top of file refers to?
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
                comp = top_comp[ic];
                I = ij_to_I[i][j];
                if (comp != sub_comp)
                    continue;
                
		/*fully implicit to preserve positivity*/
        K0 = K[ic];
		Cmu = eqn_params->Cmu;
		rhs = K0 + m_dt*Pk[ic];

		if (keps_model == REALIZABLE)
		  coeff = 1.0 + m_dt*std::max(field->eps[ic],0.0);
		else
		  coeff = 1.0 + m_dt*std::max(Cmu*K0*rho/mu_t[ic],0.0);
		
        if (isinf(coeff) || isnan(coeff) || isinf(rhs) || isnan(rhs))
		{
		    printf("In computeAdvectionK(): ");
		    printf("icoords[%d %d], index = %d\n",i,j,ic);
		    printf("coeff=%f, K=%e, E=%e, mu_t=%e, Pk=%e\n",
			coeff,K0,field->eps[ic],mu_t[ic],Pk[ic]);
		    clean_up(ERROR);
		}
                for (int l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (int l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
                
                //TODO: this looks like a bug -- adding nu!
               	D = nu + mu_t[ic]/eqn_params->sigma_k/rho;

                for (int l = 0; l < dim; ++l)
                {
                    lambda = D*m_dt/sqr(top_h[l]);
                    eta = v[l]*m_dt/(top_h[l]); //upwind difference
                    double eta_p = std::max(eta, 0.0);
                    double eta_m = std::min(eta, 0.0);
                        //eta = v[l]*m_dt/(2.0*top_h[l]);
                    coeff += 2*lambda;

                    for (int m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
                        I_nb = ij_to_I[ipn[0]][ipn[1]];
            
                        fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icoords,dir[l][m],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords);
                        
                        coeff += ((m == 0) ? eta_p : -eta_m); //upwind
                        if (!fr_crx_grid_seg) 
                        {
                            coeff_nb = -lambda;
                            coeff_nb += (m == 0) ? -eta_p : eta_m;
                                //coeff_nb = -lambda + (pow(-1,m+1)*eta);
                            solver.Add_A(I,I_nb,coeff_nb);
                        }
			/*else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                    wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                    wave_type(hs) == ELASTIC_BOUNDARY)*/
			else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                    wave_type(hs) == MOVABLE_BODY_BOUNDARY)
			{
                double v_tan[MAXD];
                setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_tan);

                double mag_vtan = Magd(v_tan,dim);

                //friction velocity
                double u_t = std::max(
                        pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                        mag_vtan/eqn_params->y_p);
                
                K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 

                //use wall function for friction velocity u_t

                /*
				setTKEatWall(icoords,l,m,comp,
						hs,intfc_state,K,&K_nb);
				rhs += lambda*K_nb - (pow(-1,m+1)*eta)*K_nb;
                */

                //TODO: Use friction velocity to compute the wall shear
                //      stress acting in opposition to the the local
                //      velocity
                
                /*
                boolean status;
                status = FT_NormalAtGridCrossing(front,icoords,dir[l][m],
                        comp,nor,&hs,crx_coords);
                
                
                vn = 0.0;
                double* vel = intfc_state->vel;
                for (int kk = 0; kk < dim; ++kk)
                    vn += (v[kk] - vel[kk])*nor[kk];
                
                double v_tan[MAXD] = {0.0};
                for (int kk = 0; kk < dim; ++kk)
                    v_tan[kk] = v[kk] - vn*nor[kk];
                double mag_vtan = Magd(v_tan,dim);

                //friction velocity
                double u_t = std::max(
                        pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                        mag_vtan/eqn_params->y_p);
                
                K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 

                //tangential stress acting opposite of flow direction
                double unit_tan[MAXD] = {0.0};
                if (mag_vtan > 0.0)
                {
                    for (int j = 0; j < dim; ++j)
                        unit_tan[j] = v_tan[j]/mag_vtan;
                }

                double mag_tanstress = rho*u_t*u_t;
                for (int j = 0; j < dim; ++j)
                    intfc_state->tan_stress[j] = -mag_tanstress*unit_tan[j];
                */
			}
			else if (wave_type(hs) == DIRICHLET_BOUNDARY)
			{
			    if (boundary_state_function(hs) &&
				strcmp(boundary_state_function_name(hs),
                            	"flowThroughBoundaryState") == 0)
			    {
                    //Outlet
				    K_nb = K0;
                    rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
				        //rhs += lambda*K_nb - (pow(-1,m+1)*eta)*K_nb;
			    }
			    else
			    {
                    //Inlet
			        K_nb = eqn_params->Cbc
				     * (sqr(intfc_state->vel[0])
				     +  sqr(intfc_state->vel[1]));
                    rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
				        //rhs += lambda*K_nb - (pow(-1,m+1)*eta)*K_nb;
			    }
			}
			else
            {
                //printf("Unknows boundary condition! \n");
                //clean_up(ERROR);
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
                ic = d_index3d(i,j,k,top_gmax);
                comp = top_comp[ic];
                I = ijk_to_I[i][j][k];
                if (comp != sub_comp) continue;

                K0 = K[ic];
		        Cmu = eqn_params->Cmu;
                rhs = K0 + m_dt*Pk[ic];

            /*fully implicit to preserve positivity*/
            if (keps_model == REALIZABLE)
                coeff = 1.0 + m_dt*std::max(field->eps[ic],0.0);
            else
                coeff = 1.0 + m_dt*std::max(Cmu*K0*rho/mu_t[ic],0.0);

            if (isinf(coeff) || isnan(coeff))
            {
                printf("In computeAdvectionK(): ");
		        printf("icoords[%d %d %d], index = %d\n",i,j,k,ic);
                printf("coeff=%f, K=%e, E=%e, mu_t=%e, Pk=%e\n",
                coeff,K0,field->eps[ic],mu_t[ic],Pk[ic]);
                clean_up(ERROR);
            }
            
            for (int l = 0; l < dim; ++l)
                v[l] = 0.0;

            if (field->vel != NULL)
            {
                for (int l = 0; l < dim; ++l)
                    v[l] = field->vel[l][ic];
            }

            D = nu + mu_t[ic]/eqn_params->sigma_k/rho;
            
            for (int l = 0; l < dim; ++l)
            {
                lambda = D*m_dt/sqr(top_h[l]);
                //eta = 0.5*v[l]*m_dt/(top_h[l]);
                eta = v[l]*m_dt/(top_h[l]); //upwind difference
                double eta_p = std::max(eta, 0.0);
                double eta_m = std::min(eta, 0.0);
                coeff += 2*lambda;

                for (int m = 0; m < 2; ++m)
                {
                    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                    icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                    I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
        

                    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                        grid_intfc,icoords,dir[l][m],comp,
                                        (POINTER*)&intfc_state,&hs,crx_coords);
        
                    coeff += eta;
                    //coeff += ((m == 0) ? eta_p : -eta_m); //upwind
                    
                    if (!fr_crx_grid_seg) 
                    {
                        coeff_nb = -lambda;
                        coeff_nb += (m == 0) ? -eta_p : eta_m;
                        solver.Add_A(I,I_nb,coeff_nb);
                    }
                    /*else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                             wave_type(hs) == ELASTIC_BOUNDARY)*/
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                    {
                        double v_tan[MAXD];
                        setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_tan);

                        double mag_vtan = Magd(v_tan,dim);

                        //friction velocity
                        double u_t = std::max(
                                pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                                mag_vtan/eqn_params->y_p);
                        
                        K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                        rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 

                        /*
                        //TODO: Use friction velocity to compute the wall shear
                        //      stress acting in opposition to the the local
                        //      velocity
                        
                        boolean status;
                        status = FT_NormalAtGridCrossing(front,icoords,dir[l][m],
                                comp,nor,&hs,crx_coords);

                        //use wall function for friction velocity u_t
                        
                        vn = 0.0;
                        double* vel = intfc_state->vel;
                        for (int kk = 0; kk < 3; ++kk)
                            vn += (v[kk] - vel[kk])*nor[kk];
                                //vn += vel[kk]*nor[kk];
                        
                        double v_tan[3];
                        for (int kk = 0; kk < 3; ++kk)
                            v_tan[kk] = v[kk] - vn*nor[kk];
                        double mag_vtan = Magd(v_tan,dim);

                        //friction velocity
                        double u_t = std::max(
                                pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                                mag_vtan/eqn_params->y_p);
                        
                        K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                        rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 

                        //tangential stress acting opposite of flow direction
                        double unit_tan[MAXD] = {0.0};
                        if (mag_vtan > 0.0)
                        {
                            for (int j = 0; j < dim; ++j)
                                unit_tan[j] = v_tan[j]/mag_vtan;
                        }

                        double mag_tanstress = rho*u_t*u_t;
                        for (int j = 0; j < dim; ++j)
                            intfc_state->tan_stress[j] = -mag_tanstress*unit_tan[j];
                        */
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (boundary_state_function(hs) &&
                        strcmp(boundary_state_function_name(hs),
                                        "flowThroughBoundaryState") == 0)
                        {
                            //Outlet
                            K_nb = K0;
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                        }
                        else
                        {
                            //Inlet
                            K_nb = eqn_params->Cbc
                             * (sqr(intfc_state->vel[0])
                             +  sqr(intfc_state->vel[1])
                             +  sqr(intfc_state->vel[2]));
                            rhs += lambda * K_nb + ((m == 0) ? eta_p*K_nb : -eta_m*K_nb); 
                        }
                    }
                    else
                    {
                        //printf("Unknown boundary condition %d! \n",wave_type(hs));
                        //clean_up(ERROR);
                    }
                }

            }

            solver.Add_A(I,I,coeff);
            solver.Add_b(I,rhs);

            }

            break;
        }

        stop_clock("set_coefficients");
        start_clock("petsc_solve");
        solver.SetMaxIter(10000);
        //solver.SetMaxIter(500);
        solver.SetTol(1e-8);
        solver.Solve();
	    solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

        if (debugging("PETSc"))
        {
            (void) printf("KE_CARTESIAN::computeAdvectionK: "
                        "num_iter = %d, rel_residual = %g \n",
                        num_iter, rel_residual);
        }

        FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
        solver.Get_x(x);
        stop_clock("petsc_solve");

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

        if (debugging("trace")) printf("Leaving computeAdvectionK()\n");
        stop_clock("computeAdvectionK");
}       /* end computeAdvectionK */

double KE_CARTESIAN::computePointFieldC2_RNG(int* icoords)
{
	int index;
	double mu_t,E0,K0,r,S;
	index = d_index(icoords,top_gmax,dim);
	mu_t = field->mu_t[index];
	K0 = field->k[index];
	E0 = field->eps[index];
	S = computePointFieldStrain(icoords);
	if (mu_t == 0.0 || E0 == 0.0)
	    r = 0.0;
	else
	    r = S*K0/E0;
	return eqn_params->C2
                        + (eqn_params->Cmu*r*r*r*(1.0-r/4.38))
                        / (1.0 + 0.012*r*r*r);
}

double KE_CARTESIAN::computePointFieldC1_REAL(int* icoords,double S)
{
	int index;
	double C1, r, K0, E0;
	index = d_index(icoords,top_gmax,dim);
	K0 = field->k[index];
        E0 = field->eps[index];
	if (E0 < 10000*MACH_EPS)
	{
	    return 1.0;
	}
	else
	{
	    r = S*K0/E0;
	    r = std::max(r,0.0);
	    return std::max(0.43,r/(r+5));
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

//TODO: Rename and implement Crank Nicholson scheme (2nd order time)
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
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
        double *E = field->eps;
	double *Pk = field->Pk;
	double *mu_t = field->mu_t;
	double sigma_eps = eqn_params->sigma_eps;
	double rho = eqn_params->rho;
        double v[MAXD],v_wall[MAXD];
        double eta;
	double Ut; /*friction velocity*/
	double nu = eqn_params->mu/eqn_params->rho;
	/*For distance y*/
	double y;
	double center[MAXD], t[MAXD], point[MAXD],crds_wall[MAXD];
	double dist;
	double vn;
	double nor[MAXD];
	/*For boundary state*/
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *intfc_state;
	INTERFACE* grid_intfc = front->grid_intfc;
	double coords[MAXD];
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
                comp = top_comp[ic];
                I = ij_to_I[i][j];
                if (comp != sub_comp)
                    continue;
                E0 = E[ic];
		K0 = field->k[ic];
		if (keps_model == REALIZABLE)
		{
		    S = computePointFieldStrain(icoords);
		    rhs = E0 + m_dt*std::max(computePointFieldC1_REAL(icoords,S)*S*E0,0.0);
		}
		else
            rhs = E0 + m_dt*std::max(Pk[ic]*eqn_params->C1*E0/K0,0.0); 

		if (keps_model == RNG)
		{
		    C2 = computePointFieldC2_RNG(icoords);
                    coeff = 1.0 + std::max(C2*E0/K0,0.0)*m_dt;
		}
		else if (keps_model == REALIZABLE)
		{
		    C2 = eqn_params->C2;
		    coeff = 1.0 + std::max(C2*E0/(K0+sqrt(std::max(nu*E0,0.0))),0.0)*m_dt;
		}
		else
		{ 
		    C2 = eqn_params->C2;
                    coeff = 1.0 + std::max(C2*E0/K0,0.0)*m_dt;
		}
                for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }

        /*        
		//set values at points adjacent to wall interface
		//skip these points after settings
		if_adj_pt = NO;
		for (l = 0; l < dim; ++l)
		for (m = 0; m < 2; ++m)
		{
		    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icoords,dir[l][m],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords);

		    //if (fr_crx_grid_seg && (wave_type(hs) == NEUMANN_BOUNDARY ||
              //          wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                //        wave_type(hs) == ELASTIC_BOUNDARY))//
		    if (fr_crx_grid_seg && (wave_type(hs) == NEUMANN_BOUNDARY ||
                        wave_type(hs) == MOVABLE_BODY_BOUNDARY))
		    {
			if_adj_pt = YES;
			for (ll = 0; ll < dim; ll++)
			    crds_wall[ll] = crx_coords[ll];
		    }
		}
		
		if (if_adj_pt == YES)
		{
		    //found adjacent point, use wall function
		    getRectangleCenter(ic,center);
		    dist = distance_between_positions(center,crds_wall,dim);
		    //set a lower bound for dist, since y+ > 11.067
		    if (field->k[ic] > 0.0)
		        dist = std::max(nu*eqn_params->y_p/(pow(eqn_params->Cmu,
				0.25)*pow(field->k[ic],0.25)),dist);
                    //found adjacent point, use wall function
		    if (field->k[ic] > 0.0)
                        rhs = pow(eqn_params->Cmu,0.75)*pow(field->k[ic],1.5)
                              /(0.41*dist);
		    else
			rhs = 0.0;
		    coeff = 1.0;
		}
		else
		{*/
		    D = nu+mu_t[ic]/eqn_params->sigma_eps/rho;
                for (l = 0; l < dim; ++l)
                {
                    lambda = D*m_dt/sqr(top_h[l]);
                    eta = v[l]*m_dt/(top_h[l]); //upwind difference
                    double eta_p = std::max(eta, 0.0);
                    double eta_m = std::min(eta, 0.0);
                        //eta = v[l]*m_dt/(2*top_h[l]);
                    coeff += 2*lambda;

                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
                        I_nb = ij_to_I[ipn[0]][ipn[1]];
            
                        fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icoords,dir[l][m],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords);
                        
                        coeff += ((m == 0) ? eta_p : -eta_m); //upwind
                        if (!fr_crx_grid_seg) 
                        {
                            coeff_nb = -lambda;
                            coeff_nb += (m == 0) ? -eta_p : eta_m;
                                //coeff_nb = -lambda + pow(-1,m+1)*eta;
                            solver.Add_A(I,I_nb,coeff_nb);
                        }
			/*else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                    wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                    wave_type(hs) == ELASTIC_BOUNDARY)*/
			else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                    wave_type(hs) == MOVABLE_BODY_BOUNDARY)
			{

                /*
                //TODO: is this impposible claim legitimate?
				    //printf("decting Neumann Boundary, impossible, check!\n");
                
                //use wall function for friction velocity u_t
                boolean status;
                status = FT_NormalAtGridCrossing(front,icoords,dir[l][m],
                        comp,nor,&hs,crx_coords);

                vn = 0.0;
                double* vel = intfc_state->vel;
                for (int kk = 0; kk < 2; ++kk)
                    vn += (v[kk] - vel[kk])*nor[kk];
                        //vn += vel[kk]*nor[kk];
                
                double v_tan[2];
                for (int kk = 0; kk < 2; ++kk)
                    v_tan[kk] = v[kk] - vn*nor[kk];
                */
                double v_tan[MAXD];
                setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_tan);

                double mag_vtan = Magd(v_tan,dim);

                //friction velocity
                double u_t = std::max(
                        pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                        mag_vtan/eqn_params->y_p);
                    
                //TODO: 0.41 is the karman constant and should not be hardcoded
                E_nb = pow(u_t,4.0)/(0.41*eqn_params->y_p*nu);
			    rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
			}
			else if (wave_type(hs) == DIRICHLET_BOUNDARY)
			{

			    if (boundary_state_function_name(hs) &&
				strcmp(boundary_state_function_name(hs),
                            	"flowThroughBoundaryState") == 0)
			    {
                    //Outlet
                    E_nb = E0;
			    	rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                        //rhs += lambda*E_nb + eta*pow(-1,m)*E_nb;
			    }
			    else
			    {
                    //Inlet
			        E_nb = eqn_params->Cmu
				     *pow(eqn_params->Cbc,1.5)
				     *pow(Mag2d(intfc_state->vel), 3.0)
				     /eqn_params->l0;
			    	rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
			        /*E_nb = eqn_params->Cmu
				     *pow(eqn_params->Cbc,1.5)
				     *pow(sqr(intfc_state->vel[0])
				     +sqr(intfc_state->vel[1]),1.5)
				     /eqn_params->l0;
                        //rhs += lambda*E_nb + eta*pow(-1,m)*E_nb;*/
			    }
			}
			else
            {
                //printf("Unknown boundary condition! \n");
                //clean_up(ERROR);
            }
                    }  /*m*/
		  }  /*l*/
                
        //}  /*if_adj_pt*/

		if (isnan(coeff) || isinf(coeff) || isnan(rhs) || isinf(rhs))
		{
		    printf("Warning: at [%d %d], coeff = %f, rhs = %f\n",i,j,coeff,rhs);
		    clean_up(ERROR);
		    coeff = 1.0;
		    rhs = E0;
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
                ic = d_index3d(i,j,k,top_gmax);
                comp = top_comp[ic];
                I = ijk_to_I[i][j][k];
                if (comp != sub_comp)
                    continue;
                E0 = E[ic];
		K0 = field->k[ic];
		if (keps_model == REALIZABLE)
		{
		    S = computePointFieldStrain(icoords);
		    rhs = E0 + m_dt*computePointFieldC1_REAL(icoords,S)*S*E0;
		}
		else
            rhs = E0 + m_dt*std::max(Pk[ic]*eqn_params->C1*E0/K0,0.0); 

		if (keps_model == RNG)
		{
		    C2 = computePointFieldC2_RNG(icoords);
                    coeff = 1.0 + std::max(C2*E0/K0,0.0)*m_dt;
		}
		else if (keps_model == REALIZABLE)
		{
		    C2 = eqn_params->C2;
		    coeff = 1.0 + std::max(C2*E0/(K0+sqrt(std::max(nu*E0,0.0))),0.0)*m_dt;
		}
		else
		{ 
		    C2 = eqn_params->C2;
                    coeff = 1.0 + std::max(C2*E0/K0,0.0)*m_dt;
		}
                for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
		D = nu+mu_t[ic]/eqn_params->sigma_eps/rho;
        for (l = 0; l < dim; ++l)
        {
                    lambda = D*m_dt/sqr(top_h[l]);
		    //eta = 0.5*v[l]*m_dt/(top_h[l]);
		    eta = v[l]*m_dt/(top_h[l]); //upwind difference
		    double eta_p = std::max(eta, 0.0);
		    double eta_m = std::min(eta, 0.0);
                    coeff += 2.0*lambda;

        for (m = 0; m < 2; ++m)
        {
            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
            icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
            I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];

			fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
				grid_intfc,icoords,dir[l][m],comp,
				(POINTER*)&intfc_state,&hs,crx_coords);

			coeff += eta;
			//coeff += ((m == 0) ? eta_p : -eta_m); //upwind

            if (!fr_crx_grid_seg) 
            {
                coeff_nb = -lambda;
                coeff_nb += (m == 0) ? -eta_p : eta_m;
                solver.Add_A(I,I_nb,coeff_nb);
            }
			/*else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                    wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                    wave_type(hs) == ELASTIC_BOUNDARY)*/
			else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                    wave_type(hs) == MOVABLE_BODY_BOUNDARY)
			{
                /*
                //TODO: should consolidate this with identical
                //      procedure in TKE transport solver.
                        
                //use wall function for friction velocity u_t
                boolean status;
                status = FT_NormalAtGridCrossing(front,icoords,dir[l][m],
                        comp,nor,&hs,crx_coords);

                vn = 0.0;
                double* vel = intfc_state->vel;
                for (int kk = 0; kk < 3; ++kk)
                    vn += (v[kk] - vel[kk])*nor[kk];
                
                double v_tan[3];
                for (int kk = 0; kk < 3; ++kk)
                    v_tan[kk] = v[kk] - vn*nor[kk];

                double u_t = std::max(
                        pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                        Mag3d(v_tan)/eqn_params->y_p);
                */

                double v_tan[MAXD];
                setSlipBoundary(icoords,l,m,comp,hs,intfc_state,field->vel,v_tan);

                double mag_vtan = Magd(v_tan,dim);

                //friction velocity
                double u_t = std::max(
                        pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                        mag_vtan/eqn_params->y_p);
                    
                E_nb = pow(u_t,4.0)/(0.41*eqn_params->y_p*nu);
			    rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
			}
			else if (wave_type(hs) == DIRICHLET_BOUNDARY)
			{

			    if (boundary_state_function_name(hs) &&
				strcmp(boundary_state_function_name(hs),
                            	"flowThroughBoundaryState") == 0)
			    {
                    //Outlet
				    E_nb = E0;
			    	rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
			    }
			    else
			    {
                    //Inlet
			        E_nb = eqn_params->Cmu
				     *pow(eqn_params->Cbc,1.5)
				     *pow(Mag3d(intfc_state->vel), 3.0)
				     /eqn_params->l0;
			    	rhs += lambda * E_nb + ((m == 0) ? eta_p*E_nb : -eta_m*E_nb); 
                }
			}
			else
            {
                //printf("Unknows boundary condition! \n");
                //clean_up(ERROR);
            }
                    }  /*m*/
		  }  /*l*/
                
                solver.Add_A(I,I,coeff);
		dbg_max_rhs = std::max(rhs, dbg_max_rhs);
		dbg_max_aii = std::max(coeff, dbg_max_aii);
                solver.Add_b(I,rhs);
            }
            break;
	}
        stop_clock("set_coefficients");
        start_clock("petsc_solve");
        solver.SetMaxIter(10000);
        //solver.SetMaxIter(500);
        solver.SetTol(1e-8);
        solver.Solve();
	solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);
        if (debugging("PETSc"))
        {
            (void) printf("KE_CARTESIAN::computeAdvectionE: "
                        "num_iter = %d, rel_residual = %g \n",
                        num_iter, rel_residual);
	    printf("max rhs = %f\n", dbg_max_rhs);
	    printf("max aii = %f\n", dbg_max_aii);
        }
	if (rel_residual > 1)
	{
	    (void) printf("Solution diverges!\n");
	    clean_up(ERROR);
	}

        FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
        solver.Get_x(x);
        stop_clock("petsc_solve");

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
    m_dt = dt;

    if (!activated)
    {
        if (front->time < eqn_params->t0)
            return;

        FT_Save(front);
        FT_Draw(front);
        
        applyInitialConditions();
        activateKE();
        
        printf("\n\nTurbulence Model Activated\n\n");
    }

	computeMuTurb();

	if (debugging("trace")) printf("Entering keps_solve()\n");
	start_clock("keps_solve");

    //TODO: can we remove this call since it is called before solve(),
    //      and all memory is statically allocated only on first call.
    //      The only thing it does is (redundently?) set the grid buffers
    //      and assigns field to eqn_params->field ... but those should
    //      already by set -- maybe be a no-op?
	
    setDomain(); //TODO: make sure this doesn't mess anything up
    //if (debugging("keps_solve")) printf("Passing setDomain()\n");

	setComponent();
	//if (debugging("keps_solve")) printf("Passing setComponent()\n");

    //TODO: remove
        //updateGamma();

    //TODO: should compute the S^2 = 0.5*|grad(u) + grad(u)^T|^2 only and
    //      Fk (aka Pk) and Feps = C1*k*S^2 can be computed later
    //
    //computes the production term Pk
    computeSource();
    //if (debugging("keps_solve")) printf("Passing computeSource()\n");


    computeKE();
    
    //TODO: rename this -- misleading
	    //computeAdvection();
	    //if (debugging("keps_solve")) printf("Passing computeAdvection()\n");

    //TODO: compute using previous time step values
	//computeMuTurb();
	//if (debugging("keps_solve")) printf("Passing computeMuTurb()\n");

	//setAdvectionDt();
	//if (debugging("keps_solve")) printf("Passing setAdvectionDt()\n");

	stop_clock("keps_solve");
	if (debugging("trace")) printf("Leaving keps_solve()\n");
}

//TODO: Write CSV file to visualize
static void printField2d(double *var,
		       const char* varname, 
		       int* ic_min, 
		       int* ic_max,
		       int* top_gmax)
{
	FILE* outfile;
	outfile = fopen(varname,"w");
	for (int j = ic_min[1]; j <= ic_max[1]; j++)
	{
	    for (int i = ic_min[0]; i <= ic_max[0]; i++)
	    {
	        int index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%e ",var[index]);
	    }
	    fprintf(outfile,"\n");
	}
	fclose(outfile);	
}

static void printField3d(double *var,
		       const char* varname, 
		       int* ic_min, 
		       int* ic_max,
		       int* top_gmax)
{
	FILE* outfile;
	outfile = fopen(varname,"w");
	int i = (ic_min[0] + ic_max[0])/2;
	for (int k = ic_min[2]; k <= ic_max[2]; k++)
	{
	    for (int j = ic_min[1]; j <= ic_max[1]; j++)
	    {
	        int index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%e ",var[index]);
	    }
	    fprintf(outfile,"\n");
	}
	fclose(outfile);	
}

//TODO: The roles of Pk and S should be switched
double KE_CARTESIAN::computePointFieldStrain(int* icoords)
{
	//S = sqrt(2*Sij*Sij)
	//Pk*rho/mu_t = 0.5*(Uij + Uji)^2 = 2*Sij*Sij
	//Sij = 0.5*(Uij + Uji) and Uij = dUi/dxj
	int index = d_index(icoords, top_gmax, dim);
	return sqrt(field->Pk[index]*eqn_params->rho/field->mu_t[index]);
}

double KE_CARTESIAN::computePointFieldCmu(int* icoords)
{
	int i,j,k,l,m,index;
	char fname[200];
	static int count = 0;
	INTERFACE *grid_intfc = front->grid_intfc;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
        STATE* intfc_state;
	COMPONENT comp;
	double crx_coords[MAXD];
	boolean Adj_Pt = NO;
	boolean fr_crx_grid_seg;
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
	double d_h[2],vel_nb[2],v_tan[MAXD];
	int index_nb,nb;

	index = d_index(icoords,top_gmax,dim);
  	comp = top_comp[index];
	getRectangleCenter(index,center);

 	if (!ifluid_comp(comp)) return 0.0;
	/*compute module of the strain rate tensor*/
	for (l = 0; l < dim; l++)
	for (m = 0; m < dim; m++)
	{		
	    for (nb = 0; nb < 2; nb++)
	    {
		fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icoords,dir[m][nb],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords); 
		d_h[nb] = top_h[m]; 
		if (!fr_crx_grid_seg)
		{
		    index_nb = next_index_in_dir(icoords,dir[m][nb],dim,top_gmax);
		    vel_nb[nb] = vel[l][index_nb];
		}
		/*else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                wave_type(hs) == ELASTIC_BOUNDARY)*/
		else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
		{
		    setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,field->vel,v_tan);
		    vel_nb[nb] = v_tan[l];
		}
		else if (wave_type(hs) == DIRICHLET_BOUNDARY)
		{
		    if (boundary_state_function_name(hs) &&
                        strcmp(boundary_state_function_name(hs),
                        "flowThroughBoundaryState") == 0)
		    {
			vel_nb[nb] = vel[l][index];
		    }
		    else
		    {
			vel_nb[nb] = intfc_state->vel[l];
		    }
		}
	    }
	    S[l][m] = 0.5*(vel_nb[1]- vel_nb[0])/(d_h[1]+d_h[0]);
	    R[l][m] = S[l][m];

	    for (nb = 0; nb < 2; nb++)
    	    {
		fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icoords,dir[l][nb],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords); 
		d_h[nb] = top_h[l]; 
		if (!fr_crx_grid_seg)
		{
		    index_nb = next_index_in_dir(icoords,dir[l][nb],dim,top_gmax);
		    vel_nb[nb] = vel[m][index_nb];
		}
		/*else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                wave_type(hs) == ELASTIC_BOUNDARY)*/
		else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
		{
			setSlipBoundary(icoords,l,nb,comp,hs,
					intfc_state,field->vel,v_tan);
		        vel_nb[nb] = v_tan[m];
		}
		else if (wave_type(hs) == DIRICHLET_BOUNDARY)
		{
			if (boundary_state_function_name(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
			{
			    vel_nb[nb] = vel[m][index];
			}
			else
			{
			    vel_nb[nb] = intfc_state->vel[m];
			}
		}
	    }
	    S[l][m] += 0.5*(vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
	    R[l][m] -= 0.5*(vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
	    if (isnan(S[l][m]))
	    printf("wave_type = %d, vel_nb = [%f %f], d_h = [%f %f]\n",
		    wave_type(hs),vel_nb[0],vel_nb[1],d_h[0],d_h[1]);
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
	    W = 0.0;

	W = W > 0? std::min(W,1.0/sqrt(6.0)) : std::max(W,-1.0/sqrt(6.0));
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
	    Cmu = 1/(A0+As*U*std::max(K[index]/E[index],0.0));
	Cmu = std::min(Cmu,eqn_params->Cmu);
	
	if (isnan(Cmu) || isinf(Cmu))
	{
	    printf("Warning: Cmu = %f, E = %f, K = %f, A0 = %f, As = %f, U = %f\n",
		    Cmu,E[index],K[index],A0,As,U);
	    printf("phi = %f, W = %f\n",phi,W);
	    clean_up(ERROR);
	}
	return Cmu;
}

void KE_CARTESIAN::computeMuTurb()
{
	int i,j,k,l,m,ll,index;
	char fname[200];
	static int count = 0;
	INTERFACE *grid_intfc = front->grid_intfc;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
        POINTER intfc_state;
	COMPONENT comp;
	int icoords[MAXD];
	double coords[MAXD];
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

    bool kill_flag = false;

    int idxKMIN, idxKMAX;
    double crdKMIN[MAXD] = {0}; double crdKMAX[MAXD] = {0};
    double KMIN = HUGE; double KMAX = -HUGE;
    
    int idxEMIN, idxEMAX;
    double crdEMIN[MAXD] = {0}; double crdEMAX[MAXD] = {0};
    double EMIN = HUGE; double EMAX = -HUGE;

    int idxMutMIN, idxMutMAX;
    double crdMutMIN[MAXD] = {0}; double crdMutMAX[MAXD] = {0};
    double MutMIN = HUGE; double MutMAX = -HUGE;

    //TODO: can probably update gammaK here

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
            if (!ifluid_comp(comp))
            {
                field->mu_t[index] = 0.0;
                continue;
            }

		    if (keps_model == REALIZABLE)
		    {
                Cmu = computePointFieldCmu(icoords);
                field->Cmu[index] = Cmu;
		    }
		    else
			    Cmu = eqn_params->Cmu;

            /*
            double limited_mix_length = lmax;
            if (Cmu*pow(field->k[index],1.5) < field->eps[index]*lmax)
            {
                limited_mix_length =
                    Cmu*pow(field->k[index],1.5)/field->eps[index];
            }
            
            //TODO: make 0.001 a variable fraction
            double nu_min = 0.001*nu;
            double nu_t = std::max(nu_min,
                    limited_mix_length*sqrt(field->k[index]));
            
            field->mu_t[index] = nu_t*rho;
            */

            if (field->k[index] < KMIN)
            {
                KMIN = field->k[index];
                idxKMIN = index;
            }
            if (field->k[index] > KMAX)
            {
                KMAX = field->k[index];
                idxKMAX = index;
            }
            if (field->eps[index] < EMIN)
            {
                EMIN = field->eps[index];
                idxEMIN = index;
            }
            if (field->eps[index] > EMAX)
            {
                EMAX = field->eps[index];
                idxEMAX = index;
            }

            double nu_t = std::max(Cmu*sqr(field->k[index])/field->eps[index],0.0001*nu);
            field->mu_t[index] = nu_t*rho;

            if (field->mu_t[index] < MutMIN)
            {
                MutMIN = field->mu_t[index];
                idxMutMIN = index;
            }
            if (field->mu_t[index] > MutMAX)
            {
                MutMAX = field->mu_t[index];
                idxMutMAX = index;
            }


            /*
            if (fabs(field->eps[index]) > MACH_EPS)
            {
                field->mu_t[index] =
                    Cmu*sqr(field->k[index])/field->eps[index]*eqn_params->rho;
            }
            else
            {
                field->mu_t[index] = 0.0001*eqn_params->mu;
            }

            
            field->mu_t[index] = std::max(field->mu_t[index],0.0001*eqn_params->mu);
            */
            
            if (isnan(field->mu_t[index]) || isinf(field->mu_t[index]))
		    {
			    printf("Warning: eddy viscosity is nan or inf\n \
                        \t mu_t=%f, Cmu=%f, k=%f, eps=%f\n",
                        field->mu_t[index],Cmu,field->k[index],field->eps[index]);
                
                getRectangleCenter(index,coords);
                printf("icoords = %d %d, index = %d, coords = %f %f\n",
                        i,j,index,coords[0],coords[1]);
                
                kill_flag = true;
                    //clean_up(EXIT_FAILURE);
		    }

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
            if (!ifluid_comp(comp))
            {
                field->mu_t[index] = 0.0;
                continue;
            }
		    
            if (keps_model == REALIZABLE)
		    {
                Cmu = computePointFieldCmu(icoords);
                field->Cmu[index] = Cmu;
		    }
		    else
                Cmu = eqn_params->Cmu;
	    	
            /*
            double limited_mix_length = lmax;
            if (Cmu*pow(field->k[index],1.5) < field->eps[index]*lmax)
            {
                limited_mix_length =
                    Cmu*pow(field->k[index],1.5)/field->eps[index];
            }
            
            //TODO: make 0.001 a variable fraction
            double nu_min = 0.001*nu;
            double nu_t = std::max(nu_min,
                    limited_mix_length*sqrt(field->k[index]));

            field->mu_t[index] = nu_t*rho;
            */

            if (field->k[index] < KMIN)
            {
                KMIN = field->k[index];
                idxKMIN = index;
            }
            if (field->k[index] > KMAX)
            {
                KMAX = field->k[index];
                idxKMAX = index;
            }
            if (field->eps[index] < EMIN)
            {
                EMIN = field->eps[index];
                idxEMIN = index;
            }
            if (field->eps[index] > EMAX)
            {
                EMAX = field->eps[index];
                idxEMAX = index;
            }

            double nu_t = std::max(Cmu*sqr(field->k[index])/field->eps[index],0.0001*nu);
            field->mu_t[index] = nu_t*rho;

            if (field->mu_t[index] < MutMIN)
            {
                MutMIN = field->mu_t[index];
                idxMutMIN = index;
            }
            if (field->mu_t[index] > MutMAX)
            {
                MutMAX = field->mu_t[index];
                idxMutMAX = index;
            }

            /*
            if (fabs(field->eps[index]) > MACH_EPS)
            {
                field->mu_t[index] =
                    Cmu*sqr(field->k[index])/field->eps[index]*eqn_params->rho;
            }
            else
            {
                field->mu_t[index] = 0.0001*eqn_params->mu;
            }


            field->mu_t[index] = std::max(field->mu_t[index],0.0001*eqn_params->mu);
            */

            if (isnan(field->mu_t[index]) || isinf(field->mu_t[index]))
		    {
			    printf("Warning: eddy viscosity is nan or inf\n \
                        \t mu_t=%f,Cmu=%f, k=%f, eps=%f\n",
                        field->mu_t[index],Cmu,field->k[index],field->eps[index]);

                getRectangleCenter(index,coords);
                printf("icoords = %d %d %d, index = %d, coords = %f %f %f\n",
                        i,j,k,index,coords[0],coords[1],coords[2]);
                
                kill_flag = true;
                    //clean_up(EXIT_FAILURE);
		    }

		}
        break;
   
        default:
            printf("Invalid Dimension\n"); LOC();
            clean_up(EXIT_FAILURE);
	}

    getRectangleCenter(idxKMIN,crdKMIN);
    getRectangleCenter(idxKMAX,crdKMAX);
    getRectangleCenter(idxEMIN,crdEMIN);
    getRectangleCenter(idxEMAX,crdEMAX);
    getRectangleCenter(idxMutMIN,crdMutMIN);
    getRectangleCenter(idxMutMAX,crdMutMAX);

    printf("\n");
    printf("KMIN = %f |  coords = %f %f %f\n",KMIN,crdKMIN[0],crdKMIN[1],crdKMIN[2]);
    printf("KMAX = %f |  coords = %f %f %f\n",KMAX,crdKMAX[0],crdKMAX[1],crdKMAX[2]);
    printf("EMIN = %f |  coords = %f %f %f\n",EMIN,crdEMIN[0],crdEMIN[1],crdEMIN[2]);
    printf("EMAX = %f |  coords = %f %f %f\n",EMAX,crdEMAX[0],crdEMAX[1],crdEMAX[2]);
    printf("MutMIN = %f |  coords = %f %f %f\n",MutMIN,crdMutMIN[0],crdMutMIN[1],crdMutMIN[2]);
    printf("MutMAX = %f |  coords = %f %f %f\n",MutMAX,crdMutMAX[0],crdMutMAX[1],crdMutMAX[2]);

        //printf("KMIN = %f, KMAX = %f\n",KMIN,KMAX);
        //printf("EMIN = %f, EMAX = %f\n",EMIN,EMAX);

    if (kill_flag)
    {
        LOC();
        clean_up(EXIT_FAILURE);
    }

    //TODO: can move some of these parallelizations into the
    //      compute K and Eps functions since we now update
    //      mu_t with the previous time step values.
    //      parallel exchange of mu_t is of course still required.
	
    FT_ParallelExchGridArrayBuffer(field->k,front,NULL);
	FT_ParallelExchGridArrayBuffer(field->eps,front,NULL);
	FT_ParallelExchGridArrayBuffer(field->mu_t,front,NULL);
	if (keps_model == REALIZABLE)
	    FT_ParallelExchGridArrayBuffer(field->Cmu,front,NULL);


	int ic_min[MAXD] = {imin,jmin,kmin};
    int ic_max[MAXD] = {imax,jmax,kmax};
	    //ic_min[0] = imin; ic_min[1] = jmin; ic_min[2] = kmin;
	    //ic_max[0] = imax; ic_max[1] = jmax; ic_max[2] = kmax;

    //TODO: Add hdf/vtk movie variables for visualization
	if (dim == 2)
	{
	    sprintf(fname,"%s/K_field",OutName(front));
	    printField2d(field->k,fname,ic_min,ic_max,top_gmax);
	    sprintf(fname,"%s/E_field",OutName(front));
	    printField2d(field->eps,fname,ic_min,ic_max,top_gmax);
	}
	else if (dim == 3)
	{
        sprintf(fname,"%s/K_field",OutName(front));
        printField3d(field->k,fname,ic_min,ic_max,top_gmax);
        sprintf(fname,"%s/E_field",OutName(front));
        printField3d(field->eps,fname,ic_min,ic_max,top_gmax);
	    if (keps_model == REALIZABLE)
	    {
            sprintf(fname,"%s/Cmu",OutName(front));
            printField3d(field->Cmu,fname,ic_min,ic_max,top_gmax);
	    }
	}
}
	
//TODO: explicit scheme unfeasable
void KE_CARTESIAN::setAdvectionDt()
{
	double D, Dl, Ds;
	double mu_max = -HUGE;

	for (int i = 0; i < comp_size; i++)
		mu_max = std::max(field->mu_t[i],mu_max);

	Dl = mu_max/eqn_params->sigma_k + eqn_params->mu;
	Ds = mu_max/eqn_params->sigma_eps + eqn_params->mu;
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
		exit(0);
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

    //TODO: Is there a bug with ilower and iupper vals? See comment at top of file.
	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];

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
	    FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(KE_FIELD));

	switch (dim)
	{
	case 1:
        //TODO: remove 1d -- makes no sense
        if (first)
        {
            comp_size = top_gmax[0]+1;
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
            FT_VectorMemoryAlloc((POINTER*)&Karray,comp_size,FLOAT);
            FT_VectorMemoryAlloc((POINTER*)&Earray,comp_size,FLOAT);
            FT_VectorMemoryAlloc((POINTER*)&field->Pk,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->gamma,comp_size,FLOAT);
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
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);//TODO:remove
	    	FT_VectorMemoryAlloc((POINTER*)&Karray,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&Earray,comp_size,FLOAT);
            FT_VectorMemoryAlloc((POINTER*)&field->Pk,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->gamma,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->k,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->eps,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->temp,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->mu_t,comp_size,FLOAT);
		    
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
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);//TODO:remove
	    	FT_VectorMemoryAlloc((POINTER*)&Karray,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&Earray,comp_size,FLOAT);
            FT_VectorMemoryAlloc((POINTER*)&field->Pk,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->gamma,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->k,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->eps,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->mu_t,comp_size,FLOAT);
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
	FILE* infile;
    char* inname = InName(front);
	infile = fopen(inname,"r");
	char string[100];

	switch (dim)
	{
	case 2:
	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,"mu_t",0,
				field->mu_t,getStateTemp,0,0);
	    break;
	case 3:
        if (CursorAfterStringOpt(infile,
                    "Type y to make eddy viscosity field movie:"))
        {
            fscanf(infile,"%s",string);
            (void)printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkScalarMovieVariable(front,"EDDYVISCOSITY",field->mu_t);
        }

        if (CursorAfterStringOpt(infile,
                    "Type y to make turbulent kinetic energy field movie:"))
        {
            fscanf(infile,"%s",string);
            (void)printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkScalarMovieVariable(front,"TKE",field->k);
        }

        if (CursorAfterStringOpt(infile,
                    "Type y to make turbulent kinetic energy dissipation field movie:"))
        {
            fscanf(infile,"%s",string);
            (void)printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkScalarMovieVariable(front,"EPS",field->eps);
        }
	    break;
	}

    fclose(infile);
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
    
    if (status == NO) return NO_PDE_BOUNDARY;

    if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE)
        return NO_PDE_BOUNDARY;
    else if (wave_type(*hs) == DIRICHLET_BOUNDARY)
        return DIRICHLET_PDE_BOUNDARY;
	else if (wave_type(*hs) == GROWING_BODY_BOUNDARY)
        return DIRICHLET_PDE_BOUNDARY;
	else if (wave_type(*hs) == NEUMANN_BOUNDARY ||
             wave_type(*hs) == ELASTIC_BOUNDARY)
        return NEUMANN_PDE_BOUNDARY;
}       /* find_state_at_crossing */

static int next_index_in_dir(int* icoords,GRID_DIRECTION dir,int dim,int* top_gmax)
{
	int index,i;
	int icrds[MAXD];
	for (i = 0; i < dim; i++)
	    icrds[i] = icoords[i];
        switch (dir)
        {
        case WEST:
            icrds[0] -= 1;
            break;
        case EAST:
            icrds[0] += 1;
            break;
        case SOUTH:
            icrds[1] -= 1;
            break;
        case NORTH:
            icrds[1] += 1;
            break;
        case LOWER:
            icrds[2] -= 1;
            break;
        case UPPER:
            icrds[2] += 1;
        }
	index = d_index(icrds,top_gmax,dim);
	return index;
}

//TODO?
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

    //TODO: Need to generate ghost value for k so that
    //      the gradient does not change in direction normal
    //      to the wall. 
    
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

//TODO: make a global utility function
static std::string dir2String(GRID_DIRECTION dir)
{
    switch (dir)
    {
        case EAST:
           return "EAST";
           break;
        case WEST:
           return "WEST";
           break;
        case NORTH:
           return "NORTH";
           break;
        case SOUTH:
           return "SOUTH";
           break;
        case UPPER:
           return "UPPER";
           break;
        case LOWER:
           return "LOWER";
           break;
        default:
           printf("not a known GRID_DIRECTION\n");
           clean_up(ERROR);
    }
}

void KE_CARTESIAN::setSlipBoundary(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double** vel,
	double* v_tan)
{
	int index;
    int ic[MAXD];
    double  coords[MAXD],coords_ref[MAXD],crx_coords[MAXD],coords_ghost[MAXD];

    double nor[MAXD];
    double v[MAXD];
    double v_tmp[MAXD];
    double vn;

    GRID_DIRECTION  ldir[3] = {WEST,SOUTH,LOWER};
    GRID_DIRECTION  rdir[3] = {EAST,NORTH,UPPER};
    GRID_DIRECTION  dir;
    double  vel_intfc[MAXD];

	index = d_index(icoords,top_gmax,dim);
	for (int i = 0; i < dim; ++i)
    {
        vel_intfc[i] = (*getStateVel[i])(state);
        coords[i] = top_L[i] + icoords[i]*top_h[i];
        ic[i] = icoords[i];
    }
	
    dir = (nb == 0) ? ldir[idir] : rdir[idir];
    ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
    
    //TODO: really don't need this call since we already know this,
    //      but should be passing in the the values computed by the calling
    //      function
    boolean status;
    status = FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);
    if (status == NO) return;

    //ghost point coords
    //TODO: can use getRectangleCenter()???
    for (int j = 0; j < dim; ++j)
    {
        coords_ghost[j] = top_L[j] + ic[j]*top_h[j];
        coords_ref[j] = coords_ghost[j];
    }

    /* Reflect ghost point through intfc-mirror at crossing */
    //first reflect across the grid line containing the intfc crossing
    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];
    vn = 0.0;
    
    for (int j = 0; j < dim; ++j)
    {
        v[j] = coords_ref[j] - crx_coords[j];
        vn += v[j]*nor[j];
    }

    //reflect v across the line containing the normal vector
    for (int j = 0; j < dim; ++j)
        v[j] = 2.0*vn*nor[j] - v[j];
  
    //desired reflected point
    for (int j = 0; j < dim; ++j)
        coords_ref[j] = crx_coords[j] + v[j];

    /* Interpolate the state at the reflected point */
    for (int j = 0; j < dim; ++j)
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,vel[j],
        	getStateVel[j],&v_tmp[j],&vel[j][index]);
	
    /* zero the relative normal velocity leaving the tangential component unchanged */
    vn = 0.0;
    for (int j = 0; j < dim; j++)
    {
        v[j] = v_tmp[j] - vel_intfc[j];
        vn += v[j]*nor[j];
    }

    //Enforce u dot nor = 0, NOT du/dn = 0
    for (int j = 0; j < dim; ++j)
    {
        v_tan[j] = v_tmp[j] - vn*nor[j];
            //v_tan[j] = v_tmp[j] - 2.0*vn*nor[j]; //TODO: NO
    }
}

//Pk wall boundary condition
double KE_CARTESIAN::computeWallPk(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER intfc_state,
	double* v_tan)
{
	int index = d_index(icoords,top_gmax,dim);

    double u_t = std::max(
            pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[index],0.0)),
                            Magd(v_tan,dim)/eqn_params->y_p);

    double nu_t = field->mu_t[index]/eqn_params->rho;
    double Pk_Wall = pow(u_t,3.0)*Magd(v_tan,dim)/(nu_t*eqn_params->y_p);
    return Pk_Wall;
}

//TODO: should probably just compute S^2 for reusability, and Pk can be computed later
void KE_CARTESIAN::computeSource()
{
	int i,j,k,l,m,nb,index;
	int index_nb,icrds[MAXD];
	double **vel = field->vel;
	double vel_nb[2],d_h[2];
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	
	double *Pk = field->Pk;
	double *mu_t = field->mu_t;
	double rho = eqn_params->rho;

	/*find crx*/
	HYPER_SURF *hs;
	STATE *intfc_state;
	INTERFACE* grid_intfc = front->grid_intfc;
	int fr_crx_grid_seg;
	int comp;
	double crx_coords[MAXD],center[MAXD];

    //TODO: make J an array and store values for later computation

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

        //TODO: factor out this common l,m,nb loop into function
        //      that takes the same args plus the icrds array.
        
        double J = 0.0;

        /*compute modulus (squared frobenius matrix norm) of the strain rate tensor*/
        for (l = 0; l < dim; l++)
        {
        for (m = 0; m < dim; m++)
        {
            //l components in the m direction
			for (nb = 0; nb < 2; nb++)
			{
			    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icrds,dir[m][nb],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords); 
			    
                d_h[nb] = top_h[m]; 

			    if (!fr_crx_grid_seg)
			    {
			        index_nb = next_index_in_dir(icrds,dir[m][nb],dim,top_gmax);
			        vel_nb[nb] = vel[l][index_nb];
			    }
			    /*else if(fr_crx_grid_seg && (wave_type(hs) == NEUMANN_BOUNDARY ||
                            wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                            wave_type(hs) == ELASTIC_BOUNDARY))*/
			    else if(fr_crx_grid_seg && (wave_type(hs) == NEUMANN_BOUNDARY ||
                            wave_type(hs) == MOVABLE_BODY_BOUNDARY))
			    {
                    double v_tan[MAXD] = {0.0};
                    setSlipBoundary(icrds,m,nb,comp,hs,intfc_state,field->vel,v_tan);
                    vel_nb[nb] = v_tan[l];
                    //TODO: should field->vel[][index] get set to v_tan (v_slip)
                    
                    /*
                    //TODO: may not be close enough to wall to apply Pk_Wall, 
                    //      and should instead compute using v_tan as before.

                    //Pk wall boundary condition
                    Pk_Wall = computeWallPk(icrds,m,nb,comp,hs,intfc_state,v_tan);
                    */
			    }
			    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
			    {
				    if (boundary_state_function_name(hs) &&
                        strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                    {
                        vel_nb[nb] = vel[l][index];
                    }
                    else
                    {
                        vel_nb[nb] = intfc_state->vel[l];
                    }
                    d_h[nb] = distance_between_positions(center,crx_coords,dim);
                }
			
            }//nb

                
            double S = (vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);

            //m components in the l direction
            for (nb = 0; nb < 2; nb++)
			{
			    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icrds,dir[l][nb],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords); 
			    d_h[nb] = top_h[l]; 
			    if (!fr_crx_grid_seg)
			    {
			        index_nb = 
				next_index_in_dir(icrds,dir[l][nb],dim,top_gmax);
			        vel_nb[nb] = vel[m][index_nb];
			    }
			    /*else if(fr_crx_grid_seg && (wave_type(hs) == NEUMANN_BOUNDARY ||
                            wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                            wave_type(hs) == ELASTIC_BOUNDARY))*/
			    else if(fr_crx_grid_seg && (wave_type(hs) == NEUMANN_BOUNDARY ||
                            wave_type(hs) == MOVABLE_BODY_BOUNDARY))
			    {
                    double v_tan[MAXD] = {0.0};
                    setSlipBoundary(icrds,l,nb,comp,hs,intfc_state,field->vel,v_tan);
                    vel_nb[nb] = v_tan[m];
                    
                    /*
                    //TODO: may not be close enough to wall to apply Pk_Wall, 
                    //      and should instead compute using v_tan as before.
                    
                    //Pk wall boundary condition
                    Pk_Wall = computeWallPk(icrds,m,nb,comp,hs,intfc_state,v_tan);
                    */
			    }
			    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
			    {
                    if (boundary_state_function_name(hs) &&
                                    strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                    {
                        vel_nb[nb] = vel[m][index];
                    }
                    else
                    {
                        vel_nb[nb] = intfc_state->vel[m];
                    }
    
                    d_h[nb] = distance_between_positions(center,crx_coords,dim);
                }

			}//nb
            
           
            S += (vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
            J += S*S;
            
        }//m


        }//l

        double nu_t = mu_t[index]/rho;
        Pk[index] = 0.5*nu_t*J; 

    
    }//j,i
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

        //TODO: factor out this common l,m,nb loop into function
        //      that takes the same args plus the icrds array.
        
        double J = 0.0;
        
        /*compute module of the strain rate tensor*/
        for (l = 0; l < dim; l++)
        {
        for (m = 0; m < dim; m++)
        {		
            //l components in the m direction
			for (nb = 0; nb < 2; nb++)
			{
			    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icrds,dir[m][nb],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords); 
			    d_h[nb] = top_h[m]; 
			    if (!fr_crx_grid_seg)
			    {
			        index_nb = 
                        next_index_in_dir(icrds,dir[m][nb],dim,top_gmax);
			        vel_nb[nb] = vel[l][index_nb];
			    }
			    /*else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                                wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                                wave_type(hs) == ELASTIC_BOUNDARY)*/
			    else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
			    {
                    double v_tan[MAXD] = {0.0};
                    setSlipBoundary(icrds,m,nb,comp,hs,intfc_state,field->vel,v_tan);
                    vel_nb[nb] = v_tan[l];
                    
                    /*
                    //TODO: may not be close enough to wall to apply Pk_Wall, 
                    //      and should instead compute using v_tan as before.
                    
                    //Pk wall boundary condition
                    Pk_Wall = computeWallPk(icrds,m,nb,comp,hs,intfc_state,v_tan);
                    */
			    }
			    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
			    {
                    if (boundary_state_function_name(hs) &&
                                    strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                    {
                        vel_nb[nb] = vel[l][index];
                    }
                    else
                    {
                        vel_nb[nb] = intfc_state->vel[l];
                    }
			    
                    d_h[nb] = distance_between_positions(center,crx_coords,dim);
			    }

			}//nb



            double S = (vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);

            //m components in the l direction
			for (nb = 0; nb < 2; nb++)
			{
			    fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                                grid_intfc,icrds,dir[l][nb],comp,
                                (POINTER*)&intfc_state,&hs,crx_coords); 
			    d_h[nb] = top_h[l]; 
			    if (!fr_crx_grid_seg)
			    {
			        index_nb = next_index_in_dir(icrds,dir[l][nb],dim,top_gmax);
			        vel_nb[nb] = vel[m][index_nb];
			    }
			    /*else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                        wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                        wave_type(hs) == ELASTIC_BOUNDARY)*/
			    else if(wave_type(hs) == NEUMANN_BOUNDARY ||
                        wave_type(hs) == MOVABLE_BODY_BOUNDARY)
			    {
                    double v_tan[MAXD] = {0.0};
                    setSlipBoundary(icrds,l,nb,comp,hs,intfc_state,field->vel,v_tan);
                    vel_nb[nb] = v_tan[m];
                    
                    /*
                    //TODO: may not be close enough to wall to apply Pk_Wall, 
                    //      and should instead compute using v_tan as before.
                    
                    //Pk wall boundary condition
                    Pk_Wall = computeWallPk(icrds,m,nb,comp,hs,intfc_state,v_tan);
                    */
			    }
			    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
			    {
                    if (boundary_state_function_name(hs) &&
                                    strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                    {
                        vel_nb[nb] = vel[m][index];
                    }
                    else
                    {
                        vel_nb[nb] = intfc_state->vel[m];
                    }
    
                    d_h[nb] = distance_between_positions(center,crx_coords,dim);
                }
            
            }//nb


            S += (vel_nb[1] - vel_nb[0])/(d_h[1]+d_h[0]);
            J += S*S;
            
        }//m
            

        }//l
        
        double nu_t = mu_t[index]/rho;
        Pk[index] = 0.5*nu_t*J; 


    }//k,j,i
    break;

    default: 
		printf("In computeSource(), Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}

    //TODO: need to scatter the array???
	return;
}

/*read k-epsilon parameters*/
void KE_CARTESIAN::read_params(
	char *inname,
	KE_PARAMS *eqn_params)
{
	char string[100];
	FILE* infile;
	infile = fopen(inname,"r");

	/*default parameter*/
	eqn_params->sigma_k = 1.0;
	eqn_params->sigma_eps = 1.3;
	eqn_params->Cmu = 0.09;
	eqn_params->C1 = 1.44;
	eqn_params->C2 = 1.92;
	eqn_params->Cbc = 0.01; /*typicaly 0.003 ~ 0.01*/
	eqn_params->rho = 1.0;
	eqn_params->B = 5.2; /*5.2 for smooth wall*/
	eqn_params->E = 9.0; /*9.0 for smooth wall*/
	eqn_params->y_p = 30.0;
	eqn_params->t0 = 0.0;
	keps_model = STANDARD;
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
	fscanf(infile,"%lf",&eqn_params->sigma_k);
	(void) printf("%f\n",eqn_params->sigma_k);

	CursorAfterStringOpt(infile,"Enter turbulent Prandtl number for epsilon:");
	fscanf(infile,"%lf",&eqn_params->sigma_eps);
	(void) printf("%f\n",eqn_params->sigma_eps);

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

	CursorAfterString(infile,"Enter l0:");
	fscanf(infile,"%lf",&eqn_params->l0);
	(void) printf("%f\n",eqn_params->l0);

	CursorAfterStringOpt(infile,"Enter mu0:");
	fscanf(infile,"%lf",&eqn_params->mu0);
	(void) printf("%f\n",eqn_params->mu0);

    CursorAfterStringOpt(infile,"Enter B:");
    fscanf(infile,"%lf",&eqn_params->B);
    (void) printf("%f\n",eqn_params->B);

    CursorAfterStringOpt(infile,"Enter E:");
    fscanf(infile,"%lf",&eqn_params->E);
    (void) printf("%f\n",eqn_params->E);

	CursorAfterStringOpt(infile,"Enter delta:");
	fscanf(infile,"%lf",&eqn_params->delta);
	(void) printf("%f\n",eqn_params->delta);

	CursorAfterStringOpt(infile,"Enter y+:");
	fscanf(infile,"%lf",&eqn_params->y_p);
	(void) printf("%f\n",eqn_params->y_p);

    CursorAfterStringOpt(infile,"Enter time to active turbulence model:");
    fscanf(infile,"%lf",&eqn_params->t0);
    (void) printf("%f\n",eqn_params->t0);
	fclose(infile);

    //initMovieVariables();
}
