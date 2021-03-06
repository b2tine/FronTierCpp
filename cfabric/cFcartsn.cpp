/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

/*******************************************************************
 * 		G_CARTESIAN.c
 *******************************************************************/

#include "cFluid.h"


struct ToFill
{
    int icoords[3];
};

EXPORT void tecplot_interface_states(const char*,INTERFACE*);

static double (*getStateMom[MAXD])(Locstate) = {getStateXmom,
                                                getStateYmom,
                                                getStateZmom};

static void printInputStencil(SWEEP,int);


void G_CARTESIAN::initMesh()
{
	int i,j,k, index;
	double coords[2];
	int num_cells;

	// init cell_center
	L_RECTANGLE       rectangle;

	if (debugging("trace"))
	    (void) printf("Entering g_cartesian.initMesh()\n");
	
    //TODO: why was this done?
    /*TMP*/
	min_dens = 0.0001;
	min_pres = 0.0001;
	
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
	    	coords[0] = top_L[0] + top_h[0]*i;
		index = d_index1d(i,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
		index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
	    	coords[2] = top_L[2] + top_h[2]*k;
		index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}
	
	setComponent();
	FT_FreeGridIntfc(front);
	if (debugging("trace"))
	    (void) printf("Leaving g_cartesian.initMesh()\n");
}

//Also assignes the ghost values for use in solver where appropriate
void G_CARTESIAN::setComponent()
{
	int		i,j, ind;
	double 		*coords;
	int 		*icoords;
	COMPONENT 	old_comp,new_comp;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	static STATE 	*state = NULL;
	double		*dens = field.dens;
	double		*engy = field.engy;
	double		*pres = field.pres;
	double		**momn = field.momn;
	int		size = (int)cell_center.size();
	
	// cell center components
	if(state == NULL)
	    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));

	for (i = 0; i < size; i++)
	{
	    icoords = cell_center[i].icoords;
	    coords = cell_center[i].m_coords;
	    old_comp = cell_center[i].comp;
	    new_comp = top_comp[i];
	    if (eqn_params->tracked && cell_center[i].comp != -1 &&
		cell_center[i].comp != top_comp[i] && gas_comp(new_comp))
	    {
		if (!FrontNearestIntfcState(front,coords,new_comp,
				(POINTER)state))
		{
		    (void) printf("In setComponent()\n");
		    (void) printf("FrontNearestIntfcState() failed\n");
		    (void) printf("old_comp = %d new_comp = %d\n",
					old_comp,new_comp);
		    clean_up(ERROR);
		}

		//GFM
		state->dim = dim;
		state->eos = &eqn_params->eos[new_comp];

		if (gas_comp(old_comp) && gas_comp(new_comp))
		{
		    if(new_comp == GAS_COMP1)
			ind = 0;
		    else
			ind = 1;

            if (Gdens[ind][i] != 0.0) // Not unset
            {
                state->dens = Gdens[ind][i];
                state->pres = Gpres[ind][i];
                for(j = 0; j < dim; ++j)
                    state->momn[j] = Gvel[ind][j][i]*Gdens[ind][i];
                state->engy = EosEnergy(state);
            }
		}

		    dens[i] = state->dens;
            pres[i] = state->pres;
            engy[i] = state->engy;
            for (j = 0; j < dim; ++j)
                momn[j][i] = state->momn[j];
	    }

	    cell_center[i].comp = top_comp[i];
	}
}	/* end setComponent() */

void G_CARTESIAN::setInitialIntfc(
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname)
{
	dim = front->rect_grid->dim;
	eqn_params = (EQN_PARAMS*)front->extra1;
	switch (eqn_params->prob_type)
	{
	case TWO_FLUID_RT:
	case TWO_FLUID_RM:
	    initSinePertIntfc(level_func_pack,inname);
	    break;
	case TWO_FLUID_RM_RAND:
	    initRandPertIntfc(level_func_pack,inname);
	    break;
	case TWO_FLUID_BUBBLE:
	case FLUID_SOLID_CIRCLE:
	    initCirclePlaneIntfc(level_func_pack,inname);
	    break;
	case IMPLOSION:
	    initImplosionIntfc(level_func_pack,inname);
	    break;
	case MT_FUSION:
	    initMTFusionIntfc(level_func_pack,inname);
	    break;
	case PROJECTILE:
	    initProjectileIntfc(level_func_pack,inname);
	    break;
	case FLUID_SOLID_RECT:
	    initRectPlaneIntfc(level_func_pack,inname);
	    break;
	case FLUID_SOLID_TRIANGLE:
	    initTrianglePlaneIntfc(level_func_pack,inname);
	    break;
	case FLUID_SOLID_CYLINDER:
             initCylinderPlaneIntfc(level_func_pack,inname);
             break;
	case RIEMANN_PROB:
	case ONED_BLAST:
	case ONED_SSINE:
	case ONED_ASINE:
	    initRiemannProb(level_func_pack,inname);
	    break;
	case OBLIQUE_SHOCK_REFLECT:
	    initObliqueIntfc(level_func_pack,inname);
	    break;
	default:
	    (void) printf("Problem type not implemented, code needed!\n");
	    clean_up(ERROR);
	}
}	/* end setInitialIntfc */

/*
void G_CARTESIAN::setProbParams(char *inname)
{
	dim = front->rect_grid->dim;
	eqn_params = (EQN_PARAMS*)front->extra1;
	switch (eqn_params->prob_type)
	{
	case TWO_FLUID_RT:
	    setRayleiTaylorParams(inname);
	    break;
	case TWO_FLUID_RM:
	case TWO_FLUID_RM_RAND:
	    setRichtmyerMeshkovParams(inname);
	    break;
	case TWO_FLUID_BUBBLE:
	    setBubbleParams(inname);
	    break;
	case IMPLOSION:
	    setImplosionParams(inname);
	    break;
	case MT_FUSION:
	    setMTFusionParams(inname);
	    break;
	case PROJECTILE:
	case FLUID_SOLID_CIRCLE:
	case FLUID_SOLID_RECT:
	case FLUID_SOLID_TRIANGLE:
	case FLUID_SOLID_CYLINDER:
	    setProjectileParams(inname);
	    break;
	case RIEMANN_PROB:
	    setRiemProbParams(inname);
	    break;
	case ONED_BLAST:
	case ONED_SSINE:
	case ONED_ASINE:
	    setOnedParams(inname);
	    break;
	case OBLIQUE_SHOCK_REFLECT:
	    setRichtmyerMeshkovParams(inname);
	    break;
	default:
	    printf("In setProbParams(), unknown problem type!\n");
	    clean_up(ERROR);
	}
}*/	/* end setProbParams */

/*
void G_CARTESIAN::setInitialStates()
{
	switch (eqn_params->prob_type)
	{
	case TWO_FLUID_RT:
	    initRayleiTaylorStates();
	    break;
	case TWO_FLUID_RM:
	case TWO_FLUID_RM_RAND:
	    initRichtmyerMeshkovStates();
	    break;
	case TWO_FLUID_BUBBLE:
	    initBubbleStates();
	    break;
	case IMPLOSION:
	    initImplosionStates();
	    break;
	case MT_FUSION:
	    initMTFusionStates();
	    break;
	case PROJECTILE:
        case FLUID_SOLID_CIRCLE:
        case FLUID_SOLID_RECT:
        case FLUID_SOLID_TRIANGLE:
        case FLUID_SOLID_CYLINDER:
	    initProjectileStates();
	    break;
	case RIEMANN_PROB:
	    initRiemProbStates();
	    break;
	case ONED_BLAST:
	    initBlastWaveStates();
	    break;
	case ONED_SSINE:
	    initShockSineWaveStates();
	    break;
	case ONED_ASINE:
	    initAccuracySineWaveStates();
	    break;
	case OBLIQUE_SHOCK_REFLECT:
	    initRichtmyerMeshkovStates();
	    break;
	default:
	    (void) printf("In setInitialStates(), case not implemented!\n");
	    clean_up(ERROR);
	}
	copyMeshStates();
}*/	/* end setInitialStates */

void G_CARTESIAN::computeConvectiveFlux()
{
    nrad = 3;
    int order = 1;

	static SWEEP *st_field;
	static FSWEEP *st_flux;
	
	if (st_flux == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&st_field,order,sizeof(SWEEP));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux,order,sizeof(FSWEEP));
	    
        for (int i = 0; i < order; ++i)
	    {
	    	allocMeshVst(&st_field[i]);
	    	allocMeshFlux(&st_flux[i]);
	    }
    }

    double delta_t = m_dt;
	copyToMeshVst(&st_field[0]);
	computeMeshFlux(st_field[0],&st_flux[0],delta_t);
        //addMeshFluxToVst(&st_field[0],st_flux[0],1.0);
	    //copyFromMeshVst(st_field[0]);
    cFlux = &st_flux[0];
}

/*
void G_CARTESIAN::computeDiffusion()
{
    printf("computeDiffusion() not yet implemented\n");
    clean_up(EXIT_FAILURE);

    //TODO: cFlux goes to RHS vector.
    //      below is iFcartsn3d.cpp implementation for template.

    COMPONENT comp;
    int index,index_nb[6],size;
    int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,rhs,U_nb[6];
    double *x;

    //TODO: index by direction and behind ahead nb = 0,1 like advection
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double aII;
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;

    //TODO: Need everything for setIndexMap() to implement (global indexing for petsc).
    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
    //TODO: fix above and pick up here.

	for (l = 0; l < dim; ++l)
	{
            PETSc solver;
            solver.Create(ilower, iupper-1, 7, 7);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ijk_to_I[i][j][k];
            	if (I == -1) continue;

            	index  = d_index3d(i,j,k,top_gmax);
            	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

            	I_nb[0] = ijk_to_I[i-1][j][k]; //west
            	I_nb[1] = ijk_to_I[i+1][j][k]; //east
            	I_nb[2] = ijk_to_I[i][j-1][k]; //south
            	I_nb[3] = ijk_to_I[i][j+1][k]; //north
            	I_nb[4] = ijk_to_I[i][j][k-1]; //lower
            	I_nb[5] = ijk_to_I[i][j][k+1]; //upper


            	mu0   = field->mu[index];
            	rho   = field->rho[index];

        //TODO: index by dir = 0,1,2 and nb = 0,1
        for (nb = 0; nb < 6; nb++)
        {
		    if ((*findStateAtCrossing)(front,icoords,dir[nb],comp,
                                &intfc_state,&hs,crx_coords))
		    {
                if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                                boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                {
                    U_nb[nb] = vel[l][index];
                }
                else
                    U_nb[nb] = getStateVel[l](intfc_state);

                if (wave_type(hs) == DIRICHLET_BOUNDARY || neumann_type_bdry(wave_type(hs)))
                    mu[nb] = mu0;
                else
                    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    
            }
            else
		    {
                U_nb[nb] = vel[l][index_nb[nb]];
    			mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
        }

                //       should be discretizing div(grad(u) + grad(u)^T)
            	coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

                //TODO: RHS should also contain the advective flux...
            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		aII = 1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5];
		rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*(vel[l][index]);

        //TODO: index by dir = 0,1,2 and nb = 0,1
		for(nb = 0; nb < 6; nb++)
		{
            //TODO: Neumann boundary at solid walls does not appear to be implemented.
            if (!(*findStateAtCrossing)(front,icoords,dir[nb],comp,
			        &intfc_state,&hs,crx_coords))
		    {
                solver.Set_A(I,I_nb[nb],-coeff[nb]);
                rhs += coeff[nb]*U_nb[nb];
		    }
		    else
		    {
                if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                                boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                {
                    aII -= coeff[nb];
                    rhs += coeff[nb]*U_nb[nb];
                }
                else
                    rhs += 2.0*coeff[nb]*U_nb[nb];
		    }
		}
	
        rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
		    //rhs -= m_dt*grad_q[l][index]/rho;
        solver.Set_A(I,I,aII);
		solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-14);

	    start_clock("Befor Petsc solve");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"computeDiffusionCN: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
        }
	FT_ParallelExchGridVectorArrayBuffer(vel,front);

        FT_FreeThese(1,x);
}
*/

void G_CARTESIAN::computeAdvection()
{
	int order;
	switch (eqn_params->num_scheme)
	{
	case TVD_FIRST_ORDER:
	case WENO_FIRST_ORDER:
	    nrad = 3;
	    order = 1;
	    break;
	case TVD_SECOND_ORDER:
	case WENO_SECOND_ORDER:
	    nrad = 3;
	    order = 2;
	    break;
	case TVD_FOURTH_ORDER:
	case WENO_FOURTH_ORDER:
	    nrad = 3;
	    order = 4;
	    break;
	default:
	    order = -1;
	}
	solveRungeKutta(order);
}	/* end computeAdvection */


void G_CARTESIAN::solveRungeKutta(int order)
{
	static SWEEP *st_field,st_tmp;
	static FSWEEP *st_flux;
	static double **a,*b;
	double delta_t;
	int i,j;

	/* Allocate memory for Runge-Kutta of order */
	start_clock("solveRungeKutta");
	if (st_flux == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&b,order,sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&a,order,order,sizeof(double));

	    FT_VectorMemoryAlloc((POINTER*)&st_field,order,sizeof(SWEEP));
	    FT_VectorMemoryAlloc((POINTER*)&st_flux,order,sizeof(FSWEEP));
	    for (i = 0; i < order; ++i)
	    {
	    	allocMeshVst(&st_tmp);
	    	allocMeshVst(&st_field[i]);
	    	allocMeshFlux(&st_flux[i]);
	    }
	    /* Set coefficient a, b, c for different order of RK method */
	    switch (order)
	    {
	    case 1:
            b[0] = 1.0;
	    	break;
	    case 2:
	    	a[0][0] = 1.0;
	    	b[0] = 0.5;  b[1] = 0.5;
	    	break;
	    case 4:
	    	a[0][0] = 0.5;
	    	a[1][0] = 0.0;  a[1][1] = 0.5;
	    	a[2][0] = 0.0;  a[2][1] = 0.0;  a[2][2] = 1.0;
	    	b[0] = 1.0/6.0;  b[1] = 1.0/3.0;
	    	b[2] = 1.0/3.0;  b[3] = 1.0/6.0;
	    	break;
	    default:
	    	(void)printf("ERROR: %d-th order RK method not implemented\n",
					order);
	    	clean_up(ERROR);
	    }
	}
	delta_t = m_dt;

	/* Compute flux and advance field */

	copyToMeshVst(&st_field[0]);
	computeMeshFlux(st_field[0],&st_flux[0],delta_t);
    
    //TODO: Add function computeMeshViscFlux(),
    //      and add to the total mesh flux appropriately.
	
	for (i = 0; i < order-1; ++i)
	{
	    copyMeshVst(st_field[0],&st_field[i+1]);
	    for (j = 0; j <= i; ++j)
	    {
            if (a[i][j] != 0.0)
                addMeshFluxToVst(&st_field[i+1],st_flux[j],a[i][j]);
	    }
	    computeMeshFlux(st_field[i+1],&st_flux[i+1],delta_t);
	}

	for (i = 0; i < order; ++i)
	{
	    if (b[i] != 0.0)
            addMeshFluxToVst(&st_field[0],st_flux[i],b[i]);
	}
	copyFromMeshVst(st_field[0]);
	stop_clock("solveRungeKutta");
}	/* end solveRungeKutta */

void G_CARTESIAN::computeMeshFlux(
	SWEEP m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	if(eqn_params->tracked)
	{
	    start_clock("get_ghost_state");
	    get_ghost_state(m_vst, 2, 0);
	    get_ghost_state(m_vst, 3, 1);
	    scatMeshGhost();
	    stop_clock("get_ghost_state");
	    start_clock("solve_exp_value");
	    solve_exp_value();//sets gnor array (interface normal vectors)
	    stop_clock("solve_exp_value");
	}

    resetFlux(m_flux);
	for (int dir = 0; dir < dim; ++dir)
	{
	    addFluxInDirection(dir,&m_vst,m_flux,delta_t);
	}
	addSourceTerm(m_vst,m_flux,delta_t);
}	/* end computeMeshFlux */

void G_CARTESIAN::resetFlux(FSWEEP *m_flux)
{
	int i,j;
	int size = (int)cell_center.size();
	for (i = 0; i < size; i++)
	{
	    m_flux->dens_flux[i] = 0.0;
	    m_flux->engy_flux[i] = 0.0;
	    for (j = 0; j < MAXD; ++j)
	    	m_flux->momn_flux[j][i] = 0.0;
	}
}	/* resetFlux */

void G_CARTESIAN::addFluxInDirection(
	int dir,
	SWEEP *m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,j,icoords[MAXD];
        int num_thread = get_num_of_thread();
	switch (dim)
	{
	case 1:
	    addFluxAlongGridLine(dir,icoords,delta_t,m_vst,m_flux);
	    break;
	case 2:
	    for (i = imin[(dir+1)%dim]; i <= imax[(dir+1)%dim]; ++i)
	    {
	    	icoords[(dir+1)%dim] = i;
	    	addFluxAlongGridLine(dir,icoords,delta_t,m_vst,m_flux);
	    }
	    break;
	case 3:
	    for (i = imin[(dir+1)%dim]; i <= imax[(dir+1)%dim]; ++i)
	    for (j = imin[(dir+2)%dim]; j <= imax[(dir+2)%dim]; ++j)
	    {
	    	icoords[(dir+1)%dim] = i;
	    	icoords[(dir+2)%dim] = j;
	    	addFluxAlongGridLine(dir,icoords,delta_t,m_vst,m_flux);
	    }
	    break;
	}
}	/* end addFluxInDirection */

void G_CARTESIAN::scatMeshFlux(FSWEEP *m_flux)
{
	int i,j,k,l,index;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index1d(i,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index1d(i,op_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index2d(i,j,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index2d(i,j,top_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = m_flux->dens_flux[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
		m_flux->dens_flux[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		array[index] = m_flux->engy_flux[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_flux->engy_flux[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (k = imin[2]; k <= imax[2]; ++k)
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index3d(i,j,k,top_gmax);
		    array[index] = m_flux->momn_flux[l][index];
	    	}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (k = 0; k <= top_gmax[2]; k++)
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
		    index = d_index3d(i,j,k,top_gmax);
		    m_flux->momn_flux[l][index] = array[index];
	    	}
	    }
	}
}	/* end scatMeshFlux */

void G_CARTESIAN::addSourceTerm(
	const SWEEP& m_vst,
	FSWEEP *m_flux,
	double delta_t)
{
	int i,j,k,l,index;
	double *gravity = eqn_params->gravity;

	switch (dim)
	{
	case 1:
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index1d(i,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst.dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst.momn[l][index];
		    }
		}
	    }
	    break;
	case 2:
            for (j = imin[1]; j <= imax[1]; j++)
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index2d(i,j,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst.dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst.momn[l][index];
		    }
		}
	    }
	    break;
	case 3:
            for (k = imin[2]; k <= imax[2]; k++)
            for (j = imin[1]; j <= imax[1]; j++)
            for (i = imin[0]; i <= imax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
		if (!gas_comp(top_comp[index]))
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] = 0.0; 
		    	m_flux->engy_flux[index] = 0.0; 
		    }
		}
		else
		{
		    for (l = 0; l < dim; ++l)
		    {
		    	m_flux->momn_flux[l][index] += 
				delta_t*gravity[l]*m_vst.dens[index];
		    	m_flux->engy_flux[index] += 
				delta_t*gravity[l]*m_vst.momn[l][index];
		    }
		}
	    }
	}
	
}	/* end addSourceTerm */

void G_CARTESIAN::solve(double dt)
{
	m_dt = dt;
	max_speed = 0.0;

	if (debugging("trace"))
	    printf("Entering solve()\n");
	start_clock("solve");
	
    setDomain();
	appendOpenEndStates(); /* open boundary test */
	scatMeshStates();

	adjustGFMStates();
	setComponent();

    ////////////////////////////////////////////////////////////////
    //TODO: part of diffusion implementation
        //setGlobalIndex();
    ////////////////////////////////////////////////////////////////
	
	if (debugging("trace"))
	    printf("Passed setComponent()\n");

    //TODO: Need to save the current solution, and keep in a seperate array
    //      the solution of the predictor step for use in step 2.
	
    // 1) Explicit Predictor Step
	start_clock("computeAdvection");
    computeAdvection();
	    //computeConvectiveFlux();
	if (debugging("trace"))
	    printf("max_speed after computeAdvection(): %20.14f\n",max_speed);
	stop_clock("computeAdvection");

    /////////////////////////////////////////////////////////////// 
    // 2) Implicit Corrector Step
    //start_clock("computeDiffusion");
        //computeDiffusion(); //TODO: diffusion implementation here
    //stop_clock("computeDiffusion");

    ///////////////////////////////////////////////////////////////

    /*if (debugging("sample_velocity"))
	{
	    sampleVelocity();
	}*/

    //TODO: Velocity and Vorticity computed in copyMeshStates(),
    //      should factor into separate functions.
	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	setAdvectionDt();
	stop_clock("solve");
	if (debugging("trace"))
	    printf("Leaving solve()\n");
}	/* end solve */


// check http://en.wikipedia.org/wiki/Bilinear_interpolation
void G_CARTESIAN::getVelocity(double *p, double *U)
{
    double **vel = eqn_params->vel;

    FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[0],getStateXvel,&U[0],
                NULL);
    if (dim > 1)
        FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[1],getStateYvel,&U[1],
                NULL);
    if (dim > 2)
        FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[2],getStateZvel,&U[2],
                NULL);
}

void G_CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void G_CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int G_CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void G_CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].m_coords[i];
}

void G_CARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(cell_center[index0].m_coords[i] +
	    		     cell_center[index1].m_coords[i]);
	}
}


double G_CARTESIAN::getDistance(double *c0, double *c1)
{
	return sqrt( (c0[0]-c1[0])*(c0[0]-c1[0])
		    +(c0[1]-c1[1])*(c0[1]-c1[1]) );
}


// input : p[]
// output: q[]

void G_CARTESIAN::getNearestInterfacePoint(
	double *p, 
	double *q)
{
	INTERFACE *intfc = front->interf;
	double t;
	HYPER_SURF_ELEMENT *phse;
	HYPER_SURF *phs;
	nearest_interface_point(p,getComponent(p),intfc,NO_BOUNDARIES,
				NULL,q,&t,&phse,&phs);
}

int G_CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int G_CARTESIAN::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	default:
	    return NO_COMP;
	}
}

void G_CARTESIAN::save(char *filename)
{
	
	INTERFACE *intfc    = front->interf;
		
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

void G_CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i,j,k;
	static int size;

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	eqn_params = (EQN_PARAMS*)front->extra1;
	
	if (first)
	{
	    first = NO;
	    dim = grid_intfc->dim;

	    hmin = HUGE;
	    size = 1;
	    
            for (i = 0; i < 3; ++i)
	    	top_gmax[i] = 0;

            for (i = 0; i < dim; ++i)
	    {
	    	lbuf[i] = front->rect_grid->lbuf[i];
	    	ubuf[i] = front->rect_grid->ubuf[i];
	    	top_gmax[i] = top_grid->gmax[i];
	    	top_L[i] = top_grid->L[i];
	    	top_U[i] = top_grid->U[i];
	    	top_h[i] = top_grid->h[i];

                if (hmin > top_h[i]) hmin = top_h[i];
	        size *= (top_gmax[i]+1);
	    	imin[i] = (lbuf[i] == 0) ? 1 : lbuf[i];
	    	imax[i] = (ubuf[i] == 0) ? top_gmax[i] - 1 : 
				top_gmax[i] - ubuf[i];
	    }

	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->dens,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->pres,size,
					sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&eqn_params->engy,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->vel,dim,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->mom,dim,size,
					sizeof(double));
	    //GFM
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->gnor,dim,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->Gdens,2,size,
					sizeof(double));
	    FT_MatrixMemoryAlloc((POINTER*)&eqn_params->Gpres,2,size,
					sizeof(double));
	    FT_TriArrayMemoryAlloc((POINTER*)&eqn_params->Gvel,2,dim,size,
					sizeof(double));

	    FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
	    
        if (dim == 2)
	    	FT_VectorMemoryAlloc((POINTER*)&eqn_params->vort,size,
                    sizeof(double));
        else if (dim == 3)
            FT_MatrixMemoryAlloc((POINTER*)&eqn_params->vorticity,
                    dim,size,sizeof(double));

	    field.dens = eqn_params->dens;
	    field.engy = eqn_params->engy;
	    field.pres = eqn_params->pres;
	    field.momn = eqn_params->mom;
	    field.vel = eqn_params->vel;
	    
        if (dim == 3)
            field.vorticity = eqn_params->vorticity;
	}

    //GFM
	for (i = 0; i < size; ++i)
	for (j = 0; j < 2; ++j)
	{
	    eqn_params->Gdens[j][i] = 0.0;
	    eqn_params->Gpres[j][i] = 0.0;
	    for (k = 0; k < dim; ++k)
		eqn_params->Gvel[j][k][i] = 0.0;
	}
}

void G_CARTESIAN::allocMeshVst(
	SWEEP *vst)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    size *= (top_gmax[i]+1);

	FT_VectorMemoryAlloc((POINTER*)&vst->dens,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->engy,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->pres,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&vst->momn,MAXD,size,sizeof(double));
}	/* end allocMeshVstFlux */

void G_CARTESIAN::allocMeshFlux(
	FSWEEP *flux)
{
	int i,size;

	size = 1;
        for (i = 0; i < dim; ++i)
	    size *= (top_gmax[i]+1);

	FT_VectorMemoryAlloc((POINTER*)&flux->dens_flux,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux->engy_flux,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&flux->momn_flux,MAXD,size,sizeof(double));
}	/* end allocMeshVstFlux */

void G_CARTESIAN::allocDirVstFlux(
        SWEEP *vst,
        FSWEEP *flux)
{
	int size = 0;
    for (int i = 0; i < dim; ++i)
    {
	    if (size < top_gmax[i]+7) 
            size = top_gmax[i]+7;
    }

	FT_VectorMemoryAlloc((POINTER*)&vst->dens,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->engy,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&vst->pres,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&vst->momn,MAXD,size,sizeof(double));

	FT_VectorMemoryAlloc((POINTER*)&flux->dens_flux,size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux->engy_flux,size,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&flux->momn_flux,MAXD,size,sizeof(double));
}	/* end allocDirMeshVstFlux */

void G_CARTESIAN::freeDirVstFlux(
        SWEEP* vst,
        FSWEEP* flux)
{
        FT_FreeThese(4,vst->dens,vst->engy,vst->pres,vst->momn);
        FT_FreeThese(3,flux->dens_flux,flux->engy_flux,flux->momn_flux);
}	/* end allocDirMeshVstFlux */

void G_CARTESIAN::checkVst(SWEEP *vst)
{
	int i,j,index;
	for (j = imin[1]; j < imax[1]; j++)
	for (i = imin[0]; i < imax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    if (isnan(vst->dens[index]))
		printf("At %d %d: dens is nan\n",i,j);
	    if (vst->dens[index] < 0.0)
		printf("At %d %d: dens is negative\n",i,j);
	}
}

void G_CARTESIAN::checkFlux(FSWEEP *flux)
{
	int i,j,index;
	//j = 140;
	for (j = imin[1]; j < imax[1]; j++)
	for (i = imin[0]; i <= imax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    printf("%d %f  %f\n",i,flux->momn_flux[1][index],
				flux->engy_flux[index]);
	}
}

void G_CARTESIAN::printGasStates(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double *dens = field.dens;
	double *engy = field.engy;
	double **momn = field.momn;

	sprintf(filename,"%s/state.ts%s",out_name,
			right_flush(front->step,7));
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
	sprintf(filename,"%s-gas",filename);
	outfile = fopen(filename,"w");

    int WID = 24; int DEC = 18;
    if (debugging("integration_test")) {DEC = 1;}

    /* Initialize states at the interface */
    fprintf(outfile,"Interface gas states:\n");
    next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
        FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
        fprintf(outfile,"%*.*f %*.*f\n",WID,DEC,getStateDens(sl),
            WID,DEC,getStateDens(sr));
        fprintf(outfile,"%*.*f %*.*f\n",WID,DEC,getStateEngy(sl),
            WID,DEC,getStateEngy(sr));
    for (i = 0; i < dim; ++i)
            fprintf(outfile,"%*.*f %*.*f\n",WID,DEC,getStateMom[i](sl),
            WID,DEC,getStateMom[i](sr));
    }
	
	fprintf(outfile,"\nInterior gas states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%*.*f\n",WID,DEC,dens[index]);
	        fprintf(outfile,"%*.*f\n",WID,DEC,engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%*.*f\n",WID,DEC,momn[l][index]);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%*.*f\n",WID,DEC,dens[index]);
	        fprintf(outfile,"%*.*f\n",WID,DEC,engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%*.*f\n",WID,DEC,momn[l][index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%*.*f\n",WID,DEC,dens[index]);
	        fprintf(outfile,"%*.*f\n",WID,DEC,engy[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%*.*f\n",WID,DEC,momn[l][index]);
	    }
	}
	fclose(outfile);
}

void G_CARTESIAN::readGasStates(char* restart_state_name)
{
    readFrontStates(front,restart_state_name);
    readInteriorStates(restart_state_name);
}

void G_CARTESIAN::readFrontStates(
        Front* frt,
        char* restart_state_name)
{
	FILE 		*infile;
	EQN_PARAMS 	*eqn_params = (EQN_PARAMS*)frt->extra1;
	INTERFACE 	*intfc = frt->interf;
        STATE 		*sl,*sr;
        POINT 		*p;
        HYPER_SURF 	*hs;
        HYPER_SURF_ELEMENT *hse;
	STATE 		*lstate,*rstate;
	char 		fname[100];
	int 		i,dim = frt->rect_grid->dim;
	int		comp;
	EOS_PARAMS	*eos = eqn_params->eos;

	sprintf(fname,"%s-gas",restart_state_name);
	infile = fopen(fname,"r");
	
	/* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface gas states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    lstate = (STATE*)sl;	rstate = (STATE*)sr;
            fscanf(infile,"%lf %lf",&lstate->dens,&rstate->dens);
            fscanf(infile,"%lf %lf",&lstate->engy,&rstate->engy);
	    for (i = 0; i < dim; ++i)
            	fscanf(infile,"%lf %lf",&lstate->momn[i],&rstate->momn[i]);
	    
	    comp = negative_component(hs);
	    lstate->eos = &eos[comp];
	    lstate->dim = dim;
	    if(gas_comp(comp))
	    	lstate->pres = EosPressure(lstate);
		
	    comp = positive_component(hs);

	    rstate->eos = &eos[comp];
	    rstate->dim = dim;
	    if(gas_comp(comp))
	    	rstate->pres = EosPressure(rstate);
	    lstate->dim = rstate->dim = dim;
        }
	FT_MakeGridIntfc(frt);
	fclose(infile);
}

void G_CARTESIAN::readInteriorStates(char* restart_state_name)
{
	FILE *infile;
	int i,j,k,l,index;
	STATE st_tmp;
	char fname[100];
	int		comp;
	EOS_PARAMS	*eos = eqn_params->eos;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;

	setDomain();
	m_dens[0] = eqn_params->rho1;		
	m_dens[1] = eqn_params->rho2;		
	m_mu[0] = eqn_params->mu1;		
	m_mu[1] = eqn_params->mu2;		
	if (eqn_params->prob_type == FLUID_SOLID_CIRCLE ||
	    eqn_params->prob_type == FLUID_RIGID_BODY ||
	    eqn_params->prob_type == FLUID_CRYSTAL)
	    m_comp[0] = SOLID_COMP;
	else
	    m_comp[0] = GAS_COMP1;
	m_comp[1] = GAS_COMP2;
	m_smoothing_radius = top_h[0] < top_h[1] ? top_h[1] : top_h[0];
	m_smoothing_radius *= 2.0;
	
	st_tmp.dim = eqn_params->dim;

	sprintf(fname,"%s-gas",restart_state_name);
	infile = fopen(fname,"r");
	
	next_output_line_containing_string(infile,"Interior gas states:");

	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		st_tmp.eos = &(eos[comp]);
	    	
		fscanf(infile,"%lf",&dens[index]);
	    	fscanf(infile,"%lf",&engy[index]);
		st_tmp.dens = dens[index];
		st_tmp.engy = engy[index];
		for (l = 0; l < dim; ++l)
		{
	    	    fscanf(infile,"%lf",&momn[l][index]);
		    st_tmp.momn[l] = momn[l][index];
		}
		pres[index] = EosPressure(&st_tmp);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		st_tmp.eos = &(eos[comp]);

	    	fscanf(infile,"%lf",&dens[index]);
	    	fscanf(infile,"%lf",&engy[index]);
		st_tmp.dens = dens[index];
		st_tmp.engy = engy[index];
		for (l = 0; l < dim; ++l)
		{
	    	    fscanf(infile,"%lf",&momn[l][index]);
		    st_tmp.momn[l] = momn[l][index];
		}
		pres[index] = EosPressure(&st_tmp);
	    }
	}
	fclose(infile);
	scatMeshStates();
	copyMeshStates();
}


void G_CARTESIAN::setAdvectionDt()
{
	double d = (double)dim;
	pp_global_max(&max_speed,1);
	if (max_speed != 0.0)
	    max_dt = hmin/max_speed/d;
	else
	    max_dt = 0.0;
	if (debugging("trace"))
	    printf("In setAdvectionDt: max_dt = %24.18g\n",max_dt);
}	/* end setAdvectionDt */


void G_CARTESIAN::augmentMovieVariables()
{
	int i;
	static HDF_MOVIE_VAR *hdf_movie_var;
	int offset,num_var;

	hdf_movie_var = front->hdf_movie_var;
	offset = front->hdf_movie_var->num_var;
	if (hdf_movie_var == NULL)
	    return initMovieVariables();
	else
	{
	    num_var = offset + dim + 3;
	    FT_ScalarMemoryAlloc((POINTER*)&hdf_movie_var,
				sizeof(HDF_MOVIE_VAR));
	    hdf_movie_var->num_var = num_var;
	    FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var->var_name,
				num_var,100,sizeof(char));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->top_var,
				num_var,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->obstacle_comp,
				num_var,sizeof(COMPONENT));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->preset_bound,
				num_var,sizeof(boolean));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_min,
				num_var,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var->var_max,
				num_var,sizeof(double));
	    for (i = 0; i < front->hdf_movie_var->num_var; ++i)
	    {
		strcpy(hdf_movie_var->var_name[i],
			front->hdf_movie_var->var_name[i]);
		hdf_movie_var->get_state_var[i] =
			front->hdf_movie_var->get_state_var[i];
		hdf_movie_var->top_var[i] = 
			front->hdf_movie_var->top_var[i];
	    }
	    sprintf(hdf_movie_var->var_name[offset+0],"dens");
	    sprintf(hdf_movie_var->var_name[offset+1],"pres");
	    sprintf(hdf_movie_var->var_name[offset+2],"vort");
	    sprintf(hdf_movie_var->var_name[offset+3],"xvel");
	    sprintf(hdf_movie_var->var_name[offset+4],"yvel");
	    hdf_movie_var->get_state_var[offset+0] = getStateDens;
	    hdf_movie_var->get_state_var[offset+1] = getStatePres;
	    hdf_movie_var->get_state_var[offset+2] = getStateVort;
	    hdf_movie_var->get_state_var[offset+3] = getStateXvel;
	    hdf_movie_var->get_state_var[offset+4] = getStateYvel;
	    if (dim == 3)
	    {
	    	sprintf(hdf_movie_var->var_name[offset+5],"zvel");
	    	hdf_movie_var->get_state_var[offset+5] = getStateZvel;
	    }
	}
	hdf_movie_var->top_var[offset+0] = eqn_params->dens;
	hdf_movie_var->top_var[offset+1] = eqn_params->pres;
	hdf_movie_var->top_var[offset+2] = eqn_params->vort;
	hdf_movie_var->top_var[offset+3] = eqn_params->vel[0];
	hdf_movie_var->top_var[offset+4] = eqn_params->vel[1];
	if (dim == 3)
	    hdf_movie_var->top_var[offset+5] = eqn_params->vel[2];
	FT_FreeThese(2,front->hdf_movie_var->var_name,
			front->hdf_movie_var->top_var);
	FT_FreeThese(1,front->hdf_movie_var);
	front->hdf_movie_var = hdf_movie_var;
	front->hdf_movie_var->num_var = num_var;

}	/* end augmentMovieVariables */

void G_CARTESIAN::initMovieVariables()
{
	boolean set_bound = NO;
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	double var_max,var_min;

	if (CursorAfterStringOpt(infile,"Type y to set movie bounds:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                set_bound = YES;
        }

	    /* Begin hdf movies */
	switch (dim)
	{
	case 1:
	    CursorAfterStringOpt(infile,"Type y to make movie of density:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterStringOpt(infile,"Enter min and max density:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"dens",0,eqn_params->dens,getStateDens,
				var_max,var_min);
	    }
	    CursorAfterStringOpt(infile,"Type y to make movie of pressure:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max pressure:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres",0,eqn_params->pres,getStatePres,
				var_max,var_min);
	    }
	    CursorAfterStringOpt(infile,"Type y to make movie of velocity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max velocity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo",0,eqn_params->vel[0],getStateXvel,
				var_max,var_min);
	    }
	    break;
	case 2:
	    CursorAfterStringOpt(infile,"Type y to make movie of density:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max density:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"dens",0,eqn_params->dens,getStateDens,
				var_max,var_min);
	    }
	    CursorAfterStringOpt(infile,"Type y to make movie of pressure:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max pressure:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres",0,eqn_params->pres,getStatePres,
				var_max,var_min);
	    }
	    CursorAfterStringOpt(infile,"Type y to make movie of vorticity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max vorticity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"vort",0,eqn_params->vort,getStateVort,
				var_max,var_min);
	    }
	    CursorAfterStringOpt(infile,"Type y to make movie of velocity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max velocity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"xvel",0,eqn_params->vel[0],getStateXvel,
				var_max,var_min);
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"yvel",0,eqn_params->vel[1],getStateYvel,
				var_max,var_min);
	    }
	    break;
	case 3:
	    CursorAfterStringOpt(infile,"Type y to make yz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterStringOpt(infile,"Type y to make movie of density:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"dens-yz",0,eqn_params->dens,getStateDens,
				0.0,0.0);
	    	CursorAfterStringOpt(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-yz",0,eqn_params->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterStringOpt(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-yz-y",0,eqn_params->vel[1],getStateYvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-yz-z",0,eqn_params->vel[2],getStateZvel,
				0.0,0.0);
		}
	    }
	    CursorAfterStringOpt(infile,"Type y to make xz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterStringOpt(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-xz",1,eqn_params->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterStringOpt(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xz-x",1,eqn_params->vel[0],getStateXvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xz-z",1,eqn_params->vel[2],getStateZvel,
				0.0,0.0);
		}
	    }
	    CursorAfterStringOpt(infile,"Type y to make xy cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterStringOpt(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-xy",2,eqn_params->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterStringOpt(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xy-x",2,eqn_params->vel[0],getStateXvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xy-y",2,eqn_params->vel[1],getStateYvel,
				0.0,0.0);
		}
	    }
	}

	if (dim != 1)
	{
	    if (CursorAfterStringOpt(infile,
               "Type y to make scalar density field movie:"))
            {
                fscanf(infile,"%s",string);
                (void)printf("%s\n",string);
                if (string[0] == 'Y' || string[0] == 'y')
                    FT_AddVtkScalarMovieVariable(front,"DENSITY",field.dens);
            }
	    if (CursorAfterStringOpt(infile,
               "Type y to make scalar pressure field movie:"))
            {
                fscanf(infile,"%s",string);
                (void)printf("%s\n",string);
                if (string[0] == 'Y' || string[0] == 'y')
                    FT_AddVtkScalarMovieVariable(front,"PRESSURE",field.pres);
            }
	    if (CursorAfterStringOpt(infile,
		"Type y to make vector velocity field movie:"))
	    {
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
                    FT_AddVtkVectorMovieVariable(front,"VELOCITY",field.vel);
	    }

        if (dim == 3)
        {
            if (CursorAfterStringOpt(infile,
            "Type y to make vector vorticity field movie:"))
            {
                fscanf(infile,"%s",string);
                (void) printf("%s\n",string);
                if (string[0] == 'Y' || string[0] == 'y')
                    FT_AddVtkVectorMovieVariable(front,"VORTICITY",field.vorticity);
            }
        }

	}

	fclose(infile);
}	/* end initMovieVariables */


void G_CARTESIAN::computeVorticity()
{
    int index;
    int icoords[MAXD];

    int index_xnb0, index_xnb1;
    int icnb_x0[MAXD], icnb_x1[MAXD];

    int index_ynb0, index_ynb1;
    int icnb_y0[MAXD], icnb_y1[MAXD];

    int index_znb0, index_znb1;
    int icnb_z0[MAXD], icnb_z1[MAXD];
	
    double **vel = field.vel;
	double **vorticity = field.vorticity;

	for (int k = imin[2]; k <= imax[2]; k++)
	for (int j = imin[1]; j <= imax[1]; j++)
    for (int i = imin[0]; i <= imax[0]; i++)
	{
        index = d_index3d(i,j,k,top_gmax);
        if (!gas_comp(top_comp[index]))//TODO: or use cell_center[index].comp ??;
        {
            for (int l = 0; l < 3; ++l)
                vorticity[l][index] = 0.0;
            continue;
        }

        icoords[0] = i;
        icoords[1] = j;
        icoords[2] = k;

        for (int l = 0; l < dim; ++l)
        {
            icnb_x0[l] = icoords[l];
            icnb_x1[l] = icoords[l];
            icnb_y0[l] = icoords[l];
            icnb_y1[l] = icoords[l];
            icnb_z0[l] = icoords[l];
            icnb_z1[l] = icoords[l];
        }

        //cell centered derivative indices wrt x
        icnb_x0[0] = icoords[0] - 1;
        icnb_x1[0] = icoords[0] + 1;
        index_xnb0 = d_index(icnb_x0,top_gmax,dim);
        index_xnb1 = d_index(icnb_x1,top_gmax,dim);

        //cell centered derivative indices wrt y
        icnb_y0[1] = icoords[1] - 1;
        icnb_y1[1] = icoords[1] + 1;
        index_ynb0 = d_index(icnb_y0,top_gmax,dim);
        index_ynb1 = d_index(icnb_y1,top_gmax,dim);

        //cell centered derivative indices wrt z
        icnb_z0[2] = icoords[2] - 1;
        icnb_z1[2] = icoords[2] + 1;
        index_znb0 = d_index(icnb_z0,top_gmax,dim);
        index_znb1 = d_index(icnb_z1,top_gmax,dim);


        //x component vorticity
        double u2_wrty = 0.5*(vel[2][index_ynb1] - vel[2][index_ynb0])/top_h[1];
        double u1_wrtz = 0.5*(vel[1][index_znb1] - vel[1][index_znb0])/top_h[2];
        vorticity[0][index] = u2_wrty - u1_wrtz;

        //y component vorticity
        double u0_wrtz = 0.5*(vel[0][index_znb1] - vel[0][index_znb0])/top_h[2];
        double u2_wrtx = 0.5*(vel[2][index_xnb1] - vel[2][index_xnb0])/top_h[0];
        vorticity[1][index] = u0_wrtz - u2_wrtx;

        //z component vorticity
        double u1_wrtx = 0.5*(vel[1][index_xnb1] - vel[1][index_xnb0])/top_h[0];
        double u0_wrty = 0.5*(vel[0][index_ynb1] - vel[0][index_ynb0])/top_h[1];
        vorticity[2][index] = u1_wrtx - u0_wrty;
    }
	
    //FT_ParallelExchGridVectorArrayBuffer(vorticity,front);
}

double G_CARTESIAN::getVorticityX(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dy,dz;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dy = top_h[1];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j-1,k,top_gmax);
	index01 = d_index3d(i,j+1,k,top_gmax);
	index10 = d_index3d(i,j,k-1,top_gmax);
	index11 = d_index3d(i,j,k+1,top_gmax);
	v00 = -momn[2][index00]/dens[index00];
	v01 =  momn[2][index01]/dens[index01];
	v10 =  momn[1][index10]/dens[index10];
	v11 = -momn[1][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dz + (v10 + v11)/2.0/dy;
	return vorticity;
}	/* end getVorticityX */

double G_CARTESIAN::getVorticityY(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dz;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dz = top_h[2];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i,j,k-1,top_gmax);
	index01 = d_index3d(i,j,k+1,top_gmax);
	index10 = d_index3d(i-1,j,k,top_gmax);
	index11 = d_index3d(i+1,j,k,top_gmax);
	v00 = -momn[0][index00]/dens[index00];
	v01 =  momn[0][index01]/dens[index01];
	v10 =  momn[2][index10]/dens[index10];
	v11 = -momn[2][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dx + (v10 + v11)/2.0/dz;
	return vorticity;
}	/* end getVorticityY */

double G_CARTESIAN::getVorticityZ(int i, int j, int k)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index3d(i,j,k,top_gmax);
	index00 = d_index3d(i-1,j,k,top_gmax);
	index01 = d_index3d(i+1,j,k,top_gmax);
	index10 = d_index3d(i,j-1,k,top_gmax);
	index11 = d_index3d(i,j+1,k,top_gmax);
	v00 = -momn[1][index00]/dens[index00];
	v01 =  momn[1][index01]/dens[index01];
	v10 =  momn[0][index10]/dens[index10];
	v11 = -momn[0][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticityZ */

double G_CARTESIAN::getVorticity(int i, int j)
{
	int index0,index00,index01,index10,index11;
	double v00,v01,v10,v11;
	double dx,dy;
	double vorticity;
	double *dens = field.dens;
	double **momn = field.momn;

	dx = top_h[0];
	dy = top_h[1];
	index0 = d_index2d(i,j,top_gmax);
	index00 = d_index2d(i-1,j,top_gmax);
	index01 = d_index2d(i+1,j,top_gmax);
	index10 = d_index2d(i,j-1,top_gmax);
	index11 = d_index2d(i,j+1,top_gmax);
	v00 = -momn[1][index00]/dens[index00];
	v01 =  momn[1][index01]/dens[index01];
	v10 =  momn[0][index10]/dens[index10];
	v11 = -momn[0][index11]/dens[index11];

	vorticity = (v00 + v01)/2.0/dy + (v10 + v11)/2.0/dx;
	return vorticity;
}	/* end getVorticity */

void G_CARTESIAN::copyMeshStates()
{
	int i,j,k,l,index;
	double **vel = eqn_params->vel;
	double **mom = eqn_params->mom;
	double *dens = eqn_params->dens;
	double *pres = eqn_params->pres;
	double *engy = eqn_params->engy;
	double *vort = eqn_params->vort;
	double **vorticity = eqn_params->vorticity;
	int symmetry[MAXD];

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front,NULL);
	    FT_ParallelExchGridArrayBuffer(pres,front,NULL);
	    FT_ParallelExchGridArrayBuffer(engy,front,NULL);
	    symmetry[0] = ODD;
	    FT_ParallelExchGridArrayBuffer(mom[0],front,symmetry);
	    FT_ParallelExchGridArrayBuffer(vel[0],front,symmetry);
	    break;
	case 2:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
		for (l = 0; l < dim; ++l)
		    vel[l][index] = mom[l][index]/dens[index];	
	    }
	    FT_ParallelExchGridArrayBuffer(dens,front,NULL);
	    FT_ParallelExchGridArrayBuffer(pres,front,NULL);
	    FT_ParallelExchGridArrayBuffer(engy,front,NULL);
	    for (l = 0; l < dim; ++l)
	    {
	    	symmetry[0] = symmetry[1] = EVEN;
		    symmetry[l] = ODD;
	    	FT_ParallelExchGridArrayBuffer(mom[l],front,symmetry);
	    	FT_ParallelExchGridArrayBuffer(vel[l],front,symmetry);
	    }
	    symmetry[0] = symmetry[1] = ODD;
	    FT_ParallelExchGridArrayBuffer(vort,front,symmetry);
	    break;
	case 3:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (k = imin[2]; k <= imax[2]; ++k)
	    {
            index = d_index3d(i,j,k,top_gmax);
            for (l = 0; l < dim; ++l)
                vel[l][index] = mom[l][index]/dens[index];	
	    }
        computeVorticity();

	    FT_ParallelExchGridArrayBuffer(dens,front,NULL);
	    FT_ParallelExchGridArrayBuffer(pres,front,NULL);
	    FT_ParallelExchGridArrayBuffer(engy,front,NULL);
        
        for (l = 0; l < dim; ++l)
	    {
	    	symmetry[0] = symmetry[1] = symmetry[2] = EVEN;
		    symmetry[l] = ODD;
	    	FT_ParallelExchGridArrayBuffer(mom[l],front,symmetry);
	    	FT_ParallelExchGridArrayBuffer(vel[l],front,symmetry);
	    }
        //TODO: check vorticity parallel communication is correct
	    symmetry[0] = symmetry[1] = symmetry[2] = ODD;
	    for (l = 0; l < dim; ++l)
            FT_ParallelExchGridArrayBuffer(vorticity[l],front,symmetry);
	    break;
	}
}	/* end copyMeshStates */

void G_CARTESIAN::compSGS(void)
{
        int i,j,k,index,index0,index1,index2,index3,index4,size;  
        double *u, *v;
        double ulx,urx,vlx,vrx;
        double uly,ury,vly,vry;
        double ux,uy,vx,vy;
        double s, *s11, *s12, *s22;
        double *ss11, *ss12, *ss22;
        double *tau00, *tau01, *tau10, *tau11;
        double *vel_u, *vel_v, *vel_uu, *vel_uv, *vel_vv;
        double sum_vel_u,sum_vel_v,sum_vel_uu,sum_vel_uv,sum_vel_vv;
        double sum_s11,sum_s12,sum_s22,sum_ss11,sum_ss12,sum_ss22,sum_s;
        double *ma11, *ma12, *la11, *la12, *la22;
        double *cs, *cs_ave, *deno, *nume, *co_coords_y;
        double coords[2];
        int    *r, num_r;
        int    ii,jj,iii,jjj;
        const int nn = pp_numnodes();
        num_r = (int)(((top_U[1]-top_L[1])/top_h[1])+1);
	double **momn = field.momn;

        size = (top_gmax[0]+1)*(top_gmax[1]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau00,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau01,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau10,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uu,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_vv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&co_coords_y,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&r,size,sizeof(int));
        FT_VectorMemoryAlloc((POINTER*)&cs,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cs_ave,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&deno,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&nume,num_r,sizeof(double));

        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            u[index] = momn[0][index];
            v[index] = momn[1][index];
            getRectangleCenter(index, coords);
            co_coords_y[index] = coords[1] + (top_h[1]/2.0);
        }

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]-2; i <= imax[0]+2; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            index1  = d_index2d(i-1,j,top_gmax);
            ulx = u[index1];
            uly = v[index1];
            index2  = d_index2d(i+1,j,top_gmax);
            urx = u[index2];
            ury = v[index2];
            index3  = d_index2d(i,j-1,top_gmax);
            vlx = u[index3];
            vly = v[index3];
            index4  = d_index2d(i,j+1,top_gmax);
            vrx = u[index4];
            vry = v[index4];

            ux = (urx - ulx) / (2.0*top_h[0]);
            uy = (ury - uly) / (2.0*top_h[1]);
            vx = (vrx - vlx) / (2.0*top_h[0]);
            vy = (vry - vly) / (2.0*top_h[1]);
            s11[index0] = ux;
            s12[index0] = (uy + vx)/2;
            s22[index0] = vy;
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                         + (2*(s12[index0]*s12[index0]))
                         + (s22[index0]*s22[index0])));
            ss11[index0] = s*s11[index0];
            ss12[index0] = s*s12[index0];
            ss22[index0] = s*s22[index0];
            vel_u[index0] = u[index0];
            vel_v[index0] = v[index0];  
            vel_uu[index0] = u[index0]*u[index0]; 
            vel_uv[index0] = u[index0]*v[index0];  
            vel_vv[index0] = v[index0]*v[index0];      
        }

        for (j = imin[1]; j <= (imax[1]/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin[0]-1; i <= (imax[0]/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                sum_vel_u = sum_vel_v = 0.0;
                sum_vel_uu = sum_vel_uv = sum_vel_vv = 0.0;
                sum_s11 = sum_s12 = sum_s22 = 0.0;
                sum_ss11 = sum_ss12 = sum_ss22 = sum_s = 0.0;
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0  = d_index2d(iii,jjj,top_gmax);
                    sum_vel_u += vel_u[index0];
                    sum_vel_v += vel_v[index0];
                    sum_vel_uu += vel_uu[index0];
                    sum_vel_uv += vel_uv[index0];
                    sum_vel_vv += vel_vv[index0];
                    sum_s11 += s11[index0];
                    sum_s12 += s12[index0];
                    sum_s22 += s22[index0];
                    sum_ss11 += ss11[index0];
                    sum_ss12 += ss12[index0];
                    sum_ss22 += ss22[index0];
                    sum_s += sqrt(2*( (s11[index0]*s11[index0]) 
                                  + (2*(s12[index0]*s12[index0])) 
                                  + (s22[index0]*s22[index0])));
                } 
                ma11[index] = (2.0*top_h[1]*top_h[1]*(sum_ss11/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s11/4.0));
                ma12[index] = (2.0*top_h[1]*top_h[1]*(sum_ss12/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s12/4.0));
                la11[index] = (sum_vel_uu/4.0)-((sum_vel_u/4.0)*
			(sum_vel_u/4.0));
                la12[index] = (sum_vel_uv/4.0)-((sum_vel_u/4.0)*
			(sum_vel_v/4.0));
                la12[index] = (sum_vel_vv/4.0)-((sum_vel_v/4.0)*
			(sum_vel_v/4.0));
                r[index] = (int)(co_coords_y[index]/(2*top_h[1]));
            }
        }

        for (k = 0; k < num_r; k++)
        {
            deno[k] = 0.0;
            nume[k] = 0.0;
        }

        for (k = 0; k < num_r; k++)
        for (j = imin[1]; j <= (imax[1]/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin[0]-1; i <= (imax[0]/2)+1; i++)
            {
                ii = (2*i)-4;
                index0 = d_index2d(ii,jj,top_gmax);
                if(k == r[index0])
                {
                    deno[k] += (ma11[index0]*ma11[index0]) + 
				(ma12[index0]*ma12[index0]);
                    nume[k] += (((la11[index0]/2.0)-(la22[index0]/2.0))*
				ma11[index0]) + (la12[index0]*ma12[index0]);
                }
            }
        }

        pp_gsync();
        
        if (nn > 1)
        {
           for (k = 0; k < num_r; k++)
           {
              pp_global_sum(&deno[k],1L);
              pp_global_sum(&nume[k],1L);
           }
        }

        for (k = 0; k < num_r; k++)
        {
            if(deno[k] < 10e-16)
                cs_ave[k] = 0.0;
            else
                cs_ave[k] = nume[k]/deno[k];
        }

        for (j = imin[1]; j <= (imax[1]/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin[0]-1; i <= (imax[0]/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0 = d_index2d(iii,jjj,top_gmax);
                    cs[index0] = cs_ave[r[index]];
                }
            }
        }

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]-1; i <= imax[0]+1; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                          + (2*(s12[index0]*s12[index0]))
                          + (s22[index0]*s22[index0])));
            tau00[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s11[index0]/2.0)-(s22[index0]/2.0));
            tau01[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau10[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau11[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s22[index0]/2.0)-(s11[index0]/2.0));
        }

        for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
        {
            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            momn[0][index0] += -m_dt*(
                              ((tau00[index2]-tau00[index1])/(2.0*top_h[0])) + 
                                ((tau01[index4]-tau01[index3])/(2.0*top_h[1])));
            momn[1][index0] += -m_dt*(
                              ((tau10[index2]-tau10[index1])/(2.0*top_h[0])) + 
                              ((tau11[index4]-tau11[index3])/(2.0*top_h[1])));
        }
        FT_FreeThese(2,u,v);
        FT_FreeThese(4,tau00,tau01,tau10,tau11);
        FT_FreeThese(6,s11,s12,s22,ss11,ss12,ss22);
        FT_FreeThese(5,vel_u,vel_v,vel_uu,vel_uv,vel_vv);
        FT_FreeThese(11,co_coords_y,ma11,ma12,la11,la12,la22,r,cs,cs_ave,
					deno,nume);
}       /* end compSGS */

void G_CARTESIAN::sampleVelocity()
{
	switch (dim)
	{
	case 2:
	    return sampleVelocity2d();
	case 3:
	    return sampleVelocity3d();
	}
}	/* end sampleVelocity */

void G_CARTESIAN::sampleVelocity3d()
{
        int i,j,k,index;
        double coords[MAXD];
        double velo1,velo2,velo_tmp1,velo_tmp2,velo;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l=-1,m=-1;
        static double lambda1,lambda2;
	SAMPLE *sample = front->sample;
	char *sample_type = sample->sample_type;
	double *sample_line = sample->sample_coords;
	char *out_name = front-> out_name;
	double dens;

	if (front->step < sample->start_step || front->step > sample->end_step)
	    return;
	if ((front->step - sample->start_step)%sample->step_interval)
	    return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
        }
        switch (sample_type[0])
        {
        case 'x':
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index3d(l,0,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[0]);
                --l;
                index = d_index3d(l,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index3d(l+1,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda1 = (sample_line[0] - x1) / (x2 - sample_line[0]);
            }

            switch (sample_type[1])
            {
                case 'y':
                    if (m == -1)
                    {
                        double y1,y2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,m,0,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[1]);
                        --m;
                        index = d_index3d(0,m,0,top_gmax);
                        getRectangleCenter(index,coords);
                        y1 = coords[1];
                        index = d_index3d(0,m+1,0,top_gmax);
                        getRectangleCenter(index,coords);
                        y2 = coords[1];
                        lambda2 = (sample_line[1] - y1)/(y2 - sample_line[1]);
                    }
                    i = l;
                    j = m;
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (k = imin[2]; k <= imax[2]; ++k)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
                    }
                    fclose(sfile);

                    printf("sample line: x = %20.14f, y = %20.14f\n",coords[0],
                        coords[1]);

                    break;

                case 'z':
                    if (m == -1)
                    {
                        double z1,z2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,0,m,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[2]);
                        --m;
                        index = d_index3d(0,0,m,top_gmax);
                        getRectangleCenter(index,coords);
                        z1 = coords[2];
                        index = d_index3d(0,0,m+1,top_gmax);
                        getRectangleCenter(index,coords);
                        z2 = coords[2];
                        lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
                    }
                    i = l;
                    k = m;
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    sprintf(sname,"%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (j = imin[1]; j <= imax[1]; ++j)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i+1,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
                    }
                    fclose(sfile);

                    printf("sample line: x = %20.14f, z = %20.14f\n",coords[0],
                        coords[2]);

                    break;

                    default:
                        printf("Incorrect input for sample velocity!\n");
                        break;

            }
            break;

        case 'y':
            if (l == -1)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index3d(0,l,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[1]);
                --l;
                index = d_index3d(0,l,0,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index3d(0,l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
                lambda1 = (sample_line[0] - y1)/(y2 - sample_line[0]);
            }

            switch (sample_type[1])
            {
                case 'z':
                    if (m == -1)
                    {
                        double z1,z2;
                        do
                        {
                            ++m;
                            index = d_index3d(0,0,m,top_gmax);
                            getRectangleCenter(index,coords);
                        }while(sample_line[1]>=coords[2]);
                        --m;
                        index = d_index3d(0,0,m,top_gmax);
                        getRectangleCenter(index,coords);
                        z1 = coords[2];
                        index = d_index3d(0,0,m+1,top_gmax);
                        getRectangleCenter(index,coords);
                        z2 = coords[2];
                        lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
                    }
                    j = l;
                    k = m;
                    sprintf(sname, "%s/x-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[0][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[0][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/y-%d-%d.xg",out_name,step,count);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[1][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[1][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    sprintf(sname, "%s/z-%d-%d.xg",out_name,step,count++);
                    sfile = fopen(sname,"w");
                    for (i = imin[0]; i <= imax[0]; ++i)
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i,j+1,k,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        index = d_index3d(i,j,k+1,top_gmax);
                        dens = field.dens[index];
                        velo1 = field.momn[2][index]/dens;
                        index = d_index3d(i,j+1,k+1,top_gmax);
                        dens = field.dens[index];
                        velo2 = field.momn[2][index]/dens;
                        velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                        velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                        getRectangleCenter(index,coords);
                        fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
                    }
                    fclose(sfile);

                    printf("sample line: y = %20.14f, z = %20.14f\n",coords[1],
                        coords[2]);

                    break;

                default:
                    printf("Incorrect input for sample velocity!\n");
                    break;
            }
        default:
            printf("Incorrect input for sample velocity!\n");
            break;
        }
}	/* end sampleVelocity3d */

void G_CARTESIAN::sampleVelocity2d()
{
	int i,j,index;
	SAMPLE *sample = front->sample;
        char *sample_type = sample->sample_type;
        double *line = sample->sample_coords;
        char *out_name = front->out_name;
        double coords[MAXD];
        double velox1,velox2,velox;
        double veloy1,veloy2,veloy;
        char dirname[200],sname[200];
        static int count = 0;
        static int step = 0;
        static int l = -1;
        static double lambda;
	double dens1,dens2,c1,c2,c,mach;
	boolean data_in_domain;
	double *L = front->rect_grid->L;
	double *U = front->rect_grid->U;
	EOS_PARAMS *eos;
	STATE state1,state2;
	COMPONENT comp1,comp2;
	double *crds,*vx,*vy,*mc;
	int n,size;

	if (front->step < sample->start_step || front->step > sample->end_step)
            return;
        if ((front->step - sample->start_step)%sample->step_interval)
            return;
        if (step != front->step)
            step = front->step;
	
	state1.dim = state2.dim = 2;
	sprintf(dirname,"%s/samples/sample-%s",
			out_name,right_flush(front->step,6));
	if (!create_directory(dirname,NO))
	{
            screen("Cannot create directory %s\n",dirname);
	    clean_up(ERROR);
	}
        switch (sample_type[0])
        {
        case 'x':
	    if (line[0] < L[0] || line[0] >= U[0])
		data_in_domain = NO;
	    else
		data_in_domain = YES;
            if (data_in_domain == YES)
            {
                double x1,x2;
		for (l = imin[0]; l <= imax[0]; ++l)
		{
                    index = d_index2d(l,0,top_gmax);
                    getRectangleCenter(index, coords);
		    if (line[0] >= coords[0])
			break;
                }
                --l;
                index = d_index2d(l,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index2d(l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda = (line[0] - x1) / (x2 - line[0]);
            	i = l;
	    	size = imax[1] - imin[1] + 1;
	    	FT_VectorMemoryAlloc((POINTER*)&crds,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&vx,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&vy,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&mc,size,sizeof(double));
            	for (j = imin[1]; j <= imax[1]; ++j)
            	{
                    index = d_index2d(i,j,top_gmax);
                    getRectangleCenter(index,coords);
		    comp1 = top_comp[index];
		    state1.eos = &eqn_params->eos[comp1];
		    state1.dens = field.dens[index];
		    state1.engy = field.engy[index];
		    state1.pres = field.pres[index];
		    state1.momn[0] = field.momn[0][index];
		    state1.momn[1] = field.momn[1][index];
		    c1 = EosSoundSpeed(&state1);
		    dens1 = field.dens[index];
                    velox1 = field.momn[0][index]/dens1;
                    veloy1 = field.momn[1][index]/dens1;

                    index = d_index2d(i+1,j,top_gmax);
		    comp2 = top_comp[index];
		    state2.eos = &eqn_params->eos[comp2];
		    state2.dens = field.dens[index];
		    state2.engy = field.engy[index];
		    state2.pres = field.pres[index];
		    state2.momn[0] = field.momn[0][index];
		    state2.momn[1] = field.momn[1][index];
		    c2 = EosSoundSpeed(&state2);
		    dens2 = field.dens[index];
                    velox2 = field.momn[0][index]/dens2;
                    veloy2 = field.momn[1][index]/dens2;

		    crds[n] = coords[1];
		    if (gas_comp(comp1) && gas_comp(comp2))
		    {
                    	vx[n] = (velox1 + lambda*velox2) / (1.0 + lambda);
                    	vy[n] = (veloy1 + lambda*veloy2) / (1.0 + lambda);
                    	c = (c1 + lambda*c1) / (1.0 + lambda);
		    	mc[n] = sqrt(sqr(vx[n]) + sqr(vy[n]))/c;
		    }
		    else if (gas_comp(comp1))
		    {
                    	vx[n] = velox1;
                    	vy[n] = veloy1;
                    	c = c1;
		    	mc[n] = sqrt(sqr(vx[n]) + sqr(vy[n]))/c;
		    }
		    else if (gas_comp(comp2))
		    {
                    	vx[n] = velox2;
                    	vy[n] = veloy2;
                    	c = c2;
		    	mc[n] = sqrt(sqr(vx[n]) + sqr(vy[n]))/c;
		    }
		    else
		    {
		    	vx[n] = vy[n] = mc[n] = 0.0;
		    }
		    n++;
            	}
            }
            sprintf(sname, "vertical-vx-%d.xg",count);
	    FT_XgraphSampleLine(dirname,sname,data_in_domain,size,crds,vx);
            sprintf(sname, "vertical-vy-%d.xg",count);
	    FT_XgraphSampleLine(dirname,sname,data_in_domain,size,crds,vy);
            sprintf(sname, "vertical-mach-%d.xg",count);
	    FT_XgraphSampleLine(dirname,sname,data_in_domain,size,crds,mc);
            break;
        case 'y':
	    if (line[0] < L[1] || line[0] >= U[1])
		data_in_domain = NO;
	    else
		data_in_domain = YES;
            if (data_in_domain == YES)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index2d(0,l,top_gmax);
                    getRectangleCenter(index, coords);
                } while (line[0] >= coords[1]);
                --l;
                index = d_index2d(0,l,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index2d(0,l+1,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
                lambda = (line[0] - y1) / (y2 - line[0]);
            	j = l;
	    	size = imax[0] - imin[0] + 1;
	    	FT_VectorMemoryAlloc((POINTER*)&crds,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&vx,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&vy,size,sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&mc,size,sizeof(double));
		n = 0;
            	for (i = imin[0]; i <= imax[0]; ++i)
            	{
                    index = d_index2d(i,j,top_gmax);
                    getRectangleCenter(index,coords);
		    comp1 = top_comp[index];
		    state1.eos = &eqn_params->eos[comp1];
		    state1.dens = field.dens[index];
		    state1.engy = field.engy[index];
		    state1.pres = field.pres[index];
		    state1.momn[0] = field.momn[0][index];
		    state1.momn[1] = field.momn[1][index];
		    c1 = EosSoundSpeed(&state1);
		    dens1 = field.dens[index];
                    velox1 = field.momn[0][index]/dens1;
                    veloy1 = field.momn[1][index]/dens1;

                    index = d_index2d(i,j+1,top_gmax);
		    comp2 = top_comp[index];
		    state2.eos = &eqn_params->eos[comp2];
		    state2.dens = field.dens[index];
		    state2.engy = field.engy[index];
		    state2.pres = field.pres[index];
		    state2.momn[0] = field.momn[0][index];
		    state2.momn[1] = field.momn[1][index];
		    c2 = EosSoundSpeed(&state2);
		    dens2 = field.dens[index];
                    velox2 = field.momn[0][index]/dens2;
                    veloy2 = field.momn[1][index]/dens2;

		    crds[n] = coords[0];
		    if (gas_comp(comp1) && gas_comp(comp2))
		    {
                    	vx[n] = (velox1 + lambda*velox2) / (1.0 + lambda);
                    	vy[n] = (veloy1 + lambda*veloy2) / (1.0 + lambda);
                    	c = (c1 + lambda*c1) / (1.0 + lambda);
		    	mc[n] = sqrt(sqr(vx[n]) + sqr(vy[n]))/c;
		    }
		    else if (gas_comp(comp1))
		    {
                    	vx[n] = velox1;
                    	vy[n] = veloy1;
                    	c = c1;
		    	mc[n] = sqrt(sqr(vx[n]) + sqr(vy[n]))/c;
		    }
		    else if (gas_comp(comp2))
		    {
                    	vx[n] = velox2;
                    	vy[n] = veloy2;
                    	c = c2;
		    	mc[n] = sqrt(sqr(vx[n]) + sqr(vy[n]))/c;
		    }
		    else
		    {
		    	vx[n] = vy[n] = mc[n] = 0.0;
		    }
		    n++;
            	}
            }
            sprintf(sname, "horizontal-vx-%d.xg",count);
	    FT_XgraphSampleLine(dirname,sname,data_in_domain,size,crds,vx);
            sprintf(sname, "horizontal-vy-%d.xg",count);
	    FT_XgraphSampleLine(dirname,sname,data_in_domain,size,crds,vy);
            sprintf(sname, "horizontal-mach-%d.xg",count);
	    FT_XgraphSampleLine(dirname,sname,data_in_domain,size,crds,mc);
            break;
        }
}	/* end sampleVelocity2d */

void G_CARTESIAN::numericalFlux(
	POINTER scheme_params,
	SWEEP *sweep,
	FSWEEP *fsweep,
	int n)
{
	switch (eqn_params->num_scheme)
	{
	case TVD_FIRST_ORDER:
	case TVD_SECOND_ORDER:
	case TVD_FOURTH_ORDER:
	    TVD_flux(scheme_params,sweep,fsweep,n);
	    break;
	case WENO_FIRST_ORDER:
	case WENO_SECOND_ORDER:
	case WENO_FOURTH_ORDER:
	    WENO_flux(scheme_params,sweep,fsweep,n);
	    break;
	default:
	    (void) printf("Unknow numerical scheme\n");
	    clean_up(ERROR);
	}
}	/* numericalFlux */


void G_CARTESIAN::scatMeshVst(SWEEP *m_vst)
{
	int i,j,k,l,index;

	FT_ParallelExchGridArrayBuffer(m_vst->dens,front,NULL);
	FT_ParallelExchGridArrayBuffer(m_vst->engy,front,NULL);
	FT_ParallelExchGridArrayBuffer(m_vst->pres,front,NULL);
	FT_ParallelExchGridVectorArrayBuffer(m_vst->momn,front);
	/*
	switch (dim)
	{
	case 1:

	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->dens[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->dens[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->engy[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->engy[index] = array[index];
	    }
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		array[index] = m_vst->pres[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index1d(i,top_gmax);
		m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index1d(i,top_gmax);
		    array[index] = m_vst->momn[l][index];
	    	}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index1d(i,top_gmax);
		    m_vst->momn[l][index] = array[index];
	    	}
	    }
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->dens[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->engy[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->engy[index] = array[index];
	    }
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		array[index] = m_vst->pres[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index2d(i,j,top_gmax);
		m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index2d(i,j,top_gmax);
		    array[index] = m_vst->momn[l][index];
	    	}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
                    index = d_index2d(i,j,top_gmax);
		    m_vst->momn[l][index] = array[index];
	    	}
	    }
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->dens[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
		index = d_index3d(i,j,k,top_gmax);
                m_vst->dens[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->engy[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                m_vst->engy[index] = array[index];
	    }
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
                array[index] = m_vst->pres[index];
	    }
	    FT_ParallelExchGridArrayBuffer(array,front,NULL);
            for (k = 0; k <= top_gmax[2]; k++)
            for (j = 0; j <= top_gmax[1]; j++)
            for (i = 0; i <= top_gmax[0]; i++)
            {
                index  = d_index3d(i,j,k,top_gmax);
                m_vst->pres[index] = array[index];
	    }
	    for (l = 0; l < dim; ++l)
	    {
	    	for (k = imin[2]; k <= imax[2]; ++k)
	    	for (j = imin[1]; j <= imax[1]; ++j)
	    	for (i = imin[0]; i <= imax[0]; ++i)
	    	{
		    index = d_index3d(i,j,k,top_gmax);
                    array[index] = m_vst->momn[l][index];
	    	}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
            	for (k = 0; k <= top_gmax[2]; k++)
            	for (j = 0; j <= top_gmax[1]; j++)
            	for (i = 0; i <= top_gmax[0]; i++)
            	{
		    index = d_index3d(i,j,k,top_gmax);
                    m_vst->momn[l][index] = array[index];
	    	}
	    }
	}
	*/
}	/* end scatMeshStates */

void G_CARTESIAN::copyMeshVst(
	const SWEEP& m_vst_orig,
	SWEEP *m_vst)
{
	int i,j,k,l,index;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		m_vst->dens[index] = m_vst_orig.dens[index];
		m_vst->engy[index] = m_vst_orig.engy[index];
		m_vst->pres[index] = m_vst_orig.pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = m_vst_orig.momn[l][index];
	    }
	}
}	/* end copyMeshVst */

void G_CARTESIAN::copyToMeshVst(
	SWEEP *m_vst)
{
	int i,j,k,l,index;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		m_vst->dens[index] = dens[index];
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		m_vst->dens[index] = dens[index];
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		m_vst->dens[index] = dens[index];
		m_vst->engy[index] = engy[index];
		m_vst->pres[index] = pres[index];
		for (l = 0; l < dim; ++l)
		    m_vst->momn[l][index] = momn[l][index];
	    }
	}
}	/* end copyToMeshVst */

void G_CARTESIAN::copyFromMeshVst(
	const SWEEP& m_vst)
{
	int i,j,k,l,index;
	STATE state;
	COMPONENT comp;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	
	//GFM
	if(eqn_params->tracked)
	{
	    start_clock("get_ghost_state");
	    get_ghost_state(m_vst, 2, 0);//GAS_COMP1 = 2
	    get_ghost_state(m_vst, 3, 1);//GAS_COMP2 = 3
	    scatMeshGhost();
	    stop_clock("get_ghost_state");
	}

	state.dim = dim;
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 2:
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		state.dens = m_vst.dens[index];
		state.engy = m_vst.engy[index];
		state.pres = m_vst.pres[index];
		for (l = 0; l < dim; ++l)
		    state.momn[l] = m_vst.momn[l][index];
		if (gas_comp(top_comp[index]))
		{
		    state.eos = &(eqn_params->eos[comp]);
		    checkCorrectForTolerance(&state);
		}
		dens[index] = state.dens;
		engy[index] = state.engy;
		pres[index] = state.pres;
		for (l = 0; l < dim; ++l)
		    momn[l][index] = state.momn[l];
	    }
	}
}	/* end copyFromMeshVst */

void G_CARTESIAN::appendStencilBuffer2d(
	SWEEP *vst,
	SWEEP *m_vst,
	int i,
	int dir)
{
	int		i1,i2,k,offset,index0,index;
	INTERFACE 	*intfc = front->interf;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	double 		crx_coords[MAXD];
	STATE 		*state;
	int		comp, icoords[3];
	INTERFACE	*grid_intfc = front->grid_intfc;

	switch (dir)
	{
	case 0:
	    i2 = i;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY ||
                rect_boundary_type(intfc,dir,0) == ELASTIC_BOUNDARY)
	    {
		i1 = imin[0];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imin[0] + k - 1;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = -m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imin[0] - k;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = m_vst->momn[0][index];
		    vst->momn[1][3-k] = m_vst->momn[1][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else
	    {
		i1 = imin[0];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			ldir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,0,0,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d: "
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }

	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY ||
                rect_boundary_type(intfc,dir,1) == ELASTIC_BOUNDARY)
	    {
		i1 = imax[0];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imax[0] - k + 1;
		    offset = imax[0] - imin[0] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = -m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (k = 1; k <= 3; ++k)
		{
		    i1 = imax[0] + k;
		    index = d_index2d(i1,i2,top_gmax);
		    offset = imax[0] - imin[0] + 3;
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = m_vst->momn[0][index];
		    vst->momn[1][offset+k] = m_vst->momn[1][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else
	    {
		i1 = imax[0];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			rdir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[0]), 1);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    offset = imax[dir] - imin[dir] + nrad;
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,1,
					offset,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }
	    break;
	case 1:
	    i1 = i;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY ||
                rect_boundary_type(intfc,dir,0) == ELASTIC_BOUNDARY)
	    {
		i2 = imin[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imin[1] + k - 1;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = -m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		i2 = imin[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imin[1] - k;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[3-k] = m_vst->dens[index];
		    vst->engy[3-k] = m_vst->engy[index];
		    vst->pres[3-k] = m_vst->pres[index];
		    vst->momn[0][3-k] = m_vst->momn[1][index];
		    vst->momn[1][3-k] = m_vst->momn[0][index];
		    vst->momn[2][3-k] = 0.0;
		}
	    }
	    else
	    {
		i2 = imin[1];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			ldir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,0,0,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }

	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY ||
                rect_boundary_type(intfc,dir,1) == ELASTIC_BOUNDARY)
	    {
		i2 = imax[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imax[1] - k + 1;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = -m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = m_vst->momn[2][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		i2 = imax[1];
		index0 = d_index2d(i1,i2,top_gmax);
		for (k = 1; k <= 3; ++k)
		{
		    i2 = imax[1] + k;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index2d(i1,i2,top_gmax);
		    vst->dens[offset+k] = m_vst->dens[index];
		    vst->engy[offset+k] = m_vst->engy[index];
		    vst->pres[offset+k] = m_vst->pres[index];
		    vst->momn[0][offset+k] = m_vst->momn[1][index];
		    vst->momn[1][offset+k] = m_vst->momn[0][index];
		    vst->momn[2][offset+k] = 0.0;
		}
	    }
	    else
	    {
		i2 = imax[1];
		index = d_index2d(i1,i2,top_gmax);
		comp = top_comp[index];
		icoords[0] = i1;
		icoords[1] = i2;
		if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
			rdir[dir],comp,(POINTER*)&state,&hs,crx_coords))
		{
		    printf("In appendStencilBuffer2d()\n");
		    printf("ERROR: No crossing found!\n");
		    print_int_vector("icoords=", icoords, 2, "\n");
		    printf("direction: %s side %d\n",
		           grid_direction_name(ldir[dir]), 0);
		    clean_up(ERROR);
		}
		switch (wave_type(hs))
		{
		case DIRICHLET_BOUNDARY:
		    offset = imax[dir] - imin[dir] + nrad;
		    setDirichletStates(state,vst,m_vst,hs,icoords,dir,1,
					offset,1);
		    break;
		default: 
		    printf("ERROR appendStencilBuffer2d"
		    	   "unknown boundary type %d\n", wave_type(hs));
		    clean_up(ERROR);
		}
	    }
	}

}	/* end appendStencilBuffer2d */

void G_CARTESIAN::appendStencilBuffer3d(
	SWEEP *vst,
	SWEEP *m_vst,
	int i1,
	int i2,
	int dir)
{
	int i,j,k,l,offset,index;
	INTERFACE *intfc = front->interf;

	switch (dir)
	{
	case 0:
	    j = i1;	k = i2;
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    i = imin[0] - l;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[0][index];
		    vst->momn[1][3-l] = m_vst->momn[1][index];
		    vst->momn[2][3-l] = m_vst->momn[2][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    i = imax[0] + l;
		    index = d_index3d(i,j,k,top_gmax);
		    offset = imax[0] - imin[0] + 3;
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[0][index];
		    vst->momn[1][offset+l] = m_vst->momn[1][index];
		    vst->momn[2][offset+l] = m_vst->momn[2][index];
		}
	    }
	    break;
	case 1:
	    k = i1;	i = i2;
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imin[1] - l;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[1][index];
		    vst->momn[1][3-l] = m_vst->momn[2][index];
		    vst->momn[2][3-l] = m_vst->momn[0][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imin[1] + l - 1;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = -m_vst->momn[1][index];
		    vst->momn[1][3-l] = m_vst->momn[2][index];
		    vst->momn[2][3-l] = m_vst->momn[0][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imax[1] + l;
		    index = d_index3d(i,j,k,top_gmax);
		    offset = imax[1] - imin[1] + 3;
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[1][index];
		    vst->momn[1][offset+l] = m_vst->momn[2][index];
		    vst->momn[2][offset+l] = m_vst->momn[0][index];
		}
	    }
	    else if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    j = imax[1] - l + 1;
		    offset = imax[1] - imin[1] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = -m_vst->momn[1][index];
		    vst->momn[1][offset+l] = m_vst->momn[2][index];
		    vst->momn[2][offset+l] = m_vst->momn[0][index];
		}
	    }
	    break;
	case 2:
	    i = i1;	j = i2;
	    if (rect_boundary_type(intfc,dir,0) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imin[2] + l - 1;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = -m_vst->momn[2][index];
		    vst->momn[1][3-l] = m_vst->momn[0][index];
		    vst->momn[2][3-l] = m_vst->momn[1][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == NEUMANN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imax[2] - l + 1;
		    offset = imax[2] - imin[2] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = -m_vst->momn[2][index];
		    vst->momn[1][offset+l] = m_vst->momn[0][index];
		    vst->momn[2][offset+l] = m_vst->momn[1][index];
		}
	    }
	    
	    if (rect_boundary_type(intfc,dir,0) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imin[2] - l; 
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[3-l] = m_vst->dens[index];
		    vst->engy[3-l] = m_vst->engy[index];
		    vst->pres[3-l] = m_vst->pres[index];
		    vst->momn[0][3-l] = m_vst->momn[2][index];
		    vst->momn[1][3-l] = m_vst->momn[0][index];
		    vst->momn[2][3-l] = m_vst->momn[1][index];
		}
	    }
	    if (rect_boundary_type(intfc,dir,1) == SUBDOMAIN_BOUNDARY)
	    {
		for (l = 1; l <= 3; ++l)
		{
		    k = imax[2] + l;
		    offset = imax[2] - imin[2] + 3;
		    index = d_index3d(i,j,k,top_gmax);
		    vst->dens[offset+l] = m_vst->dens[index];
		    vst->engy[offset+l] = m_vst->engy[index];
		    vst->pres[offset+l] = m_vst->pres[index];
		    vst->momn[0][offset+l] = m_vst->momn[2][index];
		    vst->momn[1][offset+l] = m_vst->momn[0][index];
		    vst->momn[2][offset+l] = m_vst->momn[1][index];
		}
	    }
	}
}	/* end appendStencilBuffer3d */

void G_CARTESIAN::scatMeshStates()
{
	SWEEP vst;
	allocMeshVst(&vst);
	copyToMeshVst(&vst);
	scatMeshVst(&vst);
	copyFromMeshVst(vst);
	freeVst(&vst);
}	/* end scatMeshStates */

void G_CARTESIAN::freeVst(
	SWEEP *vst)
{
	FT_FreeThese(4,vst->dens,vst->engy,vst->pres,vst->momn);
}	/* end freeVstFlux */

void G_CARTESIAN::freeFlux(
	FSWEEP *flux)
{
	FT_FreeThese(3,flux->dens_flux,flux->engy_flux,flux->momn_flux);
}

void G_CARTESIAN::addMeshFluxToVst(
	SWEEP *m_vst,
	const FSWEEP& m_flux,
	double chi)
{
	int 		i,j,k,l,index;
	double		ke,c,u;
	EOS_PARAMS	*eos;
	STATE		st;
	int		comp;
	double		temp;

	switch (dim)
	{
	case 1:
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
		comp = top_comp[index];
		if (!gas_comp(comp))
		{
		    m_vst->dens[index] = 0.0;
		    m_vst->engy[index] = 0.0;
		    for (l = 0; l < dim; ++l)
		    	m_vst->momn[l][index] = 0.0; 
		    continue;
		}
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
	    scatMeshVst(m_vst);
	    break;
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index2d(i,j,top_gmax);
		comp = top_comp[index];
		if (!gas_comp(comp))
		{
		    m_vst->dens[index] = 0.0;
		    m_vst->engy[index] = 0.0;
		    for (l = 0; l < dim; ++l)
		    	m_vst->momn[l][index] = 0.0; 
		    continue;
		}
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
	    scatMeshVst(m_vst);
	    break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
		index = d_index3d(i,j,k,top_gmax);
		comp = top_comp[index];
		eos = &(eqn_params->eos[comp]);

		m_vst->dens[index] += chi*m_flux.dens_flux[index];
		m_vst->engy[index] += chi*m_flux.engy_flux[index];
		ke = u = 0.0;
		for (l = 0; l < dim; ++l)
		{
		    m_vst->momn[l][index] += chi*m_flux.momn_flux[l][index];
		    ke += sqr(m_vst->momn[l][index]);
		    u += sqr(m_vst->momn[l][index]);
		}
		
		CovertVstToState(&st, m_vst, eos, index, dim);
		checkCorrectForTolerance(&st);
		m_vst->dens[index] = st.dens;
		m_vst->pres[index] = st.pres;
		m_vst->engy[index] = st.engy;
		u = sqrt(u)/m_vst->dens[index];
		c = EosSoundSpeed(&st);
		temp = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
                if (max_speed < temp)
                    max_speed = temp;
	    }
	    scatMeshVst(m_vst);
	}
}	/* end addMeshFluxToVst */

void G_CARTESIAN::appendGhostBuffer(
	SWEEP *vst,
	SWEEP *m_vst,
	int n,
	int *icoords,
	int idir,
	int nb)
{
	int		i,j,k,index,ic[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	HYPER_SURF	*hs1;
	COMPONENT 	comp;
	double 		crx_coords[MAXD];
	STATE 		*state,ghost_st;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic_next[MAXD];
	INTERFACE	*grid_intfc = front->grid_intfc;
	static int count = 0;
	count++;
	boolean Debug = NO;

	SURFACE **s;


	if (debugging("append_buffer"))
		printf("Entering appendGhostBuffer()\n");

	for (i = 0; i < dim; ++i)
        ic[i] = icoords[i];
	
	index = d_index(ic,top_gmax,dim);
	comp = cell_center[index].comp;

	switch(nb)
	{
	case 0:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] - i;
		index = d_index(ic,top_gmax,dim);

//	 	The following is for debugging		    
		boolean status;
//		check neighbor in ldir[idir] 
		for (k = 0; k < dim; ++k)
		    ic_next[k] = ic[k];
		ic_next[idir]++;
		status = FT_StateStructAtGridCrossing(front,grid_intfc,
			ic_next,ldir[idir],comp,(POINTER*)&state,
			&hs,crx_coords);
/*
		if (status)
		if (status && wave_type(hs) != 6)
		    printf("HI,status=%d, wave_type(*hs)=%d\n",status,wave_type(hs));
*/

		if (!needBufferFromIntfc(comp,cell_center[index].comp) && !status)
		{
		    vst->dens[nrad-i] = m_vst->dens[index];
		    vst->engy[nrad-i] = m_vst->engy[index];
		    vst->pres[nrad-i] = m_vst->pres[index];

		    for (j = 0; j < 3; j++)
			vst->momn[j][nrad-i] = 0.0;
		    if (dim == 1)
			vst->momn[0][nrad-i] = m_vst->momn[0][index];
		    else if (dim == 2)
                for(j = 0; j < 2; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3){
			for (j = 0; j < 3; j++){
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind3[idir][j]][index];
			}

		    }

		}
		else
		{
		    if (!status) /* extreme cases */
		    {
			double coords[MAXD], wtol[MAXD], tol[MAXD];
			int ic_tmp[MAXD];
			for (k = 0; k < dim; ++k)
                            tol[k] = 2.0 * IG_TOL * top_h[k];
			/* check neighbor in the opposite direction */
			for (k = 0; k < dim; ++k)
                            ic_tmp[k] = ic[k];
                        ic_tmp[idir]--;
			status = FT_StateStructAtGridCrossing(front,grid_intfc,
					ic_tmp,rdir[idir],comp,(POINTER*)&state,
					&hs,crx_coords);
			if(!status)
			{
			    /* check second neighbor in the same direction */
			    for (k = 0; k < dim; ++k)
                                ic_tmp[k] = ic[k];
                            ic_tmp[idir] += 2;
                            status = FT_StateStructAtGridCrossing(front,
                                        grid_intfc,ic_tmp,ldir[idir],comp,
                                        (POINTER*)&state,&hs,crx_coords);
			    if (!status)
			    {
				/* must be something wrong */
		    	    	printf("In appendGhostBuffer() Case 0\n");
		    	    	printf("ERROR: No crossing found!\n");
		    	    	print_int_vector("ic=",ic,dim,"\n");
		    	    	printf("direction: %s side %d\n",
		           		grid_direction_name(ldir[idir]), nb);
				clean_up(ERROR);
			    }
			    else
			    {
				/* check if crossing is close enough */
			        boolean close_enough = YES;
                                for (k = 0; k < dim; ++k)
                                {
                                    coords[k] = top_L[k]+ic_next[k]*top_h[k];
                                    wtol[k] = crx_coords[k] - coords[k];
                                    if (fabs(wtol[k]) > tol[k])
                                        close_enough = NO;
                                }
                                if (!close_enough)
                                {
                                    (void) printf("ERROR: Not close enough!\n");
                                    clean_up(ERROR);
                                }
			    }
			}
			else
			{
			    /* check if crossing is close enough */
			    boolean close_enough = YES;
			    for (k = 0; k < dim; ++k)
                {
                    coords[k] = top_L[k] + ic[k] * top_h[k];
                    wtol[k] = crx_coords[k] - coords[k];
                    if (fabs(wtol[k]) > tol[k])
                        close_enough = NO;
                }
                if (!close_enough)
                {
                    (void) printf("ERROR: Not close enough!\n");
                    clean_up(ERROR);
			    }
            }
		    
            }

		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,0,i,comp);
		    	break;
		    case ELASTIC_BOUNDARY:
		    	setElasticStates(vst,m_vst,hs,state,ic_next,idir,
						nb,0,i,comp);
			break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,
					idir,nb,0,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	//GFM
		    	GFMGhostState(ic,comp,&ghost_st);
		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[nrad-k] = ghost_st.dens;
		    	    vst->engy[nrad-k] = ghost_st.engy;
		    	    vst->pres[nrad-k] = ghost_st.pres;
			
			    for (j=0; j < 3; j++)
			    	    vst->momn[j][nrad-k] = 0.0;
			    if (dim == 1)
				vst->momn[0][nrad-k] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for (j=0; j < 2; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for (j = 0; j < 3; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
					wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=", icoords,3,"\n");
		    	clean_up(ERROR);
		    }
		    break;
		}
	    
        }
	    break;
	case 1:
	    for (i = 1; i <= nrad; ++i)  //nrad=3
	    {
		ic[idir] = icoords[idir] + i;
		index = d_index(ic,top_gmax,dim);

//		For debugging
		boolean status;
//		check neighbor in rdir[idir] 
		for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		ic_next[idir]--;
		status = FT_StateStructAtGridCrossing(front,grid_intfc,
				icoords,rdir[idir],comp,(POINTER*)&state,
				&hs,crx_coords);

//For the needBufferFromIntfc function, if the two component are different,
//YES is returned. Then for the following, the if statement is satisfied when 
//the two component are the same, which means it does not meet the rectangle
//boundary. It may meet the elastic boundary or does not meet any boundary.
//Then !status exclude the possibility of meeting elastic boundary.


		if (!needBufferFromIntfc(comp,cell_center[index].comp) && !status )
		{
//		    if (status && wave_type(hs) == 13)
//			printf("233: target boundary found.\n");
		    vst->dens[n+nrad+i-1] = m_vst->dens[index];
		    vst->engy[n+nrad+i-1] = m_vst->engy[index];
		    vst->pres[n+nrad+i-1] = m_vst->pres[index];
		    
		    for (j = 0; j < 3; j++)
			vst->momn[j][n+nrad+i-1] = 0.0;
		    if (dim == 1)
			vst->momn[0][n+nrad+i-1] = 
			         	m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    	vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{


		    if (!status) /* extreme cases */
		    {
			double coords[MAXD], wtol[MAXD], tol[MAXD];
			int ic_tmp[MAXD];
			for (k = 0; k < dim; ++k)
                            tol[k] = 2.0 * IG_TOL * top_h[k];
			/* check neighbor in the opposite direction */
			for (k = 0; k < dim; ++k)
                            ic_tmp[k] = ic[k];
                        ic_tmp[idir]++;
			status = FT_StateStructAtGridCrossing(front,grid_intfc,
					ic_tmp,ldir[idir],comp,(POINTER*)&state,
					&hs,crx_coords);
			if(!status)
			{
			    /* check second neighbor in the same direction */
			    for (k = 0; k < dim; ++k)
                                ic_tmp[k] = ic[k];
                            ic_tmp[idir] -= 2;
                            status = FT_StateStructAtGridCrossing(front,
                                        grid_intfc,ic_tmp,rdir[idir],comp,
                                        (POINTER*)&state,&hs,crx_coords);
			    if (!status)
			    {
				/* must be something wrong */
		    	    	printf("In appendGhostBuffer() Case 1\n");
		    	    	printf("ERROR: No crossing found!\n");
		    	    	print_int_vector("ic=",ic,dim,"\n");
		    	    	printf("direction: %s side %d\n",
		           		grid_direction_name(ldir[idir]), nb);
				clean_up(ERROR);
			    }
			    else
			    {
				/* check if crossing is close enough */
			        boolean close_enough = YES;
                                for (k = 0; k < dim; ++k)
                                {
                                    coords[k] = top_L[k] + ic_next[k]*top_h[k];
                                    wtol[k] = crx_coords[k] - coords[k];
                                    if (fabs(wtol[k]) > tol[k])
                                        close_enough = NO;
                                }
                                if (!close_enough)
                                {
                                    (void) printf("ERROR: Not close enough!\n");
                                    clean_up(ERROR);
                                }
			    }
			}
			else
			{
			    /* check if crossing is close enough */
			    boolean close_enough = YES;
			    for (k = 0; k < dim; ++k)
                            {
                                coords[k] = top_L[k] + ic[k] * top_h[k];
                                wtol[k] = crx_coords[k] - coords[k];
			        if (fabs(wtol[k]) > tol[k])
				    close_enough = NO;
			    }
			    if (!close_enough)
			    {
			        (void) printf("ERROR: Not close enough!\n");
			        clean_up(ERROR);
			    }
			}
		    }

		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,n,i,comp);
		    	break;
		    case ELASTIC_BOUNDARY:
		    	setElasticStates(vst,m_vst,hs,state,ic_next,idir,
						nb,n,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,idir,nb,
						n,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	//GFM
		    	GFMGhostState(ic,comp,&ghost_st);

		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[n+nrad+k-1] = ghost_st.dens;
		    	    vst->engy[n+nrad+k-1] = ghost_st.engy;
		    	    vst->pres[n+nrad+k-1] = ghost_st.pres;
			
			    for(j=0; j<3; j++)
			    	vst->momn[j][n+nrad+k-1] = 0.0;
			    if (dim == 1)
				vst->momn[0][n+nrad+k-1] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for(j = 0; j < 2; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for(j = 0; j < 3; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
				wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=",icoords,3,"\n");
		    	(void) printf("nb = %d\n",nb);
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	}
}	/* end appendGhostBuffer */

/*
void G_CARTESIAN::appendGhostBuffer(
	SWEEP *vst,
	SWEEP *m_vst,
	int n,
	int *icoords,
	int idir,
	int nb)
{
	int		i,j,k,index,ic[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	COMPONENT 	comp;
	double 		crx_coords[MAXD];
	STATE 		*state,ghost_st;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic_next[MAXD];
	INTERFACE	*grid_intfc = front->grid_intfc;

	if (debugging("append_buffer"))
		printf("Entering appendGhostBuffer()\n");

	for (i = 0; i < dim; ++i) ic[i] = icoords[i];
	
	index = d_index(ic,top_gmax,dim);
	comp = cell_center[index].comp;

	switch(nb)
	{
	case 0:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] - i;
		index = d_index(ic,top_gmax,dim);
		    
		if (!needBufferFromIntfc(comp,cell_center[index].comp))
		{
		    vst->dens[nrad-i] = m_vst->dens[index];
		    vst->engy[nrad-i] = m_vst->engy[index];
		    vst->pres[nrad-i] = m_vst->pres[index];

		    for (j = 0; j < 3; j++)
			vst->momn[j][nrad-i] = 0.0;
		    if (dim == 1)
			vst->momn[0][nrad-i] = m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{
		    boolean status;
		    // check neighbor in ldir[idir]
		    for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		    ic_next[idir]++;
		    status = FT_StateStructAtGridCrossing(front,grid_intfc,
				ic_next,ldir[idir],comp,(POINTER*)&state,
				&hs,crx_coords);
		    if (!status) // extreme cases
		    {
			double coords[MAXD], wtol[MAXD], tol[MAXD];
			int ic_tmp[MAXD];
			for (k = 0; k < dim; ++k)
                            tol[k] = 2.0 * IG_TOL * top_h[k];
			// check neighbor in the opposite direction
			for (k = 0; k < dim; ++k)
                            ic_tmp[k] = ic[k];
                        ic_tmp[idir]--;
			status = FT_StateStructAtGridCrossing(front,grid_intfc,
					ic_tmp,rdir[idir],comp,(POINTER*)&state,
					&hs,crx_coords);
			if(!status)
			{
			    // check second neighbor in the same direction
			    for (k = 0; k < dim; ++k)
                                ic_tmp[k] = ic[k];
                            ic_tmp[idir] += 2;
                            status = FT_StateStructAtGridCrossing(front,
                                        grid_intfc,ic_tmp,ldir[idir],comp,
                                        (POINTER*)&state,&hs,crx_coords);
			    if (!status)
			    {
				// must be something wrong
		    	    	printf("In appendGhostBuffer() Case 0\n");
		    	    	printf("ERROR: No crossing found!\n");
		    	    	print_int_vector("ic=",ic,dim,"\n");
		    	    	printf("direction: %s side %d\n",
		           		grid_direction_name(ldir[idir]), nb);
				clean_up(ERROR);
			    }
			    else
			    {
				// check if crossing is close enough
			        boolean close_enough = YES;
                                for (k = 0; k < dim; ++k)
                                {
                                    coords[k] = top_L[k]+ic_next[k]*top_h[k];
                                    wtol[k] = crx_coords[k] - coords[k];
                                    if (fabs(wtol[k]) > tol[k])
                                        close_enough = NO;
                                }
                                if (!close_enough)
                                {
                                    (void) printf("ERROR: Not close enough!\n");
                                    clean_up(ERROR);
                                }
			    }
			}
			else
			{
			    // check if crossing is close enough
			    boolean close_enough = YES;
			    for (k = 0; k < dim; ++k)
                            {
                                coords[k] = top_L[k] + ic[k] * top_h[k];
                                wtol[k] = crx_coords[k] - coords[k];
                                if (fabs(wtol[k]) > tol[k])
                                    close_enough = NO;
                            }
                            if (!close_enough)
                            {
                                (void) printf("ERROR: Not close enough!\n");
                                clean_up(ERROR);
			    }
                        }
		    }
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,0,i,comp);
		    	break;
            case ELASTIC_BOUNDARY:
                setElasticStates(vst,m_vst,hs,state,ic_next,idir,nb,n,i,comp);
                break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,
					idir,nb,0,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	//GFM
                        if (debugging("append_buffer"))
                            printf("Calling GFMGhostState()\n");
		    	GFMGhostState(ic,comp,&ghost_st);
		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[nrad-k] = ghost_st.dens;
		    	    vst->engy[nrad-k] = ghost_st.engy;
		    	    vst->pres[nrad-k] = ghost_st.pres;
			
			    for (j=0; j < 3; j++)
			    	    vst->momn[j][nrad-k] = 0.0;
			    if (dim == 1)
				vst->momn[0][nrad-k] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for (j=0; j < 2; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for (j = 0; j < 3; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
					wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=", icoords,3,"\n");
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	    break;
	case 1:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] + i;
		index = d_index(ic,top_gmax,dim);
		if (!needBufferFromIntfc(comp,cell_center[index].comp))
		{
		    vst->dens[n+nrad+i-1] = m_vst->dens[index];
		    vst->engy[n+nrad+i-1] = m_vst->engy[index];
		    vst->pres[n+nrad+i-1] = m_vst->pres[index];
		    
		    for (j = 0; j < 3; j++)
			vst->momn[j][n+nrad+i-1] = 0.0;
		    if (dim == 1)
			vst->momn[0][n+nrad+i-1] = 
			         	m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    	vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{
		    boolean status;
		    // check neighbor in rdir[idir]
		    for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		    ic_next[idir]--;
		    status = FT_StateStructAtGridCrossing(front,grid_intfc,
				ic_next,rdir[idir],comp,(POINTER*)&state,
				&hs,crx_coords);
		    if (!status) // extreme cases
		    {
			double coords[MAXD], wtol[MAXD], tol[MAXD];
			int ic_tmp[MAXD];
			for (k = 0; k < dim; ++k)
                            tol[k] = 2.0 * IG_TOL * top_h[k];
			// check neighbor in the opposite direction
			for (k = 0; k < dim; ++k)
                            ic_tmp[k] = ic[k];
                        ic_tmp[idir]++;
			status = FT_StateStructAtGridCrossing(front,grid_intfc,
					ic_tmp,ldir[idir],comp,(POINTER*)&state,
					&hs,crx_coords);
			if(!status)
			{
			    // check second neighbor in the same direction
			    for (k = 0; k < dim; ++k)
                                ic_tmp[k] = ic[k];
                            ic_tmp[idir] -= 2;
                            status = FT_StateStructAtGridCrossing(front,
                                        grid_intfc,ic_tmp,rdir[idir],comp,
                                        (POINTER*)&state,&hs,crx_coords);
			    if (!status)
			    {
				// must be something wrong
		    	    	printf("In appendGhostBuffer() Case 1\n");
		    	    	printf("ERROR: No crossing found!\n");
		    	    	print_int_vector("ic=",ic,dim,"\n");
		    	    	printf("direction: %s side %d\n",
		           		grid_direction_name(ldir[idir]), nb);
				clean_up(ERROR);
			    }
			    else
			    {
				// check if crossing is close enough
			        boolean close_enough = YES;
                                for (k = 0; k < dim; ++k)
                                {
                                    coords[k] = top_L[k] + ic_next[k]*top_h[k];
                                    wtol[k] = crx_coords[k] - coords[k];
                                    if (fabs(wtol[k]) > tol[k])
                                        close_enough = NO;
                                }
                                if (!close_enough)
                                {
                                    (void) printf("ERROR: Not close enough!\n");
                                    clean_up(ERROR);
                                }
			    }
			}
			else
			{
			    // check if crossing is close enough
			    boolean close_enough = YES;
			    for (k = 0; k < dim; ++k)
                            {
                                coords[k] = top_L[k] + ic[k] * top_h[k];
                                wtol[k] = crx_coords[k] - coords[k];
			        if (fabs(wtol[k]) > tol[k])
				    close_enough = NO;
			    }
			    if (!close_enough)
			    {
			        (void) printf("ERROR: Not close enough!\n");
			        clean_up(ERROR);
			    }
			}
		    }
		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,n,i,comp);
		    	break;
            case ELASTIC_BOUNDARY:
                setElasticStates(vst,m_vst,hs,state,ic_next,idir,nb,n,i,comp);
                break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,idir,nb,
						n,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	//GFM
		    	GFMGhostState(ic,comp,&ghost_st);

		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[n+nrad+k-1] = ghost_st.dens;
		    	    vst->engy[n+nrad+k-1] = ghost_st.engy;
		    	    vst->pres[n+nrad+k-1] = ghost_st.pres;
			
			    for(j=0; j<3; j++)
			    	vst->momn[j][n+nrad+k-1] = 0.0;
			    if (dim == 1)
				vst->momn[0][n+nrad+k-1] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for(j = 0; j < 2; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for(j = 0; j < 3; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
				wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=",icoords,3,"\n");
		    	(void) printf("nb = %d\n",nb);
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	}
}*/	/* end appendGhostBuffer */

//ghost fluid method.
void G_CARTESIAN::solve_exp_value()
{
	int		i, j, k, n;
	int		index;

	fflush(NULL);

	double **gnor = eqn_params->gnor;
	get_normal_from_front(); //modifies gnor

	if (dim == 1)
	{
	    for(k=0; k<dim; k++)
	    {
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	 	    index = d_index1d(i,top_gmax);
		    array[index] = gnor[k][index];
		}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index1d(i,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
	else if (dim == 2)
	{
	    for(k=0; k<dim; k++)
	    {
		for (j = imin[1]; j <= imax[1]; ++j)
		for (i = imin[0]; i <= imax[0]; ++i)
		{
	 	    index = d_index2d(i,j,top_gmax);
		    array[index] = gnor[k][index];
		}
	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        	for (j = 0; j <= top_gmax[1]; j++)
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index2d(i,j,top_gmax);
	    	    gnor[k][index] = array[index];
		}
	    }
	}
	else
	{
	    for(k=0; k<dim; k++)
	    {
            for (n = imin[2]; n <= imax[2]; ++n)
            for (j = imin[1]; j <= imax[1]; ++j)
            for (i = imin[0]; i <= imax[0]; ++i)
            {
                    index = d_index3d(i,j,n,top_gmax);
                    array[index] = gnor[k][index];
            }

	    	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        	for (n = 0; n <= top_gmax[2]; n++)
        	for (j = 0; j <= top_gmax[1]; j++)
        	for (i = 0; i <= top_gmax[0]; i++)
        	{
	    	    index  = d_index3d(i,j,n,top_gmax);
	    	    gnor[k][index] = array[index];
		    }
	    }
	}
}

void G_CARTESIAN::scatMeshGhost()
{
	int		i, j, k, n, index;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	static int count = 0;

	count++;

	if(dim == 2)
	{
	for(k=0; k<2; k++)
	{
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gdens[k][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gdens[k][index] = array[index];
	}
	
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gpres[k][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gpres[k][index] = array[index];
	}

	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gvel[k][0][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gvel[k][0][index] = array[index];
	}

	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = Gvel[k][1][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index2d(i,j,top_gmax);
	    Gvel[k][1][index] = array[index];
	}
	}    //for k
	}    //if dim == 2
	else if(dim == 3)
	{
	for(k=0; k<2; k++)
	{
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gdens[k][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gdens[k][index] = array[index];
	}

	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gpres[k][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gpres[k][index] = array[index];
	}
	
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][0][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][0][index] = array[index];
	}

	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][1][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][1][index] = array[index];
	}
	
	for (n = imin[2]; n <= imax[2]; ++n)
	for (j = imin[1]; j <= imax[1]; ++j)
	for (i = imin[0]; i <= imax[0]; ++i)
	{
	    index = d_index3d(i,j,n,top_gmax);
	    array[index] = Gvel[k][2][index];
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
        for (n = 0; n <= top_gmax[2]; n++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
	    index  = d_index3d(i,j,n,top_gmax);
	    Gvel[k][2][index] = array[index];
	}
	}    //for k
	}    //for dim == 3
}

#define	corner_index(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]-0.5))

boolean	find_block(
	double		*f,
	int		*icrds,
	double		*p,
	RECT_GRID	*gr)
{
	int	i;
	int	dim=gr->dim;

	for(i=0; i<dim; i++)
	{
	    icrds[i] = corner_index(p[i],i,gr);
	    if(icrds[i] < -gr->lbuf[i] || icrds[i] >= gr->gmax[i]+gr->ubuf[i]-1)
		return  NO;
	    f[i] = p[i] - (gr->L[i]+(0.5+icrds[i])*gr->h[i]);
	    f[i] /= gr->h[i];
	}
	return  YES;
}

boolean G_CARTESIAN::get_ave_normal(
	int		*ic,
	int		***norset)
{
	double		f;
	int		i, j, k, n, ic1[3], ic2[3], dir, num;
	boolean		found;
	int		index0, index;
	double		**gnor = eqn_params->gnor;
	
	found = NO;

	for(i=0; i<dim; i++)
	for(j=0; j<2; j++)
	{
		dir = j == 0 ? -1 : 1;
		ft_assign(ic1, ic, 3*INT);
		ic1[i] = ic[i] + dir;

		if(ic1[i] < 0 || ic1[i] > top_gmax[i])
		    continue;
		if(norset[ic1[0]][ic1[1]][ic1[2]] == 1)
		    found = YES;
	}
	if(!found)
	    return NO;
	
	index0  = d_index(ic,top_gmax,dim);

	gnor[0][index0] = 0.0;
	if (dim > 1)
	    gnor[1][index0] = 0.0;
	if (dim > 2)
	    gnor[2][index0] = 0.0;

	num = 0;
	for(i=ic[0]-1; i<=ic[0]+1; i++)
	for(j=ic[1]-1; j<=ic[1]+1; j++)
	for(k=ic[2]-1; k<=ic[2]+1; k++)
	{
	    if(i < 0 || i > top_gmax[0] || 
	       j < 0 || j > top_gmax[1] || 
	       k < 0 || k > top_gmax[2]) 
		continue;
	    if(norset[i][j][k] != 1)
		continue;

	    ic2[0] = i;
	    ic2[1] = j;
	    ic2[2] = k;
	    index  = d_index(ic2,top_gmax,dim);
		    
	    //do not use length weighted normal direction
	    gnor[0][index0] += gnor[0][index];
	    if(dim > 1)
	    	gnor[1][index0] += gnor[1][index];
	    if(dim > 2)
		gnor[2][index0] += gnor[2][index];
	    num++;
	}
	
	f = 0.0;
	for(n=0; n<dim; n++)
	    f += sqr(gnor[n][index0]);
	f = sqrt(f);

	if(f < 1.0e-6)
	{
	    gnor[0][index0] = 0.0;
	    if (dim > 1) gnor[1][index0] = 0.0;
	    if (dim > 2) gnor[2][index0] = 0.0;
	}
	else
	{
	    gnor[0][index0] /= f;
	    if (dim > 1) gnor[1][index0] /= f;
	    if (dim > 2) gnor[2][index0] /= f;
	}

	return YES;
}

boolean	find_block(double*,int*,double*,RECT_GRID*);

void	get_normal_from_front();

//it will fill gnor field by interface normals
void G_CARTESIAN::get_normal_from_front()
{
	INTERFACE               *intfc = front->interf;
	RECT_GRID		*rgr = front->rect_grid;
	HYPER_SURF              *hs;
	HYPER_SURF_ELEMENT      *hse;
	POINT                   *p;
	int			i,j,k,n, num;
	int			ic[3];
	double			curv,nor[3],d[3],f,d1,d2,d3,*pt,tol;
	int			ix, iy, iz, index;
	boolean			found;
	double			**gnor = eqn_params->gnor;
	static	int		***norset;
	int ict[3];
	double ptt[3];
	boolean status;

	if (norset == NULL)
	    FT_TriArrayMemoryAlloc((POINTER*)&norset,top_gmax[0]+1,
				top_gmax[1]+1,top_gmax[2]+1,INT);

	tol = hmin*1.0e-6;

	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k <= top_gmax[2]; k++)
	{
	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = k;
	    index = d_index(ic,top_gmax,dim);
	    gnor[0][index] = 0.0;
	    if (dim > 1)
	    	gnor[1][index] = 0.0;
	    if (dim > 2)
		gnor[2][index] = 0.0;
	    norset[i][j][k] = 0;
	}

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (Boundary_point(p))
	    {
		p->_nor[0] = 0.0;
		p->_nor[1] = 0.0;
		p->_nor[2] = 0.0;
		continue;
	    }

	    normal(p,hse,hs,nor,front);
	    curv = p->curvature;
            if (the_point(p))
            {
                printf("Testing normal of point: %f %f %f\n",
                        Coords(p)[0],Coords(p)[1],Coords(p)[2]);
                printf("nor = %f %f %f curv = %f\n",nor[0],nor[1],nor[2],curv);
            }

	    pt = Coords(p);
	   
	    status = rect_in_which(pt,ict,top_grid);
	    if (!status) continue;
	    for(i = 0; i < dim; i++)
		ptt[i] = top_grid->L[i] + ict[i]*top_grid->h[i];

	    for(i = 0; i < dim; i++)
	    {
		d[i] = fabs(pt[i]-ptt[i])/rgr->h[i];
	        if(d[i] < -tol || d[i] > 1.0 + tol)
		{
		    status = NO;
		}
	    }
	    if (status == NO) continue;

	    if (dim == 1)
	    {
	    	for(i = 0; i < 2; i++)
		{
		    ix = ict[0] + i;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    f = d1;

		    index = d_index1d(ix,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    norset[ix][0][0] = 1;
		}
	    }
	    else if (dim == 2)
	    {
	    	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
		    ix = ict[0] + i;
		    iy = ict[1] + j;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    d2 = (j == 0) ? fabs(1.0-d[1]) : d[1];
		    f = d1*d2;
		    index = d_index2d(ix,iy,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    gnor[1][index] += nor[1]*f;
		    norset[ix][iy][0] = 1;
		}
	    }
	    else if (dim == 3)
	    {
	    	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
		{
		    ix = ict[0] + i;
		    iy = ict[1] + j;
		    iz = ict[2] + k;
		    d1 = (i == 0) ? fabs(1.0-d[0]) : d[0];
		    d2 = (j == 0) ? fabs(1.0-d[1]) : d[1];
		    d3 = (k == 0) ? fabs(1.0-d[2]) : d[2];
		    f = d1*d2*d3;
		    index = d_index3d(ix,iy,iz,top_gmax);
		    gnor[0][index] += nor[0]*f;
		    gnor[1][index] += nor[1]*f;
		    gnor[2][index] += nor[2]*f;
		    norset[ix][iy][iz] = 1;
		}
	    }
	}

	for (i = 0; i <= top_gmax[0]; i++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (k = 0; k <= top_gmax[2]; k++)
	{
	    //make sure Vel(st) is assigned 
	    if (norset[i][j][k] != 1)
		continue;

	    ic[0] = i;
	    ic[1] = j;
	    ic[2] = k;
	    index = d_index(ic,top_gmax,dim);
	    f = 0.0;
	    for(n=0; n<dim; n++)
		f += sqr(gnor[n][index]);
	    f = sqrt(f);

	    if (f < 1.0e-10)
	    {
		gnor[0][index] = 0.0;
		if (dim > 1) gnor[1][index] = 0.0;
		if (dim > 2) gnor[2][index] = 0.0;
	    }
	    else
	    {
		gnor[0][index] /= f;
		if (dim > 1) gnor[1][index] /= f;
		if (dim > 2) gnor[2][index] /= f;
	    }
	}

	found = YES;
	num = 1;
	while (found && num > 0)
	{
	    found = NO;
	    num = 0;

	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
	    {
		if(norset[i][j][k] != 0)
		    continue;

		found = YES;
		ic[0] = i;
		ic[1] = j;
		ic[2] = k;

		if(get_ave_normal(ic,norset))
		{
		    num++;
		    norset[i][j][k] = 2;
		}
	    }

	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]; k++)
		if(norset[i][j][k] == 2)
		    norset[i][j][k] = 1;
	}
}

void G_CARTESIAN::tecplot_interior_states(
			char	*bname)
{
	char		s[1000];
	double		coords[3];
	int		ix, iy, iz, comp, i, imin[3], imax[3];
	double		**vel = eqn_params->vel;
	double		*dens = eqn_params->dens;
	double		*pres = eqn_params->pres;
	double		**gnor = eqn_params->gnor;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	FILE		*fp;
	int		index;

	sprintf(s,"%s-%d.plt", bname,pp_mynode());
	printf("tecplot_interior_states  file name %s \n",s);

	fp = fopen(s, "w");
	if(fp == NULL)
	{
	    printf("WARNING tecplot_interior_states, can not open file %s\n", s);
	    return; 
	}
	
	for(i=0; i<3; i++)
	{
	    imin[i] = 0;
	    imax[i] = top_gmax[i];
	}

	fprintf(fp, "TITLE = \"inner states\" ");
	if(dim == 2)
	{
	    fprintf(fp, "VARIABLES = \"x\", \"y\", \"comp\",  ");
	    fprintf(fp, "\"dens\", \"press\", \"u\", \"v\",  " );
	    fprintf(fp, "\"nx\", \"ny\",  " );
	    fprintf(fp, "\"dens1\", \"press1\", \"u1\", \"v1\",  " );
	    fprintf(fp, "\"dens2\", \"press2\", \"u2\", \"v2\"  \n" );
	}
	else
	{
	    fprintf(fp, "VARIABLES = \"x\", \"y\", \"z\", \"comp\",  ");
	    fprintf(fp, "\"dens\", \"press\", \"u\", \"v\", \"w\", " );
	    fprintf(fp, "\"nx\", \"ny\", \"nz\", " );
	    fprintf(fp, "\"dens1\", \"press1\", \"u1\", \"v1\", \"w1\" " );
	    fprintf(fp, "\"dens2\", \"press2\", \"u2\", \"v2\"  \"w2\"\n" );
	}

	if(dim == 2)
	    fprintf(fp, "ZONE i=%d, j=%d \n", imax[0]-imin[0]+1, imax[1]-imin[1]+1);
	else
	    fprintf(fp, "ZONE i=%d, j=%d, k=%d \n", imax[0]-imin[0]+1, imax[1]-imin[1]+1, imax[2]-imin[2]+1);

	if(dim == 2)
	{
	    for(iy=imin[1]; iy <= imax[1]; iy++)
		  for(ix=imin[0]; ix <= imax[0]; ix++)
		  {
			index = d_index2d(ix,iy,top_gmax);
			
			getRectangleCenter(index, coords);
			comp = cell_center[index].comp;

			fprintf(fp, "%f ", coords[0]);
			fprintf(fp, "%f ", coords[1]);
			fprintf(fp, "%d ", comp);
		
			if(!gas_comp(comp))
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e   %12.5e %12.5e  ", 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ", 
			    		0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ", 
			    		0.0, 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e   %12.5e %12.5e  ", 
					dens[index], pres[index], 
					vel[0][index],
					vel[1][index],
					gnor[0][index],
					gnor[1][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ",
			    		Gdens[0][index], Gpres[0][index],
					Gvel[0][0][index], Gvel[0][1][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e  ",
			    		Gdens[1][index], Gpres[1][index],
					Gvel[1][0][index], Gvel[1][1][index]);
			}

			fprintf(fp, "\n");
		  }
	}
	else
	{
	for(iz=imin[2]; iz <= imax[2]; iz++)
	    for(iy=imin[1]; iy <= imax[1]; iy++)
		  for(ix=imin[0]; ix <= imax[0]; ix++)
		  {
			index = d_index3d(ix,iy,iz,top_gmax);
			
			getRectangleCenter(index, coords);
			comp = cell_center[index].comp;

			fprintf(fp, "%f %f %f ", coords[0], coords[1], coords[2]);
			fprintf(fp, "%d ", comp);
		
			if(!gas_comp(comp))
			{
			    fprintf(fp, "%12.5e %12.5e  %12.5e %12.5e %12.5e  %12.5e %12.5e %12.5e ", 
					0.0,0.0, 0.0,0.0,0.0,  0.0,0.0,0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ", 
			    		0.0, 0.0, 0.0, 0.0, 0.0);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ", 
			    		0.0, 0.0, 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e  %12.5e %12.5e %12.5e ", 
					dens[index], pres[index], 
					vel[0][index],
					vel[1][index],
					vel[2][index],
					gnor[0][index],
					gnor[1][index],
					gnor[2][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ",
			    		Gdens[0][index], Gpres[0][index],
					Gvel[0][0][index], Gvel[0][1][index], 
					Gvel[0][2][index]);
			    fprintf(fp, "%12.5e %12.5e   %12.5e %12.5e %12.5e ",
			    		Gdens[1][index], Gpres[1][index],
					Gvel[1][0][index], Gvel[1][1][index], 
					Gvel[1][2][index]);
			}
			fprintf(fp, "\n");
		  }
	}

	fclose(fp);

}

EXPORT  void    tecplot_surface_states(
	const char	*bname,
	FILE		*file,
	SURFACE		*s)
{
	TRI	*tri;
	POINT	*p;
	int	i,npts,ntri;
	Locstate  sl, sr;

	if (bname != NULL)//direct call
	{
	    if ((file = fopen(bname,"w")) == NULL)
	    {
		screen("WARNING in tecplot_surface(), "
		       "can't open %s\n",bname);
		return;
	    }
	    (void) fprintf(file,"TITLE = \"tecplot surface\"\n"
	 		"VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
			"\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\" \n");
	}
	
	//called from tecplot_interface
	if (file == NULL)
	{
	    screen("ERROR, in tecplot_surface, file is NULL\n");
	    clean_up(ERROR);
	}
	if (!(first_tri(s)))
	{
	    screen("WARNING, first bond of the curve is NULL\n");
	    return;
	}

	//count number of points(npts) and number of tris(ntri)
	for (tri = first_tri(s), ntri = 0; !at_end_of_tri_list(tri,s); 
			tri = tri->next, ntri++)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri = first_tri(s), npts = 0; !at_end_of_tri_list(tri,s); 
			tri = tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		}
	    }
	}
	//end counting
	
	fprint_wave_type(file, "ZONE T=\"", wave_type(s), "\"", s->interface);
    	fprintf(file, " N=%d E=%d\nF=FEPOINT, ET=TRIANGLE\n",npts,ntri);

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri = first_tri(s), npts = 0; !at_end_of_tri_list(tri,s); 
			tri = tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		    
		    FT_GetStatesAtPoint(p,Hyper_surf_element(tri),
				Hyper_surf(s),&sl,&sr);
	            fprintf(file,"%15.8e %15.8e %15.8e  %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    	 Coords(p)[1],Coords(p)[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}
	    }
	}
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		fprintf(file,"%d ",Index_of_point(Point_of_tri(tri)[i]));
	    }
	    fprintf(file,"\n");
	}

	if (ntri != s->num_tri)
	{
	    printf("WARNING, num of tri in surface is wrong\n"); 
	}
	if (bname != NULL)
	    fclose(file);

}	/* end tecplot_surface_states */


EXPORT  void    tecplot_interface_states(const char*, INTERFACE	*);

EXPORT  void    tecplot_interface_states(
	const char	*bname,
	INTERFACE	*intfc)
{
	SURFACE	**s;
	char    bname1[200];
	FILE	*file;

	sprintf(bname1, "%s.plt", bname);
	if ((file = fopen(bname1,"w")) == NULL)
	{
	    screen("WARNING in tecplot_interface_states(), "
	           "can't open %s\n",bname1);
	    return;
	}
	(void) fprintf(file,"TITLE = \"tecplot interface\"\n"
	    "VARIABLES = \"x\", \"y\", \"z\", \"PL\", \"PR\", \"DL\", \"DR\" "
	    "\"u\", \"v\", \"w\", \"u1\", \"v1\", \"w1\" \n");

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    tecplot_surface_states(NULL,file,*s);
	}
	fclose(file);
}	/* end tecplot_interface */

boolean G_CARTESIAN::get_ave_state(
	const SWEEP& m_vst,
	int		*ic,
	int		***norset,
	int		comp,
	int		ind)
{
	int		i, j, k, l, num, ic1[3], dir;
	float		gd, gp, gvel[3];
	boolean		found = NO;
	double		**momn = m_vst.momn;
	double		*dens = m_vst.dens;
	double		*pres = m_vst.pres;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	int		index, index0;
	int		icoords[MAXD];
	int 		istart = std::max(0, ic[0] - 3);
        int 		jstart = std::max(0, ic[1] - 3);
        int 		kstart = std::max(0, ic[2] - 3);
        int 		iend = std::min(top_gmax[0],ic[0]+3);
        int 		jend = std::min(top_gmax[1],ic[1]+3);
        int 		kend = std::min(top_gmax[2],ic[2]+3);

	for (i = istart; i <= iend; i++)
        for (j = jstart; j <= jend; j++)
        for (k = kstart; k <= kend; k++)
        {
             if (norset[i][j][k] == 1)
             {
                 found = YES;
                 break;
             }
        }

	if(!found)
	    return NO;

	index0 = d_index(ic,top_gmax,dim);

	num = 0;
	gd = 0.0;
	gp = 0.0;
	gvel[0] = 0.0;
	gvel[1] = 0.0;
	gvel[2] = 0.0;

	for (i = istart; i <= iend; i++)
        for (j = jstart; j <= jend; j++)
        for (k = kstart; k <= kend; k++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    if(i < 0 || i > top_gmax[0] || 
	       j < 0 || j > top_gmax[1] ||
	       k < 0 || k > top_gmax[2])
		continue;
	    if(norset[i][j][k] != 1)
		continue;

	    index = d_index(icoords,top_gmax,dim);

	    if(cell_center[index].comp == comp)
	    {
		gd += dens[index];
		gp += pres[index];
		for (l = 0; l < dim; ++l)
		    gvel[l] += momn[l][index]/dens[index];
	    }
	    else
	    {
		gd += Gdens[ind][index];
		gp += Gpres[ind][index];
		for (l = 0; l < dim; ++l)
		    gvel[l] += Gvel[ind][l][index];
	    }
	    num++;
	}

	Gdens[ind][index0] = gd/num;
	Gpres[ind][index0] = gp/num;
	for (l = 0; l < dim; ++l)
	    Gvel[ind][l][index0] = gvel[l]/num;

	return YES;
}

void G_CARTESIAN::get_ghost_state(
	const SWEEP& m_vst,
	int		comp,
	int		ind)
{
	int			ic[3],index;
	int			c, num;
	boolean			found;
	double			**momn = m_vst.momn;
	double			*dens = m_vst.dens;
	double			*pres = m_vst.pres;
	double			***Gvel = eqn_params->Gvel;
	double			**Gdens = eqn_params->Gdens;
	double			**Gpres = eqn_params->Gpres;
	static	int		***norset;
	static 	int 		loop_count = 0;
	std::list<ToFill> resetThese;
	std::list<ToFill> fillThese;


	if (norset == NULL)
	{
	    int ft_vec_size = 0;
	    for (int i = 0; i < dim; ++i)
	    {
		if (top_gmax[i]+8 > ft_vec_size)
		    ft_vec_size = top_gmax[i]+8;
	    }
	    if(dim == 1)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   1,1,INT);

	    if(dim == 2)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   ft_vec_size,1,INT);
	    if(dim == 3)
	    	FT_TriArrayMemoryAlloc((POINTER*)&norset,ft_vec_size,
				   ft_vec_size,ft_vec_size,INT);
	}

	ToFill aghst;
	
	for (int i = 0; i <= top_gmax[0]; i++)
	for (int j = 0; j <= top_gmax[1]; j++)
	for (int k = 0; k <= top_gmax[2]; k++)
	{
        int icoords[3] = {i,j,k};
        index = d_index(icoords,top_gmax,dim);
	    c = cell_center[index].comp;

	    // for each cell that has component "comp" we
	    // set G* values and mark norset  for that cell to 1
	    if(c == comp)
	    {
            norset[i][j][k] = 1;
            Gdens[ind][index] = dens[index];
            Gpres[ind][index] = pres[index];
            for (int l = 0; l < dim; ++l)
                Gvel[ind][l][index] = momn[l][index]/dens[index];
	    }
	    else
	    {
            aghst.icoords[0] = i;
            aghst.icoords[1] = j;
            aghst.icoords[2] = k;
            //second arg = 1 corresponds to stencil size of 4.... is hardcoded right now...
            if (withinStencilLen(aghst.icoords,1))
            {
                fillThese.push_back(aghst);
            }
            norset[i][j][k] = 0;
	    }
	}

	found = YES;
	num = 1;
	while(found && (num > 0))
	{
	    std::list<ToFill>::iterator it;

	    found = NO;
	    loop_count++;
	    num = 0;

	    resetThese.clear();	
	    for (it=fillThese.begin() ; it != fillThese.end();)
	    {
            found = YES;
            ic[0] = it->icoords[0]; 
            ic[1] = it->icoords[1]; 
            ic[2] = it->icoords[2]; 

            // if no neighbors are 1, return 0.
            if (get_ave_state(m_vst,ic,norset,comp,ind))
            {
                num++;
                norset[ ic[0] ][ ic[1] ][ ic[2] ] = 2;
                aghst.icoords[0] = ic[0];
                aghst.icoords[1] = ic[1];
                aghst.icoords[2] = ic[2]; 
                resetThese.push_back(aghst);
                // erase returns the next valid entery 
                // after the one we just erased.
                it=fillThese.erase(it);
            }
            else
            {
                ++it;
            }
	    }

	    for (it=resetThese.begin(); it != resetThese.end(); it++)
            norset[it->icoords[0]][it->icoords[1]][it->icoords[2]] = 1;
	}

	fillThese.clear();
	resetThese.clear();	
	loop_count = 0;
}

void G_CARTESIAN::GFMGhostState(
	int	*ic,
	int	comp,
	STATE	*ghost_st)
{
	int		i, index;
	double		ncor;
	double		***Gvel = eqn_params->Gvel;
	double		**Gdens = eqn_params->Gdens;
	double		**Gpres = eqn_params->Gpres;
	double		**Gnor = eqn_params->gnor;
	EOS_PARAMS	*eos = eqn_params->eos;

	index = d_index(ic,top_gmax,dim);
	ghost_st->eos = &(eos[comp]);
	ghost_st->dim = dim;

	ncor = 0.0;
	for(i=0; i<dim; i++)
	    ncor += (Gvel[1][i][index] - Gvel[0][i][index])*Gnor[i][index];
		    
	if(comp == 2)
	{
	    ghost_st->pres = Gpres[1][index];
	    ghost_st->dens = Gdens[0][index];
	    for(i=0; i<dim; i++)
		ghost_st->vel[i] = Gvel[0][i][index] + ncor*Gnor[i][index];
	}
	else
	{
	    ghost_st->pres = Gpres[0][index];
	    ghost_st->dens = Gdens[1][index];
	    for(i=0; i<dim; i++)
		ghost_st->vel[i] = Gvel[1][i][index] - ncor*Gnor[i][index];
	}
	for(i=0; i<dim; i++)
	    ghost_st->momn[i] = ghost_st->dens*ghost_st->vel[i];
	
	ghost_st->engy = EosEnergy(ghost_st);
}

void G_CARTESIAN::setNeumannStates(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int 		i,j,index;
	int             ind2[2][2] = {{0,1},{1,0}};
        int             ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic[MAXD];
	double		*vel_intfc = state->vel;
	double		coords[MAXD],coords_ref[MAXD],crx_coords[MAXD];
	double		nor[MAXD],vn,v[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION  dir;
	STATE		st_tmp;

	//st_tmp.eos = state->eos;
	st_tmp.eos = &eqn_params->eos[comp];
	st_tmp.dim = dim;
	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic[i] = icoords[i];
	}
	dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	if (debugging("neumann_buffer"))
	{
	    (void) printf("Entering setNeumannStates()\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	    (void) print_general_vector("vel_intfc = ",vel_intfc,dim,"\n");
	}

	for (i = istart; i <= nrad; ++i)
	{
	    /* Find ghost point */
	    ic[idir] = (nb == 0) ? icoords[idir] - (i - istart + 1) :
                                icoords[idir] + (i - istart + 1);
	
        for (j = 0; j < dim; ++j)
            coords_ref[j] = top_L[j] + ic[j]*top_h[j];//coords_ghost[j]

	    /* Reflect ghost point through intfc-mirror at crossing */
	    vn = 0.0;
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ref[idir];

	    for (j = 0; j < dim; ++j)
	    {
            v[j] = coords_ref[j] - crx_coords[j];
            vn += v[j]*nor[j];
	    }

        //reflect v across the line containing the normal vector
	    for (j = 0; j < dim; ++j)
		    v[j] = 2.0*vn*nor[j] - v[j];

        //desired reflected point
	    for (j = 0; j < dim; ++j)
		    coords_ref[j] = crx_coords[j] + v[j];
			
	    /* Interpolate the state at the reflected point */
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		    m_vst->dens,getStateDens,&st_tmp.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
		    m_vst->pres,getStatePres,&st_tmp.pres,&m_vst->pres[index]);

	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[0],getStateXmom,&st_tmp.momn[0],
			&m_vst->momn[0][index]);
	    if (dim > 1)
		FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[1],getStateYmom,&st_tmp.momn[1],
			&m_vst->momn[1][index]);
	    if (dim > 2)
		FT_IntrpStateVarAtCoords(front,comp,coords_ref,
			m_vst->momn[2],getStateZmom,&st_tmp.momn[2],
			&m_vst->momn[2][index]);
		
        /* Galileo Transformation */
	    vn = 0.0;
        double vel_reflect[3];
	    for (j = 0; j < dim; j++)
	    {
            vel_reflect[j] = st_tmp.momn[j]/st_tmp.dens;
            vn += (vel_reflect[j] - vel_intfc[j])*nor[j];
	    }
            
        /* Only normal component is reflected, 
            relative tangent velocity is zero */
        for (j = 0; j < dim; j++)
	    {
            v[j] = vel_intfc[j] - vn*nor[j];
		    st_tmp.momn[j] = v[j]*st_tmp.dens;
	    }

	    st_tmp.engy = EosEnergy(&st_tmp);

	    /* debugging printout */
	    if (st_tmp.engy < 0.0 || st_tmp.eos->gamma < 0.001)
	    {
		printf("negative engrgy! \n");
		printf("icoords = %d %d %d \n", icoords[0],icoords[1],
						icoords[2]);
		printf("%f %f %f %f %f %f \n",st_tmp.dens,st_tmp.momn[0],
			st_tmp.momn[1],st_tmp.momn[2],st_tmp.pres,
			st_tmp.engy);
		printf("st_tmp.dim = %d, idir = %d, nb = %d \n",
			st_tmp.dim,idir,nb);
		printf("gamma = %f, einf = %f, pinf = %f \n",st_tmp.eos->gamma,
			st_tmp.eos->einf,st_tmp.eos->pinf);
		printf("coords_ref = %f %f %f \n",coords_ref[0],coords_ref[1],
						coords_ref[2]);
		clean_up(0);
	    }

	    if (nb == 0)
	    {
		vst->dens[nrad-i] = st_tmp.dens;
		vst->engy[nrad-i] = st_tmp.engy;
		vst->pres[nrad-i] = st_tmp.pres;
	    	for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-i] = 0.0;
		if (dim == 1)
		   vst->momn[0][nrad-i] = st_tmp.momn[0];
	    	else if (dim == 2)
		    for (j = 0; j < 2; j++)
		    	vst->momn[j][nrad-i] = 
				st_tmp.momn[ind2[idir][j]];
	    	else if (dim == 3)
		    for (j = 0; j < 3; j++)
		    	vst->momn[j][nrad-i] = 
				st_tmp.momn[ind3[idir][j]];
	    }
	    else
	    {
		vst->dens[n+nrad+i-1] = st_tmp.dens;
		vst->engy[n+nrad+i-1] = st_tmp.engy;
		vst->pres[n+nrad+i-1] = st_tmp.pres;
	    	for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+i-1] = 0.0;
		if (dim == 1)
		   vst->momn[0][n+nrad+i-1] = st_tmp.momn[0];
	    	else if (dim == 2)
		    for (j = 0; j < 2; j++)
		    	vst->momn[j][n+nrad+i-1] = 
				st_tmp.momn[ind2[idir][j]];
	    	else if (dim == 3)
		    for (j = 0; j < 3; j++)
		    	vst->momn[j][n+nrad+i-1] = 
				st_tmp.momn[ind3[idir][j]];
	    }
	}
	if (debugging("neumann_buffer"))
	    (void) printf("Leaving setNeumannStates()\n");
}	/* end setNeumannStates */

void G_CARTESIAN::setElasticStates(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
    if (eqn_params->poro_scheme == PORO_SCHEME::RIEMANN)
        setElasticStatesRiem(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
    else if (eqn_params->poro_scheme == PORO_SCHEME::REFLECTION)
        setElasticStatesRFB(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
    else
        setElasticStatesRFB_normal(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
}

//Reflection Boundary Formulation of Porosity -- Allow relative tangential velocity
void G_CARTESIAN::setElasticStatesRFB(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int i,j,index,index_ghost;
	int ind2[2][2] = {{0,1},{1,0}};
    int ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int ic_ghost[MAXD];

	double	coords[MAXD],coords_ref[MAXD],coords_ghost[MAXD],crx_coords[MAXD];
	double	nor[MAXD],v[MAXD],v_ghost[MAXD],v_real[MAXD];
	
	double* vel_intfc = state->vel;
	double poro = eqn_params->porosity;
	
	GRID_DIRECTION  dir;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};

	STATE st_tmp_real;
	STATE st_tmp_ghost;	

	st_tmp_real.dim = dim;
	st_tmp_real.eos = &eqn_params->eos[comp];

    st_tmp_ghost.dim = dim;
    st_tmp_ghost.eos = &eqn_params->eos[comp];
	//st_tmp_ghost.eos = state->eos;

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic_ghost[i] = icoords[i];
	}
	
    dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	if (debugging("elastic_buffer"))
	{
	    (void) printf("\nEntered setElasticStatesRFB():\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	    (void) print_general_vector("vel_intfc = ",vel_intfc,dim,"\n");
	}

	    //if nb = 0, the point is above the boundary, and we
        //           select three points below the boundary
	    //if nb = 1, the point is below the boundary, and we
        //           select three points above the boundary

	for (i = istart; i <= nrad; ++i)
	{
	    //ghost point icoords and index
	    ic_ghost[idir] = (nb == 0) ?
            icoords[idir] - (i - istart + 1) : icoords[idir] + (i - istart + 1);

	    index_ghost = d_index(ic_ghost,top_gmax,dim);

        //ghost point coords
	    for (j = 0; j < dim; ++j)
        {
            coords_ghost[j] = top_L[j] + ic_ghost[j]*top_h[j];
            coords_ref[j] = coords_ghost[j];
	    }
        
        /* Reflect ghost point through intfc-mirror at crossing */
        //first reflect across the grid line containing the intfc crossing 
	    double vn = 0.0;
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];

	    for (j = 0; j < dim; ++j)
	    {
		    v[j] = coords_ref[j] - crx_coords[j];
		    vn += v[j]*nor[j];
	    }
           
        //reflect v across the line containing the normal vector
	    for (j = 0; j < dim; ++j)
		    v[j] = 2.0*vn*nor[j] - v[j];
	    
        //desired reflected point
        for (j = 0; j < dim; ++j)
		    coords_ref[j] = crx_coords[j] + v[j];
			
        /* Interpolate the state at the reflected point */
	    
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->dens,getStateDens,&st_tmp_ghost.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->pres,getStatePres,&st_tmp_ghost.pres,&m_vst->pres[index]);
	    
        for (j = 0; j < dim; ++j)
        {
            FT_IntrpStateVarAtCoords(front,comp,coords_ref,m_vst->momn[j],
                    getStateMom[j],&st_tmp_ghost.momn[j],&m_vst->momn[j][index]);
        }
        
        vn = 0.0;
        double v_reflect[3], v_rel[3];
        for (j = 0; j < dim; j++)
	    {
            v_reflect[j] = st_tmp_ghost.momn[j]/st_tmp_ghost.dens;
            vn += (v_reflect[j] - vel_intfc[j])*nor[j];
	    }
	    
        //Ghost vel has relative normal velocity component equal in magnitude to
        //reflected point's relative normal velocity and going in the opposite direction.
        //TODO: We leave the question of relative tangential velocity wrt to the intfc
        //      to be determined
        //
        //      NEEDS TO BE TESTED
        for (j = 0; j < dim; j++)
        {
            //allow relative tangential velocity
            v_ghost[j] = v_reflect[j] - 2.0*vn*nor[j];
            
            /*
            //zero relative tangential velocity
            v_ghost[j] = vel_intfc[j] - vn*nor[j];
            */
        }

	    st_tmp_real.dens = m_vst->dens[index_ghost];
	    st_tmp_real.pres = m_vst->pres[index_ghost];
	    
        st_tmp_ghost.dens = (1.0 - poro)*st_tmp_ghost.dens + poro*st_tmp_real.dens;
	    st_tmp_ghost.pres = (1.0 - poro)*st_tmp_ghost.pres + poro*st_tmp_real.pres;
	  
        for (j = 0; j < dim; ++j)
        {
            st_tmp_real.momn[j] = m_vst->momn[j][index_ghost];
            v_real[j] = st_tmp_real.momn[j]/st_tmp_real.dens;
            v_ghost[j] = (1.0 - poro)*v_ghost[j] + poro*v_real[j];
            st_tmp_ghost.momn[j] = v_ghost[j]*st_tmp_ghost.dens;
        }
	    
	    st_tmp_ghost.engy = EosEnergy(&st_tmp_ghost);

	    /* debugging printout */
	    if (st_tmp_ghost.engy < 0.0 || st_tmp_ghost.eos->gamma < 0.001)
	    {
            printf("negative engrgy! \n");
            printf("icoords = %d %d %d \n",icoords[0],icoords[1],icoords[2]);
            printf("%f %f %f %f %f %f \n",st_tmp_ghost.dens,st_tmp_ghost.momn[0],
                st_tmp_ghost.momn[1],st_tmp_ghost.momn[2],st_tmp_ghost.pres,
                st_tmp_ghost.engy);
            printf("st_tmp_ghost.dim = %d, idir = %d, nb = %d \n",
                st_tmp_ghost.dim,idir,nb);
            printf("gamma = %f, einf = %f, pinf = %f \n",st_tmp_ghost.eos->gamma,
                st_tmp_ghost.eos->einf,st_tmp_ghost.eos->pinf);
            printf("coords_ref = %f %f %f \n",coords_ref[0],coords_ref[1],
                            coords_ref[2]);
            clean_up(EXIT_FAILURE);
	    }

	    if (nb == 0)
	    {
            vst->dens[nrad-i] = st_tmp_ghost.dens;
            vst->engy[nrad-i] = st_tmp_ghost.engy;
            vst->pres[nrad-i] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][nrad-i] = 0.0;

            if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][nrad-i] = st_tmp_ghost.momn[0];
            }
	    }
	    else
	    {
            /* Debug selectively!
            if (debugging("crx_reflection"))
            {
                    sprintf(fname,"intfc-%d-%d",count,i);
                    sprintf(fname,"intfc-xx");
                    xgraph_2d_reflection(fname,front->grid_intfc,coords,
                    crx_coords,coords_ref,nor);
            }
            */
            vst->dens[n+nrad+i-1] = st_tmp_ghost.dens;
            vst->engy[n+nrad+i-1] = st_tmp_ghost.engy;
            vst->pres[n+nrad+i-1] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][n+nrad+i-1] = 0.0;
    
	    	if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][n+nrad+i-1] = st_tmp_ghost.momn[0];
            }
	    }
	}

	if (debugging("elastic_buffer"))
        (void) printf("Leaving setElasticStatesRFB()\n");
}	/* end setElasticStatesRFB */

//Reflection Boundary Formulation of Porosity -- No relative tangential velocity
void G_CARTESIAN::setElasticStatesRFB_normal(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int i,j,index,index_ghost;
	int ind2[2][2] = {{0,1},{1,0}};
    int ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int ic_ghost[MAXD];

	double	coords[MAXD],coords_ref[MAXD],coords_ghost[MAXD],crx_coords[MAXD];
	double	nor[MAXD],v[MAXD],v_ghost[MAXD],v_real[MAXD];
	
	double* vel_intfc = state->vel;
	double poro = eqn_params->porosity;
	
	GRID_DIRECTION  dir;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};

	STATE st_tmp_real;
	STATE st_tmp_ghost;	

	st_tmp_real.dim = dim;
	st_tmp_real.eos = &eqn_params->eos[comp];

    st_tmp_ghost.dim = dim;
    st_tmp_ghost.eos = &eqn_params->eos[comp];
	//st_tmp_ghost.eos = state->eos;

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic_ghost[i] = icoords[i];
	}
	
    dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	if (debugging("elastic_buffer"))
	{
	    (void) printf("\nEntered setElasticStatesRFB_normal():\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	    (void) print_general_vector("vel_intfc = ",vel_intfc,dim,"\n");
	}

	    //if nb = 0, the point is above the boundary, and we
        //           select three points below the boundary
	    //if nb = 1, the point is below the boundary, and we
        //           select three points above the boundary

	for (i = istart; i <= nrad; ++i)
	{
	    //ghost point icoords and index
	    ic_ghost[idir] = (nb == 0) ?
            icoords[idir] - (i - istart + 1) : icoords[idir] + (i - istart + 1);

	    index_ghost = d_index(ic_ghost,top_gmax,dim);

        //ghost point coords
	    for (j = 0; j < dim; ++j)
        {
            coords_ghost[j] = top_L[j] + ic_ghost[j]*top_h[j];
            coords_ref[j] = coords_ghost[j];
	    }
        
        /* Reflect ghost point through intfc-mirror at crossing */
        //first reflect across the grid line containing the intfc crossing 
	    double vn = 0.0;
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];

	    for (j = 0; j < dim; ++j)
	    {
		    v[j] = coords_ref[j] - crx_coords[j];
		    vn += v[j]*nor[j];
	    }
           
        //reflect v across the line containing the normal vector
	    for (j = 0; j < dim; ++j)
		    v[j] = 2.0*vn*nor[j] - v[j];
	    
        //desired reflected point
        for (j = 0; j < dim; ++j)
		    coords_ref[j] = crx_coords[j] + v[j];
			
        /* Interpolate the state at the reflected point */
	    
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->dens,getStateDens,&st_tmp_ghost.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->pres,getStatePres,&st_tmp_ghost.pres,&m_vst->pres[index]);
	    
        for (j = 0; j < dim; ++j)
        {
            FT_IntrpStateVarAtCoords(front,comp,coords_ref,m_vst->momn[j],
                    getStateMom[j],&st_tmp_ghost.momn[j],&m_vst->momn[j][index]);
        }
        
        //Compute relative normal velocity in frame of interface crossing.
        vn = 0.0;
        double v_reflect[3], v_rel[3];
        for (j = 0; j < dim; j++)
	    {
            v_reflect[j] = st_tmp_ghost.momn[j]/st_tmp_ghost.dens;
            vn += (v_reflect[j] - vel_intfc[j])*nor[j];
	    }
	    
        //Ghost vel has relative normal velocity component equal in magnitude to
        //reflected point's relative normal velocity and going in the opposite direction.
        //TODO: We leave the question of relative tangential velocity wrt to the intfc
        //      to be determined
        //
        //      NEEDS TO BE TESTED
        for (j = 0; j < dim; j++)
        {
            /*
            //allow relative tangential velocity
            v_ghost[j] = v_reflect[j] - 2.0*vn*nor[j];
            */

            //zero relative tangential velocity
            v_ghost[j] = vel_intfc[j] - vn*nor[j];
        }

	    st_tmp_real.dens = m_vst->dens[index_ghost];
	    st_tmp_real.pres = m_vst->pres[index_ghost];
	    
        st_tmp_ghost.dens = (1.0 - poro)*st_tmp_ghost.dens + poro*st_tmp_real.dens;
	    st_tmp_ghost.pres = (1.0 - poro)*st_tmp_ghost.pres + poro*st_tmp_real.pres;
	  
        for (j = 0; j < dim; ++j)
        {
            st_tmp_real.momn[j] = m_vst->momn[j][index_ghost];
            v_real[j] = st_tmp_real.momn[j]/st_tmp_real.dens;
            v_ghost[j] = (1.0 - poro)*v_ghost[j] + poro*v_real[j];
            st_tmp_ghost.momn[j] = v_ghost[j]*st_tmp_ghost.dens;
        }
	    
	    st_tmp_ghost.engy = EosEnergy(&st_tmp_ghost);

	    /* debugging printout */
	    if (st_tmp_ghost.engy < 0.0 || st_tmp_ghost.eos->gamma < 0.001)
	    {
            printf("negative engrgy! \n");
            printf("icoords = %d %d %d \n",icoords[0],icoords[1],icoords[2]);
            printf("%f %f %f %f %f %f \n",st_tmp_ghost.dens,st_tmp_ghost.momn[0],
                st_tmp_ghost.momn[1],st_tmp_ghost.momn[2],st_tmp_ghost.pres,
                st_tmp_ghost.engy);
            printf("st_tmp_ghost.dim = %d, idir = %d, nb = %d \n",
                st_tmp_ghost.dim,idir,nb);
            printf("gamma = %f, einf = %f, pinf = %f \n",st_tmp_ghost.eos->gamma,
                st_tmp_ghost.eos->einf,st_tmp_ghost.eos->pinf);
            printf("coords_ref = %f %f %f \n",coords_ref[0],coords_ref[1],
                            coords_ref[2]);
            clean_up(EXIT_FAILURE);
	    }

	    if (nb == 0)
	    {
            vst->dens[nrad-i] = st_tmp_ghost.dens;
            vst->engy[nrad-i] = st_tmp_ghost.engy;
            vst->pres[nrad-i] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][nrad-i] = 0.0;

            if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][nrad-i] = st_tmp_ghost.momn[0];
            }
	    }
	    else
	    {
            /* Debug selectively!
            if (debugging("crx_reflection"))
            {
                    sprintf(fname,"intfc-%d-%d",count,i);
                    sprintf(fname,"intfc-xx");
                    xgraph_2d_reflection(fname,front->grid_intfc,coords,
                    crx_coords,coords_ref,nor);
            }
            */
            vst->dens[n+nrad+i-1] = st_tmp_ghost.dens;
            vst->engy[n+nrad+i-1] = st_tmp_ghost.engy;
            vst->pres[n+nrad+i-1] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][n+nrad+i-1] = 0.0;
    
	    	if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][n+nrad+i-1] = st_tmp_ghost.momn[0];
            }
	    }
	}

	if (debugging("elastic_buffer"))
        (void) printf("Leaving setElasticStatesRFB_normal()\n");
}	/* end setElasticStatesRFB_normal */

//Riemann Problem Formulation of Porosity
void G_CARTESIAN::setElasticStatesRiem(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int i,j,index,index_ghost;
	int ic_ghost[MAXD];

	double	coords[MAXD],coords_ref[MAXD],coords_ghost[MAXD],crx_coords[MAXD];
	double	nor[MAXD],v[MAXD],v_ghost[MAXD],v_real[MAXD];
	
    double vn, vn_intfc;
	
	GRID_DIRECTION  dir;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};

	STATE sl, sr, state_ghost;

	sl.dim = sr.dim = state_ghost.dim = dim;
	sl.eos = sr.eos = state_ghost.eos = &eqn_params->eos[comp];

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic_ghost[i] = icoords[i];
	}
    dir = (nb == 0) ? ldir[idir] : rdir[idir];

    /*
	if (debugging("elastic_buffer"))
	{
	    (void) printf("\nEntered setElasticStatesRiem():\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	}
    */

	for (i = istart; i <= nrad; ++i)
	{
	    //ghost point icoords and index
	    ic_ghost[idir] = (nb == 0) ?
            icoords[idir] - (i - istart + 1) : icoords[idir] + (i - istart + 1);

        //ghost point coords
	    for (j = 0; j < dim; ++j)
        {
            coords_ghost[j] = top_L[j] + ic_ghost[j]*top_h[j];
	    }
        
        index_ghost = d_index(ic_ghost,top_gmax,dim);
	    COMPONENT comp_ghost = cell_center[index_ghost].comp;

        //Find the closest interface point
        boolean status;
        double intrp_a[3];
        HYPER_SURF_ELEMENT* nearHse;
        HYPER_SURF* nearHs;

        status = FT_FindNearestIntfcPointInRange(front,comp_ghost,
                coords_ghost,INCLUDE_BOUNDARIES,crx_coords,intrp_a,
                &nearHse,&nearHs,5);
        
        if (!status) 
        {
            LOC();
            printf("ERROR: could not find interface point\n");
            clean_up(EXIT_FAILURE);
        }
        
        //Get 2 points straddling interface in normal direction
        double pl[MAXD], pr[MAXD], nor[MAXD];

        //TODO: Use all three points of triangle in approximation,
        //      and use the triangle normal like in compute_total_canopy_force3d()
        TRI* nearTri = Tri_of_hse(nearHse);
        FT_NormalAtPoint(Point_of_tri(nearTri)[0],front,nor,comp);
        STATE* state_intfc = (STATE*)left_state(Point_of_tri(nearTri)[0]);
            //nor = Tri_normal_vector(nearTri);//USE THIS
        double h = FT_GridSizeInDir(nor,front);

        for (j = 0; j < 3; ++j)
        {
            pl[j] = crx_coords[j] + 1.5*h*nor[j];
            pr[j] = crx_coords[j] - 1.5*h*nor[j];//on ghost point side
        }

        //Interpolate states for the 2 points
        FT_IntrpStateVarAtCoords(front,comp,pl,
                m_vst->dens,getStateDens,&sl.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,pl,
                m_vst->pres,getStatePres,&sl.pres,&m_vst->pres[index]);
	    
        FT_IntrpStateVarAtCoords(front,comp_ghost,pr,
                m_vst->dens,getStateDens,&sr.dens,&m_vst->dens[index_ghost]);
	    FT_IntrpStateVarAtCoords(front,comp_ghost,pr,
                m_vst->pres,getStatePres,&sr.pres,&m_vst->pres[index_ghost]);
        
        //TODO: retain the tangential velocities of the interpolated states
        //      and add to the normal ghost velocity obtained from the
        //      riemann problem.
        //
        //Using relative velocity wrt to interface velocity
        //TODO: can we do better than this interpolation?
        /*
        double intfc_dens;
        double intfc_momn[3], vel_intfc[3];
        FT_IntrpStateVarAtCoords(front,comp,crx_coords,
                m_vst->dens,getStateDens,&intfc_dens,&m_vst->dens[index]);
        */

        //TODO: Should the left state velocity be set to the right state velocity
        //      as in GFM for strong shock interaction part 1&2 (or is that a typo)
        //      Don't think it is a typo as results are much better.
        //      Need to understand why.
        double vl[3], vr[3];
        for (j = 0; j < dim; ++j)
        {
            /*
            FT_IntrpStateVarAtCoords(front,comp,crx_coords,m_vst->momn[j],
                    getStateMom[j],&intfc_momn[j],&m_vst->momn[j][index]);
            vel_intfc[j] = intfc_momn[j]/intfc_dens;
            */
            
            FT_IntrpStateVarAtCoords(front,comp,pl,m_vst->momn[j],
                    getStateMom[j],&sl.momn[j],&m_vst->momn[j][index]);
            vl[j] = sl.momn[j]/sl.dens - state_intfc->vel[j];//relative left state vel
            //vl[j] = sl.momn[j]/sl.dens - vel_intfc[j];//relative left state vel
                //vl[j] = sl.momn[j]/sl.dens;

            FT_IntrpStateVarAtCoords(front,comp_ghost,pr,m_vst->momn[j],
                    getStateMom[j],&sr.momn[j],&m_vst->momn[j][index_ghost]);
            vr[j] = sr.momn[j]/sr.dens - state_intfc->vel[j];//relative right state vel
            //vr[j] = sr.momn[j]/sr.dens - vel_intfc[j];//relative right state vel
                //vr[j] = sr.momn[j]/sr.dens;
        }

        double nor_vl = 0.0;
        double nor_vr = 0.0;
        for (j = 0; j < 3; ++j)
        {
            //relative normal velocities
            nor_vl += vl[j]*nor[j];
            nor_vr += vr[j]*nor[j];
        }

        /*
        //DEBUG: assuming centerline is through (0.5, 0.5, z)
        double ctrlinedist = sqrt(sqr(coords[0] - 0.5) + sqr(coords[1] - 0.5));
        double tolprint = sqrt(sqr(top_h[0]) + sqr(top_h[1]));
        if (ctrlinedist < tolprint)
        {
            printf("comp = %d\n",comp);
            printf("comp_ghost = %d\n",comp_ghost);
            print_general_vector("coords = ",coords,dim,"\n");
            print_general_vector("coords_ghost = ",coords_ghost,dim,"\n");
            print_general_vector("crx_coords = ",crx_coords,dim,"\n");
            print_general_vector("nor = ",nor,dim,"\n");
            print_general_vector("pr = ",pr,dim,"\n");
            print_general_vector("pl = ",pl,dim,"\n");

            printf("input states: sl sr\n");
            printf("\tdens: %f %f\n",sl.dens,sr.dens);
            printf("\tvn: %f %f\n",nor_vl,nor_vr);
            printf("\tpres: %f %f\n",sl.pres,sr.pres);
        }
        */

        //solve 1d riemann problem in interface normal direction
        RIEMANN_INPUT riem_input;
        RIEMANN_SOLN riem_soln;

        riem_input.left_state.d = sl.dens;
        riem_input.left_state.p = sl.pres;
        riem_input.left_state.u = nor_vr;//Typo or correct? UPDATE: not a typo
            //riem_input.left_state.u = nor_vl;
        riem_input.left_state.gamma = sl.eos->gamma;

        riem_input.right_state.d = sr.dens;
        riem_input.right_state.p = sr.pres;
        riem_input.right_state.u = nor_vr;
        riem_input.right_state.gamma = sr.eos->gamma;

        bool rp_status;
        rp_status = RiemannSolution(riem_input,&riem_soln);
        if (!rp_status)
        {
            printf("\nERROR: RiemannSolution()\n");
            printf("input states: sl sr\n");
            printf("\tpres: %f %f\n",sl.pres,sr.pres);
            printf("\tdens: %f %f\n",sl.dens,sr.dens);
            printf("\tvn: %f %f\n",nor_vl,nor_vr);
            clean_up(EXIT_FAILURE);
        }

        /*
        //DEBUG
        if (ctrlinedist < tolprint)
        {
            printf("\nsoln states, velocity\n");
            printf("right_state.u = %g\n",riem_soln.right_state.u);
            printf("right_center_state.u = %g\n",riem_soln.right_center_state.u);
            printf("left_center_state.u = %g\n",riem_soln.left_center_state.u);
            printf("left_state.u = %g\n",riem_soln.left_state.u);

            printf("\nsoln states, density\n");
            printf("right_state.d = %g\n",riem_soln.right_state.d);
            printf("right_center_state.d = %g\n",riem_soln.right_center_state.d);
            printf("left_center_state.d = %g\n",riem_soln.left_center_state.d);
            printf("left_state.d = %g\n",riem_soln.left_state.d);
           
            printf("\nsoln states, pressure\n");
            printf("right_state.p = %g\n",riem_soln.right_state.p);
            printf("right_center_state.p = %g\n",riem_soln.right_center_state.p);
            printf("left_center_state.p = %g\n",riem_soln.left_center_state.p);
            printf("left_state.p = %g\n",riem_soln.left_state.p);
        }
        */
       
        ///////////////////////////////////////////////////////////////////////
        /*
        ///////////////////////////////////////////////////////////////////////
        //TODO: This was first attempt that used the center state
        //          (centered on interface) soln.
        RIEM_STATE riem_soln_intfc;
        rp_status = RiemannSolnAtXi(&riem_soln,&riem_soln_intfc,0.0);
        if (!rp_status)
        {
            printf("ERROR: RiemannSolnAtXi()\n");
            clean_up(EXIT_FAILURE);
        }

        //printf("riem_soln_intfc:\n");
        //printf("\t(d,u,p) = %f %f %f\n",riem_soln_intfc.d,
          //      riem_soln_intfc.u,riem_soln_intfc.p);
        
        //Assign the midpoint solution state to the ghost state
        double dens_ghost = riem_soln_intfc.d;
        double pres_ghost = riem_soln_intfc.p;
        double vn_ghost = riem_soln_intfc.u;
        */

        /*
        RIEM_STATE left_center_state = riem_soln.left_center_state;
        double dens_ghost = left_center_state.d;
        double pres_ghost = left_center_state.p;
        double vn_ghost = left_center_state.u;
        */
        
        RIEM_STATE right_center_state = riem_soln.right_center_state;
        double dens_ghost = right_center_state.d;
        double pres_ghost = right_center_state.p;
        double vn_ghost = right_center_state.u;
        //TODO: Try using the interpolated left state pressure and density values,
        //      and only derive the velocity from the riemann solution.
            //double dens_ghost = sl.dens; //copy from the interpolated left state
            //double pres_ghost = sl.pres;  //copy from the interpolated left state

        //take weighted average using porosity to get the modified ghost point
	    double poro = eqn_params->porosity;
        state_ghost.dens = (1.0 - poro)*dens_ghost + poro*m_vst->dens[index_ghost];
        state_ghost.pres = (1.0 - poro)*pres_ghost + poro*m_vst->pres[index_ghost];
        
        double v_real[3];
        for (j = 0; j < dim; ++j)
        {
            v_ghost[j] = vn_ghost*nor[j] + state_intfc->vel[j];
                //v_ghost[j] = vn_ghost*nor[j] + vel_intfc[j];

            //reflect normal velocity component and convert back to world frame
                //v_ghost[j] = vel_intfc[j] - vn_ghost*nor[j];
                    //v_ghost[j] = vn_ghost*nor[j];
            v_real[j] = m_vst->momn[j][index_ghost]/m_vst->dens[index_ghost];
            v_ghost[j] = (1.0 - poro)*v_ghost[j] + poro*v_real[j];
            state_ghost.momn[j] = v_ghost[j]*state_ghost.dens;
        }
        state_ghost.engy = EosEnergy(&state_ghost);

	    // debugging printout
	    if (state_ghost.engy < 0.0 || state_ghost.eos->gamma < 0.001)
	    {
            printf("ERROR: Negative Energy! \n");
            printf("icoords = %d %d %d \n",icoords[0],icoords[1],icoords[2]);
            printf("%f %f %f %f %f %f \n",state_ghost.dens,state_ghost.momn[0],
                state_ghost.momn[1],state_ghost.momn[2],state_ghost.pres,
                state_ghost.engy);
            printf("state_ghost.dim = %d, idir = %d, nb = %d \n",
                state_ghost.dim,idir,nb);
            printf("gamma = %f, einf = %f, pinf = %f \n",state_ghost.eos->gamma,
                state_ghost.eos->einf,state_ghost.eos->pinf);
            printf("coords_ghost = %f %f %f \n",coords_ghost[0],coords_ghost[1],
                            coords_ghost[2]);
            clean_up(EXIT_FAILURE);
	    }

	    int ind2[2][2] = {{0,1},{1,0}};
        int ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};

	    if (nb == 0)
	    {
            vst->dens[nrad-i] = state_ghost.dens;
            vst->engy[nrad-i] = state_ghost.engy;
            vst->pres[nrad-i] = state_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][nrad-i] = 0.0;

            if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = state_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][nrad-i] = state_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][nrad-i] = state_ghost.momn[0];
            }
	    }
	    else
	    {
            vst->dens[n+nrad+i-1] = state_ghost.dens;
            vst->engy[n+nrad+i-1] = state_ghost.engy;
            vst->pres[n+nrad+i-1] = state_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][n+nrad+i-1] = 0.0;
    
	    	if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = state_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][n+nrad+i-1] = state_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][n+nrad+i-1] = state_ghost.momn[0];
            }
	    }
	}

	if (debugging("elastic_buffer"))
        (void) printf("Leaving setElasticStatesRiem()\n");
}	/* end setElasticStatesRiem */

void G_CARTESIAN::setDirichletStates(
	STATE		*crx_state,
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	int		*icoords,
	int		dir,
	int		nb,
	int		n,
	int		istart)
{
	int		j, k, index;
	STATE 		*state;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};

	if (nb == 0)
	{
	  if (boundary_state(hs) != NULL)
	  {
	    //preset state bdry
	    state = (STATE*)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		vst->dens[nrad-k] = state->dens;
		vst->engy[nrad-k] = state->engy;
		vst->pres[nrad-k] = state->pres;
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-k] = 0.0;
		if (dim == 1)
		    vst->momn[0][nrad-k] = state->momn[0];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][nrad-k] = state->momn[ind2[dir][j]];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-k] = state->momn[ind3[dir][j]];
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "cF_flowThroughBoundaryState") == 0)
	  {
	    //flow through bdry
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		vst->dens[nrad-k] = m_vst->dens[index];
		vst->engy[nrad-k] = m_vst->engy[index];
		vst->pres[nrad-k] = m_vst->pres[index];
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-k] = 0.0;
		if (dim == 1)
		    vst->momn[0][nrad-k] = m_vst->momn[0][index];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][nrad-k] = m_vst->momn[ind2[dir][j]][index];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][nrad-k] = m_vst->momn[ind3[dir][j]][index];
	    }
	  }
	  else
	  {
	    (void) printf("Unimplemented Dirichlet boundary type!\n");
	    clean_up(ERROR);
	  }
	}
	else
	{
	  if (boundary_state(hs) != NULL)
	  {
	    state = (STATE*)boundary_state(hs);
	    for (k = istart; k <= nrad; ++k)
	    {
		vst->dens[n+nrad+k-1] = state->dens;
		vst->engy[n+nrad+k-1] = state->engy;
		vst->pres[n+nrad+k-1] = state->pres;
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+k-1] = 0.0;
		if (dim == 1)
		    vst->momn[0][n+nrad+k-1] = state->momn[0];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][n+nrad+k-1] = state->momn[ind2[dir][j]];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+k-1] = state->momn[ind3[dir][j]];
	    }
	  }
	  else if (boundary_state_function(hs) &&
              strcmp(boundary_state_function_name(hs),
	      "cF_flowThroughBoundaryState") == 0)
	  {
	    for (k = istart; k <= nrad; ++k)
	    {
		index = d_index(icoords,top_gmax, dim);
		vst->dens[n+nrad+k-1] = m_vst->dens[index];
		vst->engy[n+nrad+k-1] = m_vst->engy[index];
		vst->pres[n+nrad+k-1] = m_vst->pres[index];
		
		for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+k-1] = 0.0;
		if (dim == 1)
		    vst->momn[0][n+nrad+k-1] = m_vst->momn[0][index];
		else if (dim == 2)
		  for (j = 0; j < 2; j++)
		    vst->momn[j][n+nrad+k-1] = 
					m_vst->momn[ind2[dir][j]][index];
		else if (dim == 3)
		  for (j = 0; j < 3; j++)
		    vst->momn[j][n+nrad+k-1] = 
					m_vst->momn[ind3[dir][j]][index];
	    }
	  }
	  else
	  {
	    (void) printf("Unimplemented Dirichlet boundary type!\n");
	    clean_up(ERROR);
	  }
	}
}

void G_CARTESIAN::initSampleVelocity(char *in_name)
{
        FILE *infile;
	static SAMPLE *sample;
	char *sample_type;
	double *sample_line;

	infile = fopen(in_name,"r");
	FT_ScalarMemoryAlloc((POINTER*)&sample,sizeof(SAMPLE));
	sample_type = sample->sample_type;
	sample_line = sample->sample_coords;
	dim = front->rect_grid->dim;

	if (dim == 2)
	{
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf",sample_line);
            (void) printf(" %f\n",sample_line[0]);
	}
	else if (dim == 3)
        {
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf %lf",sample_line,sample_line+1);
            (void) printf(" %f %f\n",sample_line[0],sample_line[1]);
        }
        CursorAfterString(infile,"Enter the start step for sample: ");
        fscanf(infile,"%d",&sample->start_step);
        (void) printf("%d\n",sample->start_step);
        CursorAfterString(infile,"Enter the end step for sample: ");
        fscanf(infile,"%d",&sample->end_step);
        (void) printf("%d\n",sample->end_step);
        CursorAfterString(infile,"Enter the step interval for sample: ");
        fscanf(infile,"%d",&sample->step_interval);
        (void) printf("%d\n",sample->step_interval);
	front->sample = sample;
        fclose(infile);
}	/* end initSampleVelocity */

void G_CARTESIAN::checkCorrectForTolerance(STATE *state)
{
	if (state->dens < min_dens)
	    state->dens = min_dens;
	if (state->pres < min_pres)
	    state->pres = min_pres;
	state->engy = EosEnergy(state);
}	/* end checkCorrectForTolerance */

//TODO: Does this still work as intended with the fabric index coating algorithm?
//      If the components on both side of the fabric are assigned values of 
//      GAS_COMP1 and GAS_COMP2 (= 2 and 3 respectively), then this function
//      will return NO. Is that what we want???
boolean G_CARTESIAN::needBufferFromIntfc(
	COMPONENT domain_comp,
	COMPONENT comp)
{
	if (eqn_params->tracked)
	    return (domain_comp != comp) ? YES : NO;
	else
	    return (gas_comp(comp)) ? NO : YES;
}	/* needBufferFromIntfc */


bool G_CARTESIAN::withinStencilLen(int *icrds, int stencil)
{
        int istart = std::max(0, icrds[0] - stencil);
        int jstart = std::max(0, icrds[1] - stencil);
        int kstart = std::max(0, icrds[2] - stencil);
        int iend = std::min(top_gmax[0],icrds[0]+stencil);
        int jend = std::min(top_gmax[1],icrds[1]+stencil);
        int kend = std::min(top_gmax[2],icrds[2]+stencil);

        int index  =  d_index(icrds,top_gmax,dim);
        int mycomp =  cell_center[index].comp;

        int i,j,k;
        if (dim == 2)
        {
            for (i = istart; i <= iend; i++)
            for (j = jstart; j <= jend; j++)
            {
                int ic[3];
                ic[0] = i; ic[1] = j; ic[2] = 0;
                index  =  d_index(ic,top_gmax,dim);
                if(mycomp != cell_center[index].comp)
                    return true;
            }
            return NO;
        }
        else if (dim == 3)
        {
            for (i = istart; i <= iend; i++)
            for (j = jstart; j <= jend; j++)
            for (k = kstart; k <= kend; k++)
            {
                int ic[3];
                ic[0] = i; ic[1] = j; ic[2] = k;
                index = d_index(ic,top_gmax,dim);
                if(mycomp != cell_center[index].comp)
                    return true;
            }
            return NO;
        }
	return YES;
}

void G_CARTESIAN::addFluxAlongGridLine(
	int idir,
	int *grid_icoords,
	double dt,
	SWEEP *m_vst,
    FSWEEP *m_flux)
{
	int i,l,n,index;
	SCHEME_PARAMS scheme_params;
	EOS_PARAMS	*eos;
	
    /*
    static SWEEP vst;       //SWEEP is a data structure storing state values
	static FSWEEP vflux;    //FSWEEP is a data structure storing flux of state values
	static boolean first = YES;
	*/

    COMPONENT comp;
	int seg_min,seg_max;
	static int icoords[MAXD];
	int icoords_next[MAXD];
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)(front->extra1);
	INTERFACE *grid_intfc=front->grid_intfc;
	double crx_coords[MAXD];
	HYPER_SURF *hs;
	SURFACE **s;
	STATE *state;
	
    GRID_DIRECTION	ldir[3]={WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3]={EAST,NORTH,UPPER};

//	printf("icoords={%d,%d,%d}\n",icoords[0],icoords[1],icoords[2]);
//	printf("icoords=%d,%d,%d, wave_type=%d\n",grid_icoords[0],grid_icoords[1],grid_icoords[2],wave_type(hs));
	
    SWEEP vst;
    FSWEEP vflux;
    allocDirVstFlux(&vst,&vflux);
	/*if (first)
    {
        //TODO: more efficient to free and release?
        //      freeing and reallocating could keep from fragmenting/page faulting
        //      and improve efficiency...
        first = NO;
        allocDirVstFlux(&vst,&vflux);
    }*/

	scheme_params.lambda = dt/top_h[idir];
    scheme_params.beta = 0.0;
	scheme_params.artificial_compression = eqn_params->articomp;

    for (i = 0; i < dim; ++i)
	    icoords[i] = grid_icoords[i];


    //TODO: use these below if checking 2 crossings..
    //      if fabric not seperated well, will be able to tell from 
    //      distance between the crossings
	double ldir_crx_coords[MAXD];
	double rdir_crx_coords[MAXD];
    
    //TODO: See the TODO needBufferFromIntfc() regarding gas_comp().
    //      Need to make sure both work correctly with index coating
    //      algorithm used for ELASTIC_BOUNDARYs
    
    seg_min = imin[idir];	
	while (seg_min <= imax[idir])
	{
	    for (; seg_min <= imax[idir]; ++seg_min)
	    {
		    icoords[idir] = seg_min;
	    	index = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[index];
	    	if (gas_comp(comp)) break;
	    }

	    if (seg_min > imax[idir]) break;
	    
        //TODO: is this good enough zeroing?
        //      what about the +7 in size value in allocDirVstFlux()??
        //      For now should be safe since allocating and freeing everytime.
        for (i = 0; i <= top_gmax[idir]; ++i)
	    {
	    	vst.dens[i] = 0.0; 
	    	vst.pres[i] = 0.0; 
	    	vst.engy[i] = 0.0; 
	    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
	    }
	    
        i = seg_min;
	    icoords[idir] = i;
	    index = d_index(icoords,top_gmax,dim);
	    comp = top_comp[index];
	    n = 0;

        vst.dens[n+nrad] = m_vst->dens[index];
        vst.engy[n+nrad] = m_vst->engy[index];
        vst.pres[n+nrad] = m_vst->pres[index];
	    
        //index cycling so 0-th component aligned along the idir-th gridline,
        //and remaining chosen to form right hand coordinate system
        for (l = 0; l < dim; ++l)
            vst.momn[l][n+nrad] = m_vst->momn[(l+idir)%dim][index];
	    for (l = dim; l < 3; ++l)
            vst.momn[l][n+nrad] = 0.0;
	    
        seg_max = i;
	    n++;

//	    printf("Component=%d\n",comp);

	    for (i = seg_min + 1; i <= imax[idir]; i++)
	    {
            icoords[idir] = i;
            index = d_index(icoords,top_gmax,dim);
	
            for (int ii = 0; ii < dim; ++ii)
                icoords_next[ii] = icoords[ii];
            icoords_next[idir]++;
            
            boolean status1, status2;
            
            //TODO: use ldir_crx_coords and rdir_crx_coords to differentiate crossings
            status1 = FT_StateStructAtGridCrossing(front,grid_intfc,
                    icoords,rdir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);

            //TODO: status2 never gets checked...
            //      Guessing is for case that fabric is folded, but the bugs never
            //      quite got worked out so it was disabled.
            status2 = FT_StateStructAtGridCrossing(front,grid_intfc,
                    icoords_next,ldir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);

            //TODO: why did cxxu abandon/comment out these 2 blocks??
            //      Dead code or unfinished code?

//		For the following part, if needBufferFromIntfc is true, which means
// 		it meets a boundary, we have a break. If needBufferFromIntfc is not true,
// 		we check the two statuses. If one of them is true, we still have a break.
// 		If all three of them are false, we assume that it does not meet a boundary
   
            //if (needBufferFromIntfc(comp,top_comp[index]))
            //{
            //    printf("get boundary \n");
            //    break;
		    //}
//		if (status1){
//		    
//		    printf("Outside::: The wave_type is %d\n", wave_type(hs));
//		    if (!needBufferFromIntfc(comp,top_comp[index])){
//			seg_max=i;
//			n++;
//		    }
//		    break;
//		}
//
            if (needBufferFromIntfc(comp,top_comp[index]))
            {
                printf("get boundary \n");
                break;
            }
		    else
            {
                vst.dens[n+nrad] = m_vst->dens[index];
                vst.engy[n+nrad] = m_vst->engy[index];
                vst.pres[n+nrad] = m_vst->pres[index];

                for (l = 0; l < dim; ++l)
                    vst.momn[l][n+nrad] = m_vst->momn[(l+idir)%dim][index];
                for (l = dim; l < 3; ++l)
                    vst.momn[l][n+nrad] = 0.0;

                n++;
                if (status1)
                {

    //			printf("icoords[0]=%d icoords[1]=%d, icoords[2]=%d, wave_type=%d\n",icoords[0],icoords[1], icoords[2],wave_type(hs));

   //			if (wave_type(hs)==7) //NOTE: 7 is NUEMANN_BOUNDARY
   //  			{
   //                 for (int ii=0;ii<dim;ii++)
   //                 {
   //                    printf("icoords[%d]=%d,",ii,icoords[ii]);
   //                 }
   //                 printf("    found\n");
   //                 
   //                 printf("crx_coords[0]=%f\n",crx_coords[0]);
   //                 for (int ii=0;ii<dim;ii++)
   //                 {
   //                     printf("crx_coords[%d]=%f\n",ii,crx_coords[ii]);
   //                 }

   //             }
   // 
                    
                    seg_max = i++;
                    break;
                 }
            }
    
            seg_max = i;
        }

        //Elastic Boundary and Porosity accounted for in appendGhostBuffer()
	    icoords[idir] = seg_min;
	    appendGhostBuffer(&vst,m_vst,0,icoords,idir,0);
	        //appendGhostBuffer(&vst,m_vst,n,icoords,idir,0);
	    icoords[idir] = seg_max;
	    appendGhostBuffer(&vst,m_vst,n,icoords,idir,1);
	    
	    eos = &(eqn_params->eos[comp]);
	    EosSetTVDParams(&scheme_params, eos);
	    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
	    n = 0;
	    for (i = seg_min; i <= seg_max; ++i)
	    {
            icoords[idir] = i;
            index = d_index(icoords,top_gmax,dim);

            m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
            m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];

            for (l = 0; l < dim; ++l)
            {
                m_flux->momn_flux[(l+idir)%dim][index] +=
                    vflux.momn_flux[l][n+nrad];
            }
            for (l = dim; l < 3; ++l)
                m_flux->momn_flux[l][index] = 0.0;
            
            n++;
	    }

	    seg_min = seg_max + 1;
	}
        
    //TODO: see above todo in if (first) block
    freeDirVstFlux(&vst,&vflux);
}	/* end addFluxAlongGridLine */



/*
void G_CARTESIAN::addFluxAlongGridLine(
	int idir,
	int *grid_icoords,
	double dt,
	SWEEP *m_vst,
        FSWEEP *m_flux)
{
	int i,l,n,index;
	SCHEME_PARAMS scheme_params;
	EOS_PARAMS	*eos;
	SWEEP vst;
	FSWEEP vflux;
	COMPONENT comp;
	int seg_min,seg_max;
	int icoords[MAXD];
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)(front->extra1);
	
        allocDirVstFlux(&vst,&vflux);
	scheme_params.lambda = dt/top_h[idir];
        scheme_params.beta = 0.0;
	scheme_params.artificial_compression = eqn_params->articomp;
	for (i = 0; i < dim; ++i)
	    icoords[i] = grid_icoords[i];
	seg_min = imin[idir];
	while (seg_min <= imax[idir])
	{
	    for (; seg_min <= imax[idir]; ++seg_min)
	    {
		icoords[idir] = seg_min;
	    	index = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[index];
	    	if (gas_comp(comp)) break;
	    }
	    if (seg_min > imax[idir]) break;
	    for (i = 0; i <= top_gmax[idir]; ++i)
	    {
	    	vst.dens[i] = 0.0; 
	    	vst.pres[i] = 0.0; 
	    	vst.engy[i] = 0.0; 
	    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
	    }
	    i = seg_min;
	    icoords[idir] = i;
	    index = d_index(icoords,top_gmax,dim);
	    comp = top_comp[index];
	    n = 0;
	    vst.dens[n+nrad] = m_vst->dens[index];
            vst.engy[n+nrad] = m_vst->engy[index];
            vst.pres[n+nrad] = m_vst->pres[index];
	    for (l = 0; l < dim; ++l)
            	vst.momn[l][n+nrad] = m_vst->momn[(l+idir)%dim][index];
	    for (l = dim; l < 3; ++l)
            	vst.momn[l][n+nrad] = 0.0;
	    seg_max = i;
	    n++;
	    for (i = seg_min+1; i <= imax[idir]; i++)
	    {
		icoords[idir] = i;
		index = d_index(icoords,top_gmax,dim);
		if (needBufferFromIntfc(comp,top_comp[index]))
		    break;
		else
		{
	    	    vst.dens[n+nrad] = m_vst->dens[index];
	    	    vst.engy[n+nrad] = m_vst->engy[index];
	    	    vst.pres[n+nrad] = m_vst->pres[index];
		    for (l = 0; l < dim; ++l)
	    	    	vst.momn[l][n+nrad] = m_vst->momn[(l+idir)%dim][index];
		    for (l = dim; l < 3; ++l)
	    	    	vst.momn[l][n+nrad] = 0.0;
		    n++;
		}
		seg_max = i;
	    }
	    icoords[idir] = seg_min;
	    appendGhostBuffer(&vst,m_vst,n,icoords,idir,0);
	    icoords[idir] = seg_max;
	    appendGhostBuffer(&vst,m_vst,n,icoords,idir,1);
	    
	    eos = &(eqn_params->eos[comp]);
	    EosSetTVDParams(&scheme_params, eos);
	    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
	    n = 0;
	    for (i = seg_min; i <= seg_max; ++i)
	    {
		icoords[idir] = i;
	    	index = d_index(icoords,top_gmax,dim);
	    	m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
	    	m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];
		for (l = 0; l < dim; ++l)
	    	    m_flux->momn_flux[(l+idir)%dim][index] += 
				vflux.momn_flux[l][n+nrad];
		for (l = dim; l < 3; ++l)
	    	    m_flux->momn_flux[l][index] = 0.0;
		n++;
	    }
	    seg_min = seg_max + 1;
	}
        freeDirVstFlux(&vst,&vflux);
}*/	/* end addFluxAlongGridLine */

static void printInputStencil(
        SWEEP vst,
        int n)
{
        int i;
        printf("  density    momn[0]    momn[1]    momn[2]    energy\n");
        for (i = 0; i < n+6; ++i)
        {
            printf("  %7.3g    %7.3g    %7.3g    %7.3g    %7.3g\n",
                    vst.dens[i],vst.momn[0][i],vst.momn[1][i],vst.momn[2][i],
                    vst.engy[i]);
        }    
}       /* end printInputStencil */

void G_CARTESIAN::errFunction()
{
	int i,index;
	double arg,rho,tmp;
	double **vel = field.vel;
	double *dens = field.dens;
	char string[100];
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	double err[3];
	double time = front->time;

	if (FT_Dimension() != 1) return;
	printf("\n");
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'O' || string[0] =='o')
	{
	    if (string[5] == 'A' || string[5] =='a')
	    {
		err[0]=0.0;
		err[1]=0.0;
		err[2]=0.0;
		for (i = 3; i <= top_gmax[0]-3; ++i)
		{
		    index = d_index1d(i,top_gmax);
		    tmp = vel[0][index];
		    arg = top_L[0]+top_h[0]*index - tmp*time;
		    rho = 1.0 + 0.2 * sin(PI*arg);
		    err[0] += fabs(rho-dens[index]);
		    err[1] += pow(rho-dens[index],2);
		    err[2] = (err[2]>fabs(rho-dens[index])) ? err[2] : 
				fabs(rho-dens[index]);
		}
		err[0] /= top_gmax[0]-6;
		err[1] = sqrt(err[1]/(top_gmax[0]-6));
		(void) printf("\n %e \t %e \t %e \n",err[0],err[1],err[2]);
	    }
	}
}	/* end errFunction, check the accuracy in AccuracySineWave case */


//Flood filling
void G_CARTESIAN::adjustGFMStates()
{
    if(eqn_params->tracked)
    {
       double ***Gvel = eqn_params->Gvel;
       double **Gdens = eqn_params->Gdens;
       double **Gpres = eqn_params->Gpres;
       int i,j,k,ii,jj,kk,index,index1,ind;
       int icoords[3];

       for (i = 0; i <= top_gmax[0]; i++)
       for (j = 0; j <= top_gmax[1]; j++)
       for (k = 0; k <= top_gmax[2]; k++)
       {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);

            if (cell_center[index].comp != top_comp[index])
            {
                if (top_comp[index] == GAS_COMP1)
                    ind = 0;
                else
                    ind = 1;

                for (jj = 0; (j+jj) <= top_gmax[1]; jj++)
                {
                     icoords[1] = j+jj;
                     index1 = d_index(icoords,top_gmax,dim);
                     Gdens[ind][index] = Gdens[ind][index1];
                     for (kk = 0; kk < dim; ++kk)
                         Gvel[ind][kk][index] = Gvel[ind][kk][index1];
                     Gpres[ind][index] = Gpres[ind][index1];

                     if (Gdens[ind][index] != 0)
                         break;
                }

                for (jj = 0; (j+jj) >= 0; jj--)
                {
                     icoords[1] = j+jj;
                     index1 = d_index(icoords,top_gmax,dim);
                     Gdens[ind][index] = Gdens[ind][index1];
                     for (kk = 0; kk < dim; ++kk)
                         Gvel[ind][kk][index] = Gvel[ind][kk][index1];
                     Gpres[ind][index] = Gpres[ind][index1];

                     if (Gdens[ind][index] != 0)
                         break;
                }

                for (ii = 0; (i+ii) >= 0; ii--)
                {
                     icoords[0] = i+ii;
                     index1 = d_index(icoords,top_gmax,dim);
                     Gdens[ind][index] = Gdens[ind][index1];
                     for (kk = 0; kk < dim; ++kk)
                         Gvel[ind][kk][index] = Gvel[ind][kk][index1];
                     Gpres[ind][index] = Gpres[ind][index1];

                     if (Gdens[ind][index] != 0)
                         break;
                }

                for (ii = 0; (i+ii) <= top_gmax[0]; ii++)
                {
                     icoords[0] = i+ii;
                     index1 = d_index(icoords,top_gmax,dim);
                     Gdens[ind][index] = Gdens[ind][index1];
                     for (kk = 0; kk < dim; ++kk)
                         Gvel[ind][kk][index] = Gvel[ind][kk][index1];
                     Gpres[ind][index] = Gpres[ind][index1];

                     if (Gdens[ind][index] != 0)
                         break;
                }
            }
       }
    }
}	/* end adjustGFMState */

void G_CARTESIAN::appendOpenEndStates()
{
	INTERFACE *intfc = front->interf;
	int dim = front->rect_grid->dim;
	int i,j,k,idir,side,ii;
	int ic[MAXD], comp, index, bdry_type;
	STATE state;

	if (debugging("trace"))
	    printf("Entering appendOpenEndStates() \n");
	if (dim != 3) return;
	for (idir = 0; idir < dim; ++idir)
        for (side = 0; side < 2; ++side)
        {
            if (rect_boundary_type(intfc,idir,side) == OPEN_BOUNDARY &&
                front->open_end_func != NULL)
            {
		int count = 0;
                for (i = 0; i < top_gmax[0]; ++i)
                for (j = 0; j < top_gmax[1]; ++j)
                for (k = 0; k < top_gmax[2]; ++k)
                {
                    ic[0] = i; ic[1] = j; ic[2] = k;
                    if ((side == 0 && ic[idir] < lbuf[idir]) ||
                        (side == 1 && ic[idir] > top_gmax[idir]-ubuf[idir]))
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        front->open_end_func(front,front->open_end_params,
                                ic,comp,idir,side,&bdry_type,&state);
                        field.dens[index] = state.dens;
                        for (ii = 0; ii < 3; ++ii)
                        {
                            field.vel[ii][index] = state.vel[ii];
                            field.momn[ii][index] = state.momn[ii];
                        }
                        field.pres[index] = state.pres;
                        field.engy[index] = state.engy;
                    }
                }
            }
        }
	if (debugging("trace"))
	    printf("Leaving appendOpenEndStates() \n");
        return;
}	/* end appendOpenEndStates */
