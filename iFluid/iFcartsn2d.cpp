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
 * 	               iFcartsn2d.cpp	
 *******************************************************************/

#include "iFluid.h"

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};
static double (*getStateOldVel[3])(POINTER) = {getStateOldXvel,getStateOldYvel,getStateOldZvel};

//For adding the convective flux to the RHS of the diffusion solver system
void Incompress_Solver_Smooth_2D_Cartesian::computeAdvectionTerm()
{
	int index;
	static HYPERB_SOLVER hyperb_solver(*front);
	static double *rho;
    static int count = 0;

	if (rho == nullptr)
	{
	    int size = (top_gmax[0]+1)*(top_gmax[1]+1);
	    FT_VectorMemoryAlloc((POINTER*)&rho,size,sizeof(double));
	}

	for (int j = 0; j <= top_gmax[1]; j++)
	for (int i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    rho[index] = field->rho[index];
        for (int l = 0; l < dim; l++)
        {
            field->adv_term_old[l][index] = field->adv_term[l][index];
            field->adv_term[l][index] = 0.0;
        }
	}
	hyperb_solver.rho = rho;
    
    count++;
	if (debugging("field_var"))
	{
	    for (int j = jmin; j <= jmax; j++)
	    for (int i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
            field->old_var[0][index] = field->vel[0][index];
            field->old_var[1][index] = field->vel[1][index];
	    }
	}

	hyperb_solver.obst_comp = SOLID_COMP;
	switch (iFparams->adv_order)
	{
	case 1:
	    hyperb_solver.order = 1;
	    hyperb_solver.numericalFlux = upwind_flux;
	    break;
	case 4:
	    hyperb_solver.order = 4;
	    hyperb_solver.numericalFlux = weno5_flux;
	    break;
	default:
	    (void) printf("Advection order %d not implemented!\n",
					iFparams->adv_order);
	    clean_up(ERROR);
	}


    hyperb_solver.adv_term = field->adv_term;
	
    hyperb_solver.var = field->vel;
	hyperb_solver.soln = field->vel;
	
    hyperb_solver.soln_comp1 = LIQUID_COMP1;
	hyperb_solver.soln_comp2 = LIQUID_COMP2;
	
    hyperb_solver.rho1 = iFparams->rho1;
	hyperb_solver.rho2 = iFparams->rho2;

	hyperb_solver.findStateAtCrossing = findStateAtCrossing;
	hyperb_solver.getStateVel[0] = getStateXvel;
	hyperb_solver.getStateVel[1] = getStateYvel;
	
	hyperb_solver.dt = m_dt;
    
    hyperb_solver.computeAdvectionTerm();

	
    if (debugging("field_var"))
	{
	    (void) printf("\nIn computeAdvection(), \n");
	    (void) printf("one step increment for v[0]:\n");
	    computeVarIncrement(field->old_var[0],field->vel[0],NO);
	    (void) printf("one step increment for v[1]:\n");
	    computeVarIncrement(field->old_var[1],field->vel[1],NO);
	    (void) printf("\n");
	}
}


void Incompress_Solver_Smooth_2D_Cartesian::computeAdvection(void)
{
	int i,j,index;
	static HYPERB_SOLVER hyperb_solver(*front);
    static int count = 0;
	
	static double *rho;

	if (rho == nullptr)
	{
	    int size = (top_gmax[0]+1)*(top_gmax[1]+1);
	    FT_VectorMemoryAlloc((POINTER*)&rho,size,sizeof(double));
	}

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    rho[index] = field->rho[index];
	}
	hyperb_solver.rho = rho;
    
    count++;
	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
            field->old_var[0][index] = field->vel[0][index];
            field->old_var[1][index] = field->vel[1][index];
	    }
	}

	hyperb_solver.obst_comp = SOLID_COMP;
	switch (iFparams->adv_order)
	{
	case 1:
	    hyperb_solver.order = 1;
	    hyperb_solver.numericalFlux = upwind_flux;
	    break;
	case 4:
	    hyperb_solver.order = 4;
	    hyperb_solver.numericalFlux = weno5_flux;
	    break;
	default:
	    (void) printf("Advection order %d not implemented!\n",
					iFparams->adv_order);
	    clean_up(ERROR);
	}
	
    //TODO: Don't overwrite soln this way -- limits our flexibility
    hyperb_solver.var = field->vel;
	hyperb_solver.soln = field->vel;
	
    hyperb_solver.soln_comp1 = LIQUID_COMP1;
	hyperb_solver.soln_comp2 = LIQUID_COMP2;
	
    hyperb_solver.rho1 = iFparams->rho1;
	hyperb_solver.rho2 = iFparams->rho2;

	hyperb_solver.findStateAtCrossing = findStateAtCrossing;
	hyperb_solver.getStateVel[0] = getStateXvel;
	hyperb_solver.getStateVel[1] = getStateYvel;
	
	hyperb_solver.dt = m_dt;
    hyperb_solver.solveRungeKutta();

	if (debugging("field_var"))
	{
	    (void) printf("\nIn computeAdvection(), \n");
	    (void) printf("one step increment for v[0]:\n");
	    computeVarIncrement(field->old_var[0],field->vel[0],NO);
	    (void) printf("one step increment for v[1]:\n");
	    computeVarIncrement(field->old_var[1],field->vel[1],NO);
	    (void) printf("\n");
	}
}

void Incompress_Solver_Smooth_2D_Cartesian::computeProjection(void)
{
	switch (iFparams->num_scheme.ellip_method)
	{
	case SIMPLE_ELLIP:
	    if (debugging("check_div"))
                printf("Use computeProjectionSimple()\n");
	    computeProjectionSimple();
	    return;
	case DOUBLE_ELLIP:
	    if (debugging("check_div"))
                printf("Use computeProjectionDouble()\n");
	    computeProjectionDouble();
	    return;
    default:
        printf("Elliptic Method Not Implemented\n");
        clean_up(1);
	}
}	/* end computeProjection */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionSimple(void)
{
	static ELLIPTIC_SOLVER elliptic_solver(front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double **vel_star = field->vel_star;
	double **prev_vel = field->prev_vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double sum_div;
	double value;
	int num_colors;

	sum_div = 0.0;
	max_value = 0.0;
	for (l = 0; l < dim; ++l)
	{
	    vmin[l] = HUGE;
	    vmax[l] = -HUGE;
	}

	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
    		field->old_var[0][index] = field->phi[index];
	    }
	}
	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
        if (!ifluid_comp(top_comp[index]))
        {
            source[index] = 0.0;
            diff_coeff[index] = 0.0;
            continue;
        }

        //source[index] = computeFieldPointDiv(icoords,vel);
        source[index] = computeFieldPointDiv(icoords,vel_star);
        diff_coeff[index] = 1.0/field->rho[index];
        
        if (debugging("check_div"))
	    {
		    for (l = 0; l < dim; ++l)
            {
                if (vmin[l] > field->vel[l][index])
                    vmin[l] = field->vel[l][index];
                if (vmax[l] < field->vel[l][index])
                    vmax[l] = field->vel[l][index];
            }
	    }
	}
    
	if (debugging("field_var"))
	{
	    (void) printf("\nCheck one step increment of div_U:\n");
	    computeVarIncrement(field->div_U,source,NO);
	    (void) printf("\n");
	}
	
    FT_ParallelExchGridArrayBuffer(source,front,NULL);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
	
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
        if (!ifluid_comp(top_comp[index])) continue;
	    
        div_U[index] = source[index];
        source[index] /= accum_dt;

        // Compute pressure jump due to porosity
        icoords[0] = i;
        icoords[1] = j;
        source[index] += computeFieldPointPressureJump(icoords,
                         iFparams->porous_coeff[0],
                         iFparams->porous_coeff[1]);
	    
        array[index] = phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		    index = d_index2d(i,j,top_gmax);
	        if (!ifluid_comp(top_comp[index])) continue;

	        value = fabs(div_U[index]);
            sum_div = sum_div + div_U[index];
            if (max_value < value)
                max_value = value;
	    }

	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}

	if (debugging("check_div"))
    {
        checkVelocityDiv("Before computeProjection()");
    }
        
    //TODO: add vel and prev_vel arrays to elliptic_solver
    //
    elliptic_solver.dt = accum_dt;
    elliptic_solver.D = diff_coeff;
    elliptic_solver.mu = field->mu;
    elliptic_solver.rho = field->rho;
    elliptic_solver.vel = field->vel_star;
    elliptic_solver.prev_vel = field->prev_vel;
    elliptic_solver.source = source;
    elliptic_solver.soln = array;
	elliptic_solver.set_solver_domain();
	elliptic_solver.getStateVar = getStatePhi;
	elliptic_solver.getStateVel[0] = getStateVel[0];
	elliptic_solver.getStateVel[1] = getStateVel[1];
	elliptic_solver.getStateOldVel[0] = getStateOldVel[0];
	elliptic_solver.getStateOldVel[1] = getStateOldVel[1];
	elliptic_solver.findStateAtCrossing = findStateAtCrossing;
	elliptic_solver.skip_neumann_solver = skip_neumann_solver;

	paintAllGridPoint(NOT_SOLVED);
    num_colors = drawColorMap();

	for (i = 1; i < num_colors; ++i)
	{
	    paintToSolveGridPoint2(i);
	    setGlobalIndex();
        setIndexMap();
        elliptic_solver.ij_to_I = ij_to_I;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;

        elliptic_solver.solve(array);
	    paintSolvedGridPoint();
	}

    /*
    while (paintToSolveGridPoint())
    {
        setGlobalIndex();
        setIndexMap();
        elliptic_solver.ij_to_I = ij_to_I;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        
        elliptic_solver.solve(array);
        paintSolvedGridPoint();
    }
    */

	FT_ParallelExchGridArrayBuffer(array,front,NULL);

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    phi[index] = array[index];
	}

	if (debugging("field_var"))
	{
	    printf("\nIn computeProjectionSimple()\n");
	    printf("Check one step increment of phi:\n");
	    computeVarIncrement(field->old_var[0],phi,NO);
	    printf("\n");
	    
	}
}	/* end computeProjectionSimple */

void Incompress_Solver_Smooth_2D_Cartesian::computeProjectionDouble(void)
{
	static DOUBLE_ELLIPTIC_SOLVER elliptic_solver(*front);
	int index;
	int i,j,l,icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double *div_U = field->div_U;
	double sum_div;
	double value;
	int num_colors;

        printf("Entering computeProjectionDouble()\n");
	sum_div = 0.0;
	max_value = 0.0;
	for (l = 0; l < dim; ++l)
	{
	    vmin[l] = HUGE;
	    vmax[l] = -HUGE;
	}

	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index  = d_index2d(i,j,top_gmax);
		field->old_var[0][index] = field->phi[index];
	    }
	}
	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    source[index] = computeFieldPointDiv(icoords,vel);

	    diff_coeff[index] = 1.0/field->rho[index];

	    if (debugging("check_div"))
	    {
		for (l = 0; l < dim; ++l)
		{
		    if (vmin[l] > field->vel[l][index])
			vmin[l] = field->vel[l][index];
		    if (vmax[l] < field->vel[l][index])
			vmax[l] = field->vel[l][index];
		}
	    }
	}
	if (debugging("field_var"))
	{
	    (void) printf("\nCheck one step increment of div_U:\n");
	    computeVarIncrement(field->div_U,source,NO);
	    (void) printf("\n");
	}
	
    FT_ParallelExchGridArrayBuffer(source,front,NULL);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
	
    for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
		continue;
	    source[index] = (source[index])/accum_dt;
            /*Compute pressure jump due to porosity*/
            icoords[0] = i; icoords[1] = j;
            source[index] += computeFieldPointPressureJump(icoords,
                             iFparams->porous_coeff[0],
                             iFparams->porous_coeff[1]);
            /*end of computing pressure jump*/
	    array[index] = phi[index];
	}

	if(debugging("step_size"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		index = d_index2d(i,j,top_gmax);
	        if (!ifluid_comp(top_comp[index]))
		    continue;
	        value = fabs(div_U[index]);
		sum_div = sum_div + div_U[index];
		if (max_value < value)
		    max_value = value;
	    }
	    pp_global_sum(&sum_div,1);
	    (void) printf("\nThe summation of divergence of U is %.16g\n",
					sum_div);
	    pp_global_max(&max_value,1);
	    (void) printf("\nThe max value of divergence of U is %.16g\n",
					max_value);
	    max_value = 0.0;
	}
	if (debugging("check_div"))
        {
	    checkVelocityDiv("Before computeProjection()");
        }
        elliptic_solver.dt = accum_dt;
        elliptic_solver.D = diff_coeff;
        elliptic_solver.source = source;
        elliptic_solver.soln = array;
        elliptic_solver.ext_gmax = ext_gmax;
	elliptic_solver.set_solver_domain();
	elliptic_solver.getStateVar = getStatePhi;
	elliptic_solver.getStateVel[0] = getStateXvel;
	elliptic_solver.getStateVel[1] = getStateYvel;
	elliptic_solver.getStateVel[2] = getStateZvel;
	elliptic_solver.findStateAtCrossing = findStateAtCrossing;
	elliptic_solver.skip_neumann_solver = skip_neumann_solver;
	num_colors = drawColorMap();
	paintAllGridPoint(NOT_SOLVED);
	for (i = 1; i < num_colors; ++i)
	{
	    paintToSolveGridPoint2(i);
	    setDoubleGlobalIndex();
            setDoubleIndexMap();
            elliptic_solver.dij_to_I = dij_to_I;
            elliptic_solver.eilower = eilower;
            elliptic_solver.eiupper = eiupper;
            elliptic_solver.ext_l = ext_l;
            elliptic_solver.ext_u = ext_u;
            elliptic_solver.ext_imin = ext_imin;
            elliptic_solver.ext_imax = ext_imax;
            elliptic_solver.ext_comp = ext_comp;
	    elliptic_solver.set_extension();
	    elliptic_solver.dsolve(array);
	    paintSolvedGridPoint();
	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    phi[index] = array[index];
	}
	if (debugging("field_var"))
	{
	    printf("\nIn computeProjectionSimple()\n");
	    printf("Check one step increment of phi:\n");
	    computeVarIncrement(field->old_var[0],phi,NO);
	    printf("\n");
	    
	}
	return;
}	/* end computeProjectionDouble */

//u^{n+1} = u^{*} - dt*grad(phi^{n+1})
void Incompress_Solver_Smooth_2D_Cartesian::computeNewVelocity(void)
{
	int i,j,k,index;
    double rho;
	COMPONENT comp;
	int icoords[MAXD];
	
    double **vel = field->vel;
	double **vel_star = field->vel_star;
	double **grad_phi = field->grad_phi;
	double *phi = field->phi;

    double point_grad_phi[MAXD];

    /*
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = phi[index];
	}
    */

	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
            vel[0][index] = 0.0;
            vel[1][index] = 0.0;
            continue;
	    }
	    
	    icoords[0] = i;
	    icoords[1] = j;
        rho = field->rho[index];

        computeFieldPointGradJump(icoords,phi,point_grad_phi);
            //computeFieldPointGradJump(icoords,array,point_grad_phi);
 
        for (int l = 0; l < dim; ++l)
        {
            //vel[l][index] -= accum_dt*point_grad_phi[l]/rho;
            vel[l][index] = vel_star[l][index] - accum_dt*point_grad_phi[l]/rho;
            grad_phi[l][index] = point_grad_phi[l];
        }

	}
	
        //TODO: May need to explicitly enforce some tangential boundary
        //      conditions following this velocity update...

    //extractFlowThroughVelocity();
	//computeVelDivergence();
	
    FT_ParallelExchGridVectorArrayBuffer(vel,front);
	FT_ParallelExchGridVectorArrayBuffer(grad_phi,front);

	if (debugging("check_div"))
    {
	    checkVelocityDiv("After computeNewVelocity()");
    }
}	/* end computeNewVelocity */

void Incompress_Solver_Smooth_2D_Cartesian::computeSourceTerm(
	double *coords, 
	double *source) 
{
    for (int i = 0; i < dim; ++i)
        source[i] = iFparams->gravity[i];

	if(iFparams->if_buoyancy)
	{
	    int ic[MAXD];
        rect_in_which(coords,ic,top_grid);
        int index = d_index(ic,top_gmax,dim);
        for (int i = 0; i < dim; ++i)
        {
            source[i] += field->ext_accel[i][index];
        }
	}
}	/* end computeSourceTerm */

void Incompress_Solver_Smooth_2D_Cartesian::solve(double dt)
{
	static boolean first = YES;
	if (first)
	{
	    accum_dt = 0.0;
        if (iFparams->num_scheme.projc_method == PMI ||
            iFparams->num_scheme.projc_method == PMII)
        {
            computeInitialPressure();
            computeGradientQ();
        }
	    first = NO;
	}
	m_dt = dt;

    setFreeStreamVelocity();
	
    start_clock("solve");
    
    setDomain();
	setComponent();
	
    //////////////////////////////////////////////////////
    //TEMP DEBUG
    if (debugging("print_grids"))
    {
        debug_print_grids();
        LOC(); exit(0);
    }
    //////////////////////////////////////////////////////

	paintAllGridPoint(TO_SOLVE);
	setGlobalIndex();

    clearGhostData();

    start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");

	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
    
    if (iFparams->extrapolate_advection)
    {
	    computeAdvectionTerm();
    }
    else
    {
        computeAdvection();
    }

	stop_clock("computeAdvection");
    if (debugging("check_div") || debugging("step_size"))
	{
	    computeMaxSpeed();
	    (void) printf("max_speed after   computeAdvection(): %20.14f ",
				max_speed);
	    (void) printf("occured at (%d, %d)\n",icrds_max[0],icrds_max[1]);
	}

	if (debugging("sample_velocity"))
	    sampleVelocity();
	
    start_clock("computeDiffusion");
	
    computeDiffusion();
    old_dt = m_dt;

	if (debugging("check_div") || debugging("step_size"))
	{
	    computeMaxSpeed();
	    (void) printf("max_speed after   computeDiffusion(): %20.14f ",
				max_speed);
	    (void) printf("occured at (%d, %d)\n",icrds_max[0],icrds_max[1]);
	}
	stop_clock("computeDiffusion");

	if (debugging("sample_velocity"))
	    sampleVelocity();

	// 2) projection step
	accum_dt += m_dt;
	
    if (debugging("step_size"))
	    (void) printf("min_dt = %f  accum_dt = %f\n",min_dt,accum_dt);
	
    if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");

	    start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");
	    
        accum_dt = 0.0;
	}

    recordVelocity();
	computeMaxSpeed();

	if (debugging("sample_velocity"))
	    sampleVelocity();
	if (debugging("step_size"))
	{
	    (void) printf("max_speed after computeNewVelocity(): %20.14f ",
				max_speed);
	    (void) printf("occured at (%d, %d)\n",icrds_max[0],icrds_max[1]);
	}

	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

	/* The following is for climate modelling */
        if(iFparams->if_ref_pres == YES)
            setReferencePressure();

	setAdvectionDt();
	stop_clock("solve");
}	/* end solve */


//TODO: improve boundary handling
double Incompress_Solver_Smooth_2D_Cartesian::getVorticity(int i, int j)
{
        int icoords[MAXD],icnb[MAXD];
        int index,index_nb;
        COMPONENT comp;
        double div,u_edge[2][2];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        int status;
        GRID_DIRECTION dir[2][2] = {{WEST,EAST},{SOUTH,NORTH}};
        int k,idir,nb;
        double u0;
        double dx,dy;
        double vorticity;
        double **vel = field->vel;
        int dim = 2;

        icoords[0] = i;
        icoords[1] = j;

        dx = top_h[0];
        dy = top_h[1];

        index = d_index(icoords,top_gmax,dim);
        comp = top_comp[index];

        for (idir = 0; idir < dim; idir++)
        {
            u0 = vel[(idir+1)%dim][index];
            for (k = 0; k < dim; ++k)
                icnb[k] = icoords[k];
            
            for (nb = 0; nb < 2; nb++)
            {
                icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                index_nb = d_index(icnb,top_gmax,dim);
        
                status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                                comp,&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
                    u_edge[idir][nb] = vel[(idir+1)%dim][index_nb];
                else if (status ==CONST_P_PDE_BOUNDARY)
                    u_edge[idir][nb] = u0;//TODO: should this use the intfc_state vel?
                else if (status ==CONST_V_PDE_BOUNDARY &&
                        wave_type(hs) == DIRICHLET_BOUNDARY)
                    u_edge[idir][nb] = getStateVel[(idir+1)%dim](intfc_state);
                else
                {
                    //TODO: Use a higher order approximation for the
                    //      NEUMANN/MOVABLE_BODY_BOUNDARY (3 point one sided)??
                    //      See computeFieldPointDivSimple().
                    u_edge[idir][nb] = u0;
                }
            }
        }
                    
        vorticity = 0.5*(u_edge[0][1] - u_edge[0][0])/dx
                    - 0.5*(u_edge[1][1] - u_edge[1][0])/dy;

	return vorticity;
}	/* end getVorticity */

void Incompress_Solver_Smooth_2D_Cartesian::copyMeshStates(void)
{
	int i,j,index;
	double **vel = field->vel;
	double *vort = field->vort;
	int symmetry[MAXD];

	double *pres = field->pres;
    double *phi = field->phi;
    double *q = field->q;

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	{
	    index = d_index2d(i,j,top_gmax);
	    if (!ifluid_comp(top_comp[index]))
        {
            pres[index] = 0.0;
            phi[index] = 0.0;
            q[index] = 0.0;
            vort[index] = 0.0;
        }
	    else
        {
            vort[index] = getVorticity(i,j);
        }
	}

    FT_ParallelExchGridArrayBuffer(pres,front,NULL);
    FT_ParallelExchGridArrayBuffer(phi,front,NULL);
    FT_ParallelExchGridArrayBuffer(q,front,NULL);
	
    symmetry[0] = symmetry[1] = ODD;
	FT_ParallelExchGridArrayBuffer(vort,front,symmetry);
}	/* end copyMeshStates */

void Incompress_Solver_Smooth_2D_Cartesian::computeDiffusion(void)
{
    return computeDiffusionCN();
}

//TODO: Write ADI version of this.
void Incompress_Solver_Smooth_2D_Cartesian::computeDiffusionCN(void)
{
    COMPONENT comp;
    int index,index_nb[4],size;
    int I,I_nb[4];
    int i,j,k,l,nb,icoords[MAXD];
    double coords[MAXD], crx_coords[MAXD];
    double coeff[4],mu[4],mu0,rho,rhs;
    
    double source[MAXD];
    double U_nb[4];
    double U_nb_prev[4];
    
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    POINTER intfc_state;
    HYPER_SURF *hs;
    INTERFACE *grid_intfc = front->grid_intfc;
	int status;
    
    PetscInt num_iter;
    double residual;
    double aII;
    double *x;
    
    double** vel = field->vel;
    double** vel_star = field->vel_star;
    double** prev_vel = field->prev_vel;
    double** f_surf = field->f_surf;
    double** grad_q = field->grad_q;
    double** grad_phi = field->grad_phi;

    double **adv_flux = field->adv_term;
    double **adv_flux_old = field->adv_term_old;


    if (debugging("trace"))
        (void) printf("Entering Incompress_Solver_Smooth_2D_Cartesian::"
                    "computeDiffusionCN()\n");
	
    start_clock("computeDiffusionCN");

    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

    for (l = 0; l < dim; ++l)
	{
        PETSc solver;
        solver.Create(ilower, iupper-1, 5, 5);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ij_to_I[i][j];
            if (I == -1) continue;

            index = d_index2d(i,j,top_gmax);
            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);
            
            icoords[0] = i;
            icoords[1] = j;
            comp = top_comp[index];

            I_nb[0] = ij_to_I[i-1][j]; // left or west
            I_nb[1] = ij_to_I[i+1][j]; // right or east
            I_nb[2] = ij_to_I[i][j-1]; // down or south
            I_nb[3] = ij_to_I[i][j+1]; // up or north

            mu0 = field->mu[index];
            rho = field->rho[index];

            for (nb = 0; nb < 4; nb++)
            {
                U_nb_prev[nb] = 0.0;

                int intfc_crx = (*findStateAtCrossing)(front,icoords,
                        dir[nb],comp,&intfc_state,&hs,crx_coords);
                
                if (intfc_crx && wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            U_nb[nb] = getStateVel[l](intfc_state);
                            //TODO: use only tangential component of grad_phi
                            //      at the outlet, or use the whole thing?
                            U_nb[nb] += m_dt*grad_phi[l][index]/rho;
                        }
                        else
                        {
                            //INLET
                            U_nb[nb] = getStateVel[l](intfc_state);
                            //TODO: What about grad_phi for the inlet?
                            /*
                            auto grad_phi_tangent = computeGradPhiTangential(
                                    icoords,dir[nb],comp,hs,crx_coords);
                            U_nb[nb] += m_dt*grad_phi_tangent[l]/rho;
                            */
                        }
                            
                        /*
                        auto grad_phi_tangent = computeGradPhiTangential(
                                icoords,dir[nb],comp,hs,crx_coords);
                        U_nb[nb] += m_dt*grad_phi_tangent[l]/rho;
                        */
                    }
                    else if (neumann_type_bdry(wave_type(hs)))
                    {
                        //TODO: Use flag added to hypersurface
                        //      data structure, no_slip(hs), instead
                        if (iFparams->use_no_slip)
                        {
                            U_nb[nb] = getStateVel[l](intfc_state);
                        }
                        else
                        {
                            /*
                            //TODO: Temporary work around in the case we do not use eddy viscosity
                            //      (do not call computeVremanEddyViscosity() specifically) and the
                            //      ghost_data array does not get populated when setSmoothedProperties()
                            //      is called at beginning of time step.
                            //
                            //      NEW WORKAROUND: now calling computeVelocityGradient() in
                            //      setSmoothedProperties() when eddy viscosity isn't used so that
                            //      the ghost_data array gets populated.
                            
                            if (iFparams->use_eddy_visc == YES)
                            {
                                U_nb[nb] = ghost_data[nb][index].vel[l];
                            }
                            else
                            {
                                //Apply slip boundary condition
                                //nb = 0; idir = 0, nbr = 0;
                                //nb = 1; idir = 0, nbr = 1;
                                //nb = 2; idir = 1, nbr = 0;
                                //nb = 3; idir = 1, nbr = 1;
                                double v_slip[MAXD] = {0.0};
                                int idir = nb/2; int nbr = nb%2;
                                setSlipBoundary(icoords,idir,nbr,comp,hs,intfc_state,field->vel,v_slip);
                                U_nb[nb] = v_slip[l];
                            }
                            */

                            U_nb[nb] = ghost_data[nb][index].vel[l];
                        }
                        
                        auto grad_phi_tangent = computeGradPhiTangential(
                                icoords,dir[nb],comp,hs,crx_coords);
                        U_nb[nb] += m_dt*grad_phi_tangent[l]/rho;
                    }
                    else
                    {
                        printf("Unknown Boundary Type!\n");
                        LOC(); clean_up(EXIT_FAILURE);
                    }


                    if (wave_type(hs) == DIRICHLET_BOUNDARY || neumann_type_bdry(wave_type(hs)))
                        mu[nb] = mu0;
                    else
                        mu[nb] = 0.5*(mu0 + field->mu[index_nb[nb]]);
                
                }
                else
                {
                    U_nb[nb] = vel[l][index_nb[nb]];
                    mu[nb] = 0.5*(mu0 + field->mu[index_nb[nb]]);
                }
            
            }


            coeff[0] = 0.5*m_dt*mu[0]/rho/(top_h[0]*top_h[0]);
            coeff[1] = 0.5*m_dt*mu[1]/rho/(top_h[0]*top_h[0]);
            coeff[2] = 0.5*m_dt*mu[2]/rho/(top_h[1]*top_h[1]);
            coeff[3] = 0.5*m_dt*mu[3]/rho/(top_h[1]*top_h[1]);

            getRectangleCenter(index,coords);
            computeSourceTerm(coords,source);

            //first equation  decoupled, some terms may be lost
            aII = 1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3];
            rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3])*vel[l][index];

            for (nb = 0; nb < 4; nb++)
            {
                status = (*findStateAtCrossing)(front,icoords,dir[nb],comp,
                            &intfc_state,&hs,crx_coords);
       
                if (status == NO_PDE_BOUNDARY)
                {
                    solver.Set_A(I,I_nb[nb],-1.0*coeff[nb]);
                    rhs += coeff[nb]*U_nb[nb];
                }
                else
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {    
                        if (status == CONST_P_PDE_BOUNDARY)
                        {
                            //TODO: outlet not the same at n and n+1
                            //OUTLET
                            rhs += 2.0*coeff[nb]*U_nb[nb];
                        }
                        else
                        {
                            //INLET
                            rhs += 2.0*coeff[nb]*U_nb[nb];
                        }
                    }
                    else if (neumann_type_bdry(wave_type(hs)))
                    {
                        //TODO: This is may be incorrect if a point from
                        //      the previous time step switches component domains
                        //      when the interface was propagated. The coeff and
                        //      the velocity could both potentially be wrong.
                        //
                        //      E.g. When a point goes from fluid comp to solid comp,
                        //      as it is covered by a moving rigid body.
                        //      May need to retain old top_comp array, or find a way
                        //      to access it if that functionality already exists
                        
                        rhs += 2.0*coeff[nb]*U_nb[nb];
                    }
                    else
                    {
                        printf("Unknown Boundary Type!\n");
                        LOC(); clean_up(EXIT_FAILURE);
                    }
                }
            }
          
            rhs += m_dt*source[l];
            rhs += m_dt*f_surf[l][index];

            if (iFparams->num_scheme.projc_method == PMI ||
                iFparams->num_scheme.projc_method == PMII)
            {
                rhs -= m_dt*grad_q[l][index]/rho;
            }
        
            //extrapolate to t^{n+1/2}
            if (iFparams->extrapolate_advection)
            {
                if (old_dt != 0.0)
                {
                    double W0 = -0.5*m_dt/old_dt;
                    double W1 = 1.0 + 0.5*m_dt/old_dt;
                    rhs += W0*adv_flux_old[l][index] + W1*adv_flux[l][index];
                }
                else
                {
                    rhs += adv_flux[l][index];
                }
            }

            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTolerances(1.0e-14,1.0e-12,1.0e06);

	    start_clock("Before Petsc solve");
        solver.Solve();
        solver.GetNumIterations(&num_iter);
        solver.GetResidualNorm(&residual);

	    stop_clock("After Petsc solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            printf("Incompress_Solver_Smooth_2D_Cartesian::"
                    "computeDiffusion: num_iter = %d, residual = %g\n",
                    num_iter,residual);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            index = d_index2d(i,j,top_gmax);
            if (I >= 0)
            {
                vel_star[l][index] = x[I-ilower];
                //vel[l][index] = x[I-ilower];
            }
            else
            {
                vel_star[l][index] = 0.0;
                //vel[l][index] = 0.0;
            }
        }
    
    }
	
    FT_ParallelExchGridVectorArrayBuffer(vel,front);
    FT_FreeThese(1,x);
	
    stop_clock("computeDiffusionCN");
	
    if (debugging("field_var"))
	{
	    (void) printf("\nIn computeDiffusionCN(), \n");
	    (void) printf("one step increment for v[0]:\n");
	    computeVarIncrement(field->old_var[0],vel[0],NO);
	    (void) printf("one step increment for v[1]:\n");
	    computeVarIncrement(field->old_var[1],vel[1],NO);
	    (void) printf("\n");
	}

    if (debugging("trace"))
        (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                    "computeDiffusionCN()\n");
}	/* end computeDiffusionCN */

//TODO: solve equation resulting from taking
//      the divergence of the momentum balance equation ...
//      \nabla^{2} P = 0 ????  Or is there a RHS???
void Incompress_Solver_Smooth_2D_Cartesian::computeInitialPressure()
{
    /*
    printf("ERROR computeInitialPressure(): function not implemented yet!\n");
    LOC(); clean_up(EXIT_FAILURE);
    */
    
    if (debugging("trace"))
            printf("Enterng computeInitialPressure()\n");

    setGlobalIndex();
    setIndexMap();

	int index,index_nb[4];
	double rhs,coeff[4];
	int I,I_nb[4];
	int i,j,l,icoords[MAXD];
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	GRID_DIRECTION opp_dir[4] = {EAST,WEST,NORTH,SOUTH};
	boolean refl_side[4];
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
    int status;
	POINTER intfc_state;
	int icrds_max[MAXD],icrds_min[MAXD];

    double* pres = field->pres;
    double* q = field->q;
    double* rho = field->rho;

	double *x;
	int size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
    for (int ii = 0; ii < size; ++ii) x[ii] = 0.0;

    PETSc solver;
    solver.Create(ilower, iupper-1, 5, 5);

    solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	double max_soln = -HUGE;
	double min_soln = HUGE;

    //set source term
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
        
        if (!ifluid_comp(top_comp[index]))
        {
            source[index] = 0.0;
            diff_coeff[index] = 0.0;
            continue;
        }

        source[index] = 0.0;
        diff_coeff[index] = 1.0/field->rho[index];
    }
    
    FT_ParallelExchGridArrayBuffer(source,front,NULL);
	FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);
	
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
        if (!ifluid_comp(top_comp[index])) continue;
	    
        // Compute pressure jump due to porosity
        icoords[0] = i;
        icoords[1] = j;
        source[index] += computeFieldPointPressureJumpQ(icoords,
                         iFparams->porous_coeff[0],
                         iFparams->porous_coeff[1]);
	}
    //end set source term

    
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    comp = top_comp[index];
	    I = ij_to_I[i][j];
	    if (I == -1) continue;

	    index_nb[0] = d_index2d(i-1,j,top_gmax);
	    index_nb[1] = d_index2d(i+1,j,top_gmax);
	    index_nb[2] = d_index2d(i,j-1,top_gmax);
	    index_nb[3] = d_index2d(i,j+1,top_gmax);
	    
        I_nb[0] = ij_to_I[i-1][j];
	    I_nb[1] = ij_to_I[i+1][j];
	    I_nb[2] = ij_to_I[i][j-1];
	    I_nb[3] = ij_to_I[i][j+1];

        icoords[0] = i;
	    icoords[1] = j;
	
	    num_nb = 0;
	    for (l = 0; l < 4; ++l)
	    {
            status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                    &intfc_state,&hs,crx_coords);

            if (status != CONST_V_PDE_BOUNDARY)
                num_nb++;
                    
            if (status == CONST_V_PDE_BOUNDARY || status == CONST_P_PDE_BOUNDARY)
            {
                index_nb[l] = index;
            }

            coeff[l] = 1.0/(top_h[l/2]*top_h[l/2]);
	    }

	    aII = 0.0;
	    rhs = source[index]*rho[index];

        std::set<int> SetIndices;

	    for (l = 0; l < 4; ++l)
        {
            refl_side[l] = NO;
            if (num_nb == 0) break;

            status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                    &intfc_state,&hs,crx_coords);
            
            if (status == NO_PDE_BOUNDARY)
            {
                if (SetIndices.count(I_nb[l]) == 0)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    SetIndices.insert(I_nb[l]);
                }
                else
                {
                    solver.FlushMatAssembly_A();
                    solver.Add_A(I,I_nb[l],coeff[l]);
                    solver.FlushMatAssembly_A();
                }
                aII -= coeff[l];
            }
            else if (is_bdry_hs(hs) && wave_type(hs) == NEUMANN_BOUNDARY)
            {
                //NEUMANN_BOUNDARY on domain hypersurface bdry
                //  do-nothing
            }
            else if (!is_bdry_hs(hs) && 
                     (wave_type(hs) == NEUMANN_BOUNDARY ||
                      wave_type(hs) == MOVABLE_BODY_BOUNDARY))
            { 
                //grad(phi) dot normal = 0
                int icoords_ghost[MAXD];
                for (int m = 0; m < dim; ++m)
                    icoords_ghost[m] = icoords[m];
                
                int idir = l/2; int nb = l%2;
                icoords_ghost[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                
                double coords_ghost[MAXD];
                double coords_reflect[MAXD];
                
                ////////////////////////////////////////////////////////////////////////
                ///  matches Incompress_Solver_Smooth_Basis::setSlipBoundaryGNOR()  ///
                //////////////////////////////////////////////////////////////////////
                for (int m = 0; m < dim; ++m)
                {
                    coords_ghost[m] = top_L[m] + icoords_ghost[m]*top_h[m];
                    coords_reflect[m] = coords_ghost[m];
                }

                double nor[MAXD];
                FT_NormalAtGridCrossing(front,icoords,
                        dir[l],comp,nor,&hs,crx_coords);
                        
                //Reflect the ghost point through intfc-mirror at crossing.
                //first reflect across the grid line containing intfc crossing,
                coords_reflect[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];
                
                //Reflect the displacement vector across the line
                //containing the intfc normal vector
                double v[MAXD];
                double vn = 0.0;

                for (int m = 0; m < dim; ++m)
                {
                    v[m] = coords_reflect[m] - crx_coords[m];
                    vn += v[m]*nor[m];
                }

                for (int m = 0; m < dim; ++m)
                    v[m] = 2.0*vn*nor[m] - v[m];

                //The desired reflected point
                for (int m = 0; m < dim; ++m)
                    coords_reflect[m] = crx_coords[m] + v[m];
                ////////////////////////////////////////////////////////////////////////
                
                //Interpolate pressure at the reflected point,
                static INTRP_CELL blk_cell;
                static bool first_pres_reflect = true;

                //TODO: This function should only get called once, so checking for the
                //      first call is unneccesary.
                if (first_pres_reflect)
                {
                    const int MAX_NUM_VERTEX_IN_CELL = 20;
                    uni_array(&blk_cell.var,MAX_NUM_VERTEX_IN_CELL,sizeof(double));
                    uni_array(&blk_cell.dist,MAX_NUM_VERTEX_IN_CELL,sizeof(double));
                    uni_array(&blk_cell.coeffs,MAX_NUM_VERTEX_IN_CELL,sizeof(double));
                    bi_array(&blk_cell.coords,MAX_NUM_VERTEX_IN_CELL,MAXD,sizeof(double));
                    bi_array(&blk_cell.icoords,MAX_NUM_VERTEX_IN_CELL,MAXD,sizeof(int));
                    bi_array(&blk_cell.p_lin,MAXD+1,MAXD,sizeof(double));
                    bi_array(&blk_cell.icoords_lin,MAXD+1,MAXD,sizeof(double));
                    uni_array(&blk_cell.var_lin,MAXD+1,sizeof(double));
                    first_pres_reflect = false;
                }
                blk_cell.is_linear = NO;
                blk_cell.is_bilinear = NO;
                
                double pres_reflect;
                FT_IntrpStateVarAtCoordsWithIntrpCoefs(front,&blk_cell,comp,
                        coords_reflect,pres,getStatePres,&pres_reflect,&pres[index]);
                    /*FT_IntrpStateVarAtCoords(front,comp,coords_reflect,pres,
                            getStatePres,&pres_reflect,&pres[index]);*/

                //Place interpolation coefficients of the points used in the
                //approximiation into the system matrix.
                if (blk_cell.is_bilinear)
                {
                    for (int m = 0; m < blk_cell.nv; ++m)
                    {
                        int* ic_intrp = blk_cell.icoords[m];
                        int I_intrp = ij_to_I[ic_intrp[0]][ic_intrp[1]];
                        if (I_intrp == I)
                        {
                           aII += coeff[l]*blk_cell.coeffs[m];
                        }
                        else if (SetIndices.count(I_intrp) == 0)
                        {
                            solver.Set_A(I,I_intrp,coeff[l]*blk_cell.coeffs[m]);
                            SetIndices.insert(I_intrp);
                        }
                        else
                        {
                            solver.FlushMatAssembly_A();
                            solver.Add_A(I,I_intrp,coeff[l]*blk_cell.coeffs[m]);
                            solver.FlushMatAssembly_A();
                        }
                    }
                }
                else if (blk_cell.is_linear)
                {
                    for (int m = 0; m < blk_cell.nv_lin; ++m)
                    {
                        int* ic_intrp = blk_cell.icoords_lin[m];
                        if (ic_intrp[0] == -1)
                        {
                            //move contribution of interface points to the RHS.
                            rhs -= coeff[l]*blk_cell.coeffs[m]*blk_cell.var_lin[m];
                            //TODO: the value at the interface is itself computed
                            //      via interpolation -- should use the interpolating
                            //      points and corresponding coefficients in the matrix
                            //      instead of using this.
                        }
                        else
                        {
                            int I_intrp = ij_to_I[ic_intrp[0]][ic_intrp[1]];
                            if (I_intrp == I)
                            {
                               aII += coeff[l]*blk_cell.coeffs[m];
                            }
                            else if (SetIndices.count(I_intrp) == 0)
                            {
                                solver.Set_A(I,I_intrp,coeff[l]*blk_cell.coeffs[m]);
                                SetIndices.insert(I_intrp);
                            }
                            else
                            {
                                solver.FlushMatAssembly_A();
                                solver.Add_A(I,I_intrp,coeff[l]*blk_cell.coeffs[m]);
                                solver.FlushMatAssembly_A();
                            }
                        }
                    }
                }
                else
                {
                    rhs -= coeff[l]*pres_reflect; 
                }
                
                aII -= coeff[l];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                if (status == CONST_V_PDE_BOUNDARY)
                {
                    //INLET
                    // do-nothing
                
                    /*
                    rhs -= coeff[l]*getStatePres(intfc_state);
                    aII -= coeff[l];
                    use_neumann_solver = NO;
                    */
                }
                else if (status == CONST_P_PDE_BOUNDARY)
                {
                    //OUTLET
                    rhs -= coeff[l]*getStatePres(intfc_state);
                    aII -= coeff[l];
                    use_neumann_solver = NO;
                }
            }
        }

        //TODO: investigate this comment and below if (num_nb > 0) block
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
	    
        if(num_nb > 0)
	    {
            solver.Set_A(I,I,aII);
	    }
        else
        {
            if (debugging("linear_solver"))
                (void) printf("WARNING: isolated value!\n");

            solver.Set_A(I,I,1.0);
            rhs = pres[index];
        }
        solver.Set_b(I,rhs);
	
        SetIndices.clear();
    }


	solver.SetMaxIter(40000);
    solver.SetTolerances(1.0e-10,1.0e-12,1.0e06);
        //solver.SetTolerances(1.0e-14,1.0e-12,1.0e06);

    use_neumann_solver = pp_min_status(use_neumann_solver);
	
	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    if (skip_neumann_solver) 
	    {
            stop_clock("Petsc Solver");
            return;
	    }
	    
        if (debugging("linear_solver"))
	    	(void) printf("\nUsing Neumann Solver!\n");
	    
        //TODO: What is the purpose of this?
        if (size < 20)
	    {
	    	if (debugging("linear_solver"))
                printf("Isolated small region for solve2d()\n");
            stop_clock("Petsc Solver");
            return;
	    }
	    
        solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetResidualNorm(&residual);
	    
        //TODO: skip residual check? GMRES often worse
        if(residual > 1)
	    {
            printf("\n The solution diverges! The residual "
                   "is %g. Solve again using GMRES!\n",residual);
            
            solver.Reset_x();
            solver.Solve_withPureNeumann_GMRES();
            solver.GetNumIterations(&num_iter);
            solver.GetResidualNorm(&residual);
                
            if(residual > 1)
            {
                printf("\n The solution diverges using GMRES. \
                        The residual is %g after %d iterations. Exiting ...\n",
                        residual,num_iter);
                LOC(); clean_up(EXIT_FAILURE);
            }
	    }

	}
	else
	{
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetResidualNorm(&residual);

        //TODO: skip residual check? GMRES often worse
	    if(residual > 1)
	    {
            printf("\n The solution diverges! The residual "
                   "is %g. Solve again using GMRES!\n",residual);
            
            solver.Reset_x();
            solver.Solve_GMRES();
            solver.GetNumIterations(&num_iter);
            solver.GetResidualNorm(&residual);

            if(residual > 1)
            {
                printf("\n The solution diverges using GMRES. \
                        The residual is %g after %d iterations. Exiting ...\n",
                        residual,num_iter);
                LOC(); clean_up(EXIT_FAILURE);
            }
	    }

	}
	stop_clock("Petsc Solver");

	solver.Get_x(x);

	if (debugging("PETSc"))
    {
        printf("In poisson_solver(): "
                "num_iter = %d, residual = %g\n",
                num_iter, residual);
    }

	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];
	    if (I == -1) continue;
	    
        pres[index] = x[I-ilower];
        q[index] = x[I-ilower];
	    
        if (max_soln < pres[index]) 
	    {
            icrds_max[0] = i;
            icrds_max[1] = j;
            max_soln = pres[index];
	    }

	    if (min_soln > pres[index]) 
	    {
            icrds_min[0] = i;
            icrds_min[1] = j;
            min_soln = pres[index];
	    }
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
        printf("Max solution = %20.14f occuring at: %d %d\n",
                max_soln,icrds_max[0],icrds_max[1]);
        //checkSolver(icrds_max,YES);
        
        printf("Min solution = %20.14f occuring at: %d %d\n",
                min_soln,icrds_min[0],icrds_min[1]);
        //checkSolver(icrds_min,YES);
	}

    /*
    if (debugging("elliptic_error"))
    {
        double error,max_error = 0.0;
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            icoords[0] = i;
            icoords[1] = j;
            if (ij_to_I[i][j] == -1) continue;

            error = checkSolver(icoords,NO);
            
            if (error > max_error)
            {
                max_error = error;
                icrds_max[0] = i;
                icrds_max[1] = j;
            }
        }

        printf("In elliptic solver:\n");
        printf("Max relative elliptic error: %20.14f\n",max_error);
        printf("Occuring at (%d %d)\n",icrds_max[0],icrds_max[1]);
        error = checkSolver(icrds_max,YES);
	}
    */

	if (debugging("trace"))
            printf("Leaving computeInitialPressure()\n");

    FT_FreeThese(1,x);
}

//TODO: Just Specify PmI, PmII, PmIII. It's confusing otherwise.
void Incompress_Solver_Smooth_2D_Cartesian::computePressure(void)
{
	switch (iFparams->num_scheme.projc_method)
	{
	case PMI:
	    computePressurePmI();
	    break;
	case PMII:
	    computePressurePmII();
	    break;
	case PMIII:
    case SIMPLE:
	    computePressurePmIII();
	    break;
    /*case SIMPLE:
	    computePressureSimple();
	    break;*/
	case ERROR_PROJC_SCHEME:
	default:
	    (void) printf("Unknown computePressure scheme!\n");
	    clean_up(ERROR);
	}

    if (iFparams->num_scheme.projc_method == PMIII ||
        iFparams->num_scheme.projc_method == SIMPLE) return;

    computeGradientQ();
}

//TODO: PmI and PmII use a lagged pressure term (q) and require solving
//      a poisson problem for the pressure as a startup step.
//      The poisson problem is obtained by taking the divergence of
//      the momentum equation.

// q = p^{n-1/2} , L = I
//  --> p^{n+1/2} = p^{n-1/2} + phi^{n+1}
void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmI(void)
{
    int i,j,index;
	double *rho = field->rho;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;

	for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
	{
        index = d_index2d(i,j,top_gmax);
        pres[index] = q[index] + phi[index];
	    q[index] = pres[index];
	}
    
    //TODO: need to scatter q?
    FT_ParallelExchGridArrayBuffer(pres,front,NULL);
}        /* end computePressurePmI2d */


// q = p^{n-1/2} , L = I - 0.5*nu*dt*grad^2 
//  --> p^{n+1/2} = p^{n-1/2} + phi^{n+1} - 0.5*nu*dt*grad^2(phi^{n+1})
void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmII(void)
{
    int i,j,index;
    double mu0;
	double *rho = field->rho;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
	{
        index = d_index2d(i,j,top_gmax);
        mu0 = field->mu[index];
        pres[index] = q[index] + phi[index] - 0.5*mu0*div_U[index];//If use computeDiffusionCN()
            //pres[index] = q[index] + phi[index] - mu0*div_U[index];//If use computeDiffusionImplicit()
	    q[index] = pres[index];
	}
    
	FT_ParallelExchGridArrayBuffer(pres,front,NULL);
}        /* end computePressurePmII2d */

// q = 0 , L = I - 0.5*nu*dt*grad^2 
//  --> p^{n+1/2} = phi^{n+1} - 0.5*nu*dt*grad^2(phi^{n+1})
void Incompress_Solver_Smooth_2D_Cartesian::computePressurePmIII(void)
{
	double *rho = field->rho;
	double *pres = field->pres;
	double *movie_pres = field->movie_pres;
	double *phi = field->phi;
	double *q = field->q;
	double **grad_q = field->grad_q;
	double *div_U = field->div_U;
    double *mu = field->mu;

	if (debugging("field_var"))
	{
	    for (int j = jmin; j <= jmax; j++)
        for (int i = imin; i <= imax; i++)
	    {
            int index = d_index2d(i,j,top_gmax);
            field->old_var[0][index] = pres[index];
	    }
	}
	
    for (int j = 0; j <= top_gmax[1]; j++)
    for (int i = 0; i <= top_gmax[0]; i++)
	{
        int index = d_index2d(i,j,top_gmax);
        
        //save p^{n-1/2}
        double old_pres = pres[index];

        //compute p^{n+1/2}
        pres[index] = phi[index] - 0.5*mu[index]*div_U[index];//If use computeDiffusionCN()
        //pres[index] = phi[index] - mu[index]*div_U[index];//If use computeDiffusionImplicit()
        
        //record q -- For PmIII q = 0 (for PmI and PmII q = p^{n+1/2})
        q[index] = 0.0;
        for (int l = 0; l < dim; ++l)
        {
            grad_q[l][index] = 0.0;
        }

        //linear extrapolation to get p^{n+1}
        if (m_dt + old_dt != 0.0)
        {
            double W0 = -1.0*m_dt/(m_dt + old_dt);
            double W1 = 1.0 + m_dt/(m_dt + old_dt);
            movie_pres[index] = W0*old_pres + W1*pres[index];
        }
	}

	FT_ParallelExchGridArrayBuffer(pres,front,NULL);
	FT_ParallelExchGridArrayBuffer(q,front,NULL);

	if (debugging("field_var"))
	{
	    (void) printf("\nCheck one step increment of Pressure:\n");
	    computeVarIncrement(field->old_var[0],pres,NO);
	    (void) printf("\n");
	}
}        /* end computePressurePmIII */

/*
void Incompress_Solver_Smooth_2D_Cartesian::computePressureSimple(void)
{
    int i,j,index;
    double mu0;
	double *rho = field->rho;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	    {
            index = d_index2d(i,j,top_gmax);
            field->old_var[0][index] = pres[index];
	    }
	}
	
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
	{
        index = d_index2d(i,j,top_gmax);
        mu0 = field->mu[index];
        pres[index] = phi[index] - mu0*div_U[index];//If use computeDiffusionImplicit()
            //pres[index] = phi[index] - 0.5*mu0*div_U[index];//If use computeDiffusionCN()
                //pres[index] = phi[index];
        q[index] = 0.0;
	}

	FT_ParallelExchGridArrayBuffer(pres,front,NULL);
	FT_ParallelExchGridArrayBuffer(q,front,NULL);

	if (debugging("field_var"))
	{
	    (void) printf("\nCheck one step increment of Pressure:\n");
	    computeVarIncrement(field->old_var[0],pres,NO);
	    (void) printf("\n");
	}
}*/        /* end computePressureSimple */

void Incompress_Solver_Smooth_2D_Cartesian::computeGradientQ()
{
	int i,j,l,index;
	int icoords[MAXD];
	double point_grad_q[MAXD];
	double **grad_q = field->grad_q;
    double *q = field->q;

	if (debugging("field_var"))
	{
	    int comp;
	    for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
	    	comp = top_comp[index];
	    	if (!ifluid_comp(comp)) continue;

	    	icoords[0] = i;
	    	icoords[1] = j;
		
            field->old_var[0][index] = field->grad_q[0][index];
		    field->old_var[1][index] = field->grad_q[1][index];
	    }
	}

	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index2d(i,j,top_gmax);
	    array[index] = q[index];
	}
	
    for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    
        computeFieldPointGradJumpQ(icoords,array,point_grad_q);
	    
        for (l = 0; l < dim; ++l)
	    	grad_q[l][index] = point_grad_q[l];
	}
	
    FT_ParallelExchGridVectorArrayBuffer(grad_q,front);
	
    if (debugging("field_var"))
	{
	    (void) printf("\nIn computeGradientQ(), \n");
	    (void) printf("one step increment for grad_phi[0]:\n");
	    computeVarIncrement(field->old_var[0],field->grad_q[0],NO);
	    (void) printf("one step increment for grad_phi[1]:\n");
	    computeVarIncrement(field->old_var[1],field->grad_q[1],NO);
	    printf("\n");
	}
}	/* end computeGradientQ */

void Incompress_Solver_Smooth_2D_Cartesian::surfaceTension(
	double *coords,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma)
{
	int i,k,nb;
	BOND *pb,*bonds[100];
	double kappa,nor[MAXD];
	double kappa0,nor0[MAXD];
	double kappa1,nor1[MAXD];
	double len,delta;
	HYPER_SURF_ELEMENT *phse;
	double p[MAXD];
	

	BondAndNeighbors(hse,hs,&nb,bonds,3);

	for (i = 0; i < dim; ++i) force[i] = 0.0;
	for (k = 0; k < nb; ++k)
	{
	    pb = bonds[k];
	    for (i = 0; i < dim; ++i) 
		p[i] = 0.5*(Coords(pb->start)[i] + Coords(pb->end)[i]);
	    delta = smoothedDeltaFunction(coords,p);
	    if (delta == 0.0) continue;

	    len = bond_length(pb);
	    phse = Hyper_surf_element(pb);
	    GetFrontNormal(pb->start,phse,hs,nor0,front);
	    GetFrontNormal(pb->end,phse,hs,nor1,front);
	    for (i = 0; i < dim; ++i) nor[i] = 0.5*(nor0[i] + nor1[i]);
	    GetFrontCurvature(pb->start,phse,hs,&kappa0,front);
	    GetFrontCurvature(pb->end,phse,hs,&kappa1,front);
	    kappa = 0.5*(kappa0 + kappa1);
	    for (i = 0; i < dim; ++i) 
	    {
		force[i] += delta*sigma*len*kappa*nor[i];
	    }
	}
}	/* end surfaceTension2d */

void Incompress_Solver_Smooth_2D_Cartesian::setInitialCondition()
{
	COMPONENT comp;
	double coords[MAXD];
	int size = (int)cell_center.size();

	FT_MakeGridIntfc(front);
	setDomain();

    m_rho[0] = iFparams->rho1;
    m_rho[1] = iFparams->rho2;
    m_mu[0] = iFparams->mu1;
    m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;

	for (int i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

    double* mu = field->mu;
    double* rho = field->rho;

    for (int index = 0; index < size; index++)
    {
        mu[index] = 0.0;
        rho[index] = 0.0;
        
        comp = top_comp[index];
        if (!ifluid_comp(comp)) continue;

        switch (comp)
        {
            case LIQUID_COMP1:
                mu[index] = m_mu[0];
                rho[index] = m_rho[0];
                break;
            case LIQUID_COMP2:
                mu[index] = m_mu[1];
                rho[index] = m_rho[1];
                break;
        }
    }

    double *pres = field->pres;
    double *phi = field->phi;
    double *q = field->q;

	// Initialize state at cell_center
    for (int i = 0; i < size; i++)
    {
        getRectangleCenter(i,coords);
        comp = top_comp[i];
        
        if (getInitialState != NULL)
        {
            (*getInitialState)(comp,coords,field,i,dim,iFparams);
            q[i] = getQFromPres(front,pres[i]);
        }

        //TODO: see comments in getPressure() function
            //pres[i] = getPressure(front,coords,NULL);
            //phi[i] = getPhiFromPres(front,pres[i]);
    }

	computeGradientQ();
    copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

void Incompress_Solver_Smooth_2D_Cartesian::setParallelVelocity()
{
        FILE *infile;
        int i,j,id,k,l,index,G_index;
        char fname[100];
        COMPONENT comp;
        double coords[MAXD];
        int size = (int)cell_center.size();
	int myid = pp_mynode();
	int numprocs = pp_numnodes();

        int G_icoords[MAXD],pp_icoords[MAXD],icoords[MAXD];
        int local_gmax[MAXD], global_gmax[MAXD];
        int G_size, L_size;
	PP_GRID *pp_grid = front->pp_grid;
	double *local_L = pp_grid->Zoom_grid.L;
	double *local_U = pp_grid->Zoom_grid.U;
	double *GU_buff,*GV_buff, *U_buff, *V_buff;
        for (i = 0; i < dim; i++)
        {
            global_gmax[i] = pp_grid->Global_grid.gmax[i]-1;
            local_gmax[i] = pp_grid->Zoom_grid.gmax[i]-1;
        }
        FT_MakeGridIntfc(front);
        setDomain();
        G_size = 1;
        L_size = 1;
	for (i = 0; i < dim; i++)
	{
	    G_size = G_size * (global_gmax[i]+1);
	    L_size = L_size * (top_gmax[i]+1);
	}
	uni_array(&U_buff,L_size,sizeof(double));
	uni_array(&V_buff,L_size,sizeof(double));
	if (myid == 0)
	{
	    uni_array(&GU_buff,G_size,sizeof(double));
	    uni_array(&GV_buff,G_size,sizeof(double));
	 
	    if (setInitialVelocity != NULL)
                (*setInitialVelocity)(comp,pp_grid->Global_grid.gmax,
				   GU_buff,GV_buff,NULL,
				    &(pp_grid->Global_grid),iFparams);
	    for (id = 0; id < numprocs; id++)
	    {            
		find_Cartesian_coordinates(id,pp_grid,pp_icoords);
		for (j = jmin; j <= jmax; ++j)
            	for (i = imin; i <= imax; ++i)
            	{
	            icoords[0] = i;
	            icoords[1] = j;
	            G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
	            G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
	            G_index = d_index(G_icoords,global_gmax,dim);
	            index = d_index(icoords,top_gmax,dim);
	            U_buff[index] = GU_buff[G_index];
	            V_buff[index] = GV_buff[G_index];
	        }
		if (id == 0)
		{	
		    for (i = 0; i < L_size; i++)
		    {
			field->vel[0][i] = U_buff[i];
			field->vel[1][i] = V_buff[i];
		    }
		}
		else
		{
	            pp_send(1,(POINTER)(U_buff),sizeof(double)*L_size,id);
	            pp_send(2,(POINTER)(V_buff),sizeof(double)*L_size,id);
		}
	    }
	    FT_FreeThese(2,GU_buff,GV_buff);
	}
	else
	{
	    pp_recv(1,0,(POINTER)(U_buff),sizeof(double)*L_size);
	    pp_recv(2,0,(POINTER)(V_buff),sizeof(double)*L_size);
	    for (i = 0; i < L_size; i++)
	    {
	        field->vel[0][i] = U_buff[i];
		field->vel[1][i] = V_buff[i];
	    }
	}


        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
        m_comp[0] = iFparams->m_comp1;
        m_comp[1] = iFparams->m_comp2;
        m_smoothing_radius = iFparams->smoothing_radius;
        m_sigma = iFparams->surf_tension;
        mu_min = rho_min = HUGE;
        for (i = 0; i < 2; ++i)
        {
            if (ifluid_comp(m_comp[i]))
            {
                mu_min = std::min(mu_min,m_mu[i]);
                rho_min = std::min(rho_min,m_rho[i]);
            }
        }
	FT_FreeThese(2,U_buff,V_buff);

	computeGradientQ();
        copyMeshStates();
        setAdvectionDt();
}       /* end setParallelVelocity */

/*
void Incompress_Solver_Smooth_2D_Cartesian::
        computeDiffusionImplicit(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
        int i,j,k,l,nb,icoords[MAXD];
        double coords[MAXD], crx_coords[MAXD];
        double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double *x;
        GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
        POINTER intfc_state;
        HYPER_SURF *hs;
        PetscInt num_iter;
        double rel_residual;
        double aII;
        double source[MAXD];
        double **vel = field->vel;
        double **f_surf = field->f_surf;
        INTERFACE *grid_intfc = front->grid_intfc;

        if (debugging("trace"))
            (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionImplicit()\n");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

        for (l = 0; l < dim; ++l)
        {
            PETSc solver;
            solver.Create(ilower, iupper-1, 5, 5);
            solver.Reset_A();
            solver.Reset_b();
            solver.Reset_x();

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I  = ij_to_I[i][j];
                if (I == -1) continue;

                index  = d_index2d(i,j,top_gmax);
                index_nb[0] = d_index2d(i-1,j,top_gmax);
                index_nb[1] = d_index2d(i+1,j,top_gmax);
                index_nb[2] = d_index2d(i,j-1,top_gmax);
                index_nb[3] = d_index2d(i,j+1,top_gmax);

                icoords[0] = i;
                icoords[1] = j;
                comp = top_comp[index];

                I_nb[0] = ij_to_I[i-1][j]; //west
                I_nb[1] = ij_to_I[i+1][j]; //east
                I_nb[2] = ij_to_I[i][j-1]; //south
                I_nb[3] = ij_to_I[i][j+1]; //north


                mu0   = field->mu[index];
                rho   = field->rho[index];

                for (nb = 0; nb < 4; nb++)
                {
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
                                dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                    {
                        if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                        {
                            U_nb[nb] = vel[l][index];
                        }
                        else
                        {
                            U_nb[nb] = getStateVel[l](intfc_state);
                        }
                        if (wave_type(hs) == DIRICHLET_BOUNDARY ||
                            neumann_type_bdry(wave_type(hs)))
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

                coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
                coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
                coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
                coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);

                getRectangleCenter(index, coords);
                computeSourceTerm(coords, source);

                aII = 1 + 2.0*(coeff[0]+coeff[1]+coeff[2]+coeff[3]);
                rhs = vel[l][index];

                for(nb = 0; nb < 4; nb++)
                {
                    if (!(*findStateAtCrossing)(front,icoords,dir[nb],comp,
                                &intfc_state,&hs,crx_coords))
                        solver.Set_A(I,I_nb[nb],-coeff[nb]);
                    else
                    {
                        if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                            boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                            aII -= coeff[nb];
                        else
                            rhs += coeff[nb]*U_nb[nb];
                    }
                }
                rhs += m_dt*source[l];
                rhs += m_dt*f_surf[l][index];
                //rhs -= m_dt*grad_q[l][index]/rho;
                solver.Set_A(I,I,aII);
                solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1.e-10);

            start_clock("Befor Petsc solve");
            //solver.Solve_GMRES();
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            stop_clock("After Petsc solve");
            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
                        "computeDiffusionImplicit: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
        FT_FreeThese(1,x);

        if (debugging("trace"))
            (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionImplicit()\n");
}*/       /* end computeDiffusionImplicit */

void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusionImplicit(void)
{
    COMPONENT comp;
    int index,index_nb[4],size;
    int I,I_nb[4];
    int i,j,k,l,nb,icoords[MAXD];
    double coords[MAXD], crx_coords[MAXD];
    double coeff[4],mu[4],mu0,rho,rhs;
    
    double source[MAXD];
    double U_nb[4];
    
    GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
    POINTER intfc_state;
    HYPER_SURF *hs;
    INTERFACE *grid_intfc = front->grid_intfc;
	int status;
    
    PetscInt num_iter;
    double residual;
    double aII;
    double *x;
    
    double **vel = field->vel;
    double **f_surf = field->f_surf;
    double **grad_q = field->grad_q;

    if (debugging("trace"))
        (void) printf("Entering Incompress_Solver_Smooth_2D_Cartesian::"
                    "computeDiffusionImplicit()\n");
	
    start_clock("computeDiffusionImplicit");

    setIndexMap();

    size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	if (debugging("field_var"))
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	index = d_index2d(i,j,top_gmax);
		    field->old_var[0][index] = vel[0][index];
		    field->old_var[1][index] = vel[1][index];
	    }
	}
	
    for (l = 0; l < dim; ++l)
	{
        PETSc solver;
        solver.Create(ilower, iupper-1, 5, 5);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ij_to_I[i][j];
            if (I == -1) continue;

            index  = d_index2d(i,j,top_gmax);
            index_nb[0] = d_index2d(i-1,j,top_gmax);
            index_nb[1] = d_index2d(i+1,j,top_gmax);
            index_nb[2] = d_index2d(i,j-1,top_gmax);
            index_nb[3] = d_index2d(i,j+1,top_gmax);
            
            icoords[0] = i;
            icoords[1] = j;
            comp = top_comp[index];

            I_nb[0] = ij_to_I[i-1][j]; // left or west
            I_nb[1] = ij_to_I[i+1][j]; // right or east
            I_nb[2] = ij_to_I[i][j-1]; // down or south
            I_nb[3] = ij_to_I[i][j+1]; // up or north

            mu0 = field->mu[index];
            rho = field->rho[index];

            for (nb = 0; nb < 4; nb++)
            {
                int intfc_crx = (*findStateAtCrossing)(front,icoords,
                        dir[nb],comp,&intfc_state,&hs,crx_coords);
                
                if (intfc_crx && wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            U_nb[nb] = getStateVel[l](intfc_state);
                        }
                        else
                        {
                            //INLET
                            U_nb[nb] = getStateVel[l](intfc_state);
                        }
                    }
                    else if (neumann_type_bdry(wave_type(hs)))
                    {
                        if (!is_bdry_hs(hs))//TODO: handle another way -- we want to include these (see below)
                        {
                            //TODO: Add option to use no-slip boundary instead where
                            //          vel_nb[nb] = intfc_state->vel[l];
                            //      Need to add a flag to hypersurface data structure
                            
                            //Apply slip boundary condition
                            //nb = 0; idir = 0, nbr = 0;
                            //nb = 1; idir = 0, nbr = 1;
                            //nb = 2; idir = 1, nbr = 0;
                            //nb = 3; idir = 1, nbr = 1;
                            double v_slip[MAXD] = {0.0};
                            int idir = nb/2; int nbr = nb%2;
                            setSlipBoundary(icoords,idir,nbr,comp,hs,intfc_state,field->vel,v_slip);
                            U_nb[nb] = v_slip[l];
                        }
                        else
                        {
                            //TODO: Without this rayleigh-taylor with NEUMANN boundaries
                            //      crashes for some reason.
                            U_nb[nb] = getStateVel[l](intfc_state);
                            //NOTE: This is just a no-slip boundary condition.
                        }
                    }
                    else if (wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        //Same as NO_PDE_BOUNDARY
                        U_nb[nb] = vel[l][index_nb[nb]];
                        mu[nb] = 0.5*(mu0 + field->mu[index_nb[nb]]);
                    }
                    else
                    {
                        printf("Unknown Boundary Type!\n");
                        LOC(); clean_up(EXIT_FAILURE);
                    }

                    if (wave_type(hs) == DIRICHLET_BOUNDARY || neumann_type_bdry(wave_type(hs)))
                        mu[nb] = mu0;
                    else
                        mu[nb] = 0.5*(mu0 + field->mu[index_nb[nb]]);
                }
                else
                {
                    //NO_PDE_BOUNDARY
                    U_nb[nb] = vel[l][index_nb[nb]];
                    mu[nb] = 0.5*(mu0 + field->mu[index_nb[nb]]);
                }
            }


            coeff[0] = m_dt*mu[0]/rho/(top_h[0]*top_h[0]);
            coeff[1] = m_dt*mu[1]/rho/(top_h[0]*top_h[0]);
            coeff[2] = m_dt*mu[2]/rho/(top_h[1]*top_h[1]);
            coeff[3] = m_dt*mu[3]/rho/(top_h[1]*top_h[1]);

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, source);

            //first equation  decoupled, some terms may be lost
            aII = 1.0 + 2.0*(coeff[0] + coeff[1] + coeff[2] + coeff[3]);
            rhs = vel[l][index];

            for (nb = 0; nb < 4; nb++)
            {
                status = (*findStateAtCrossing)(front,icoords,dir[nb],comp,
                            &intfc_state,&hs,crx_coords);
       
                if (status == NO_PDE_BOUNDARY)
                {
                    solver.Set_A(I,I_nb[nb],-1.0*coeff[nb]);
                }
                else
                {
                    if (wave_type(hs) == ELASTIC_BOUNDARY)
                    {
                        //Same as NO_PDE_BOUNDARY
                        solver.Set_A(I,I_nb[nb],-1.0*coeff[nb]);
                    }
                    else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {    
                        if (status == CONST_P_PDE_BOUNDARY)
                        {
                            //OUTLET
                            rhs += coeff[nb]*U_nb[nb];
                                //rhs += coeff[nb]*U_nb[nb]; //u^n val
                                //solver.Set_A(I,I_nb[nb],-coeff[nb]);
                                //aII -= coeff[nb];
                                //rhs += coeff[nb]*U_nb[nb];
                        }
                        else
                        {
                            //INLET
                            rhs += coeff[nb]*U_nb[nb];
                        }
                    }
                    else if (neumann_type_bdry(wave_type(hs)))
                    {
                        //NEUMANN
                        rhs += coeff[nb]*U_nb[nb];
                    }
                    else
                    {
                        printf("Unknown Boundary Type!\n");
                        LOC(); clean_up(EXIT_FAILURE);
                    }
                }
            }
          
            rhs += m_dt*source[l];
            rhs += m_dt*f_surf[l][index];

            if (iFparams->num_scheme.projc_method != PMIII &&
                iFparams->num_scheme.projc_method != SIMPLE)
            {
                rhs -= m_dt*grad_q[l][index]/rho;
            }

            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTolerances(1.0e-14,1.0e-12,1.0e06);

	    start_clock("Before Petsc solve");
        solver.Solve();
        solver.GetNumIterations(&num_iter);
        solver.GetResidualNorm(&residual);

	    stop_clock("After Petsc solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            printf("Incompress_Solver_Smooth_2D_Cartesian::"
                    "computeDiffusion: num_iter = %d, residual = %g\n",
                    num_iter,residual);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ij_to_I[i][j];
            index = d_index2d(i,j,top_gmax);
            if (I >= 0)
                vel[l][index] = x[I-ilower];
            else
                vel[l][index] = 0.0;
        }
    
    }
	
    FT_ParallelExchGridVectorArrayBuffer(vel,front);
    FT_FreeThese(1,x);
	
    stop_clock("computeDiffusionImplicit");
	
    if (debugging("field_var"))
	{
	    (void) printf("\nIn computeDiffusionCN(), \n");
	    (void) printf("one step increment for v[0]:\n");
	    computeVarIncrement(field->old_var[0],vel[0],NO);
	    (void) printf("one step increment for v[1]:\n");
	    computeVarIncrement(field->old_var[1],vel[1],NO);
	    (void) printf("\n");
	}

    if (debugging("trace"))
        (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                    "computeDiffusionImplicit()\n");
}	/* end computeDiffusionImplicit */

static int parab_find_state_at_crossing(
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

        switch (wave_type(*hs))
        {
        case FIRST_PHYSICS_WAVE_TYPE:
            return NO_PDE_BOUNDARY;
        case DIRICHLET_BOUNDARY:
        case GROWING_BODY_BOUNDARY:
        case MOVABLE_BODY_BOUNDARY:
        case NEUMANN_BOUNDARY:
        case ICE_PARTICLE_BOUNDARY:
            return DIRICHLET_PDE_BOUNDARY;
        default:
            (void) printf("In parab_find_state_at_crossing()\n");
            (void) printf("Unknown wave type %s\n",
                                f_wave_type_as_string(wave_type(*hs)));
            clean_up(ERROR);
        }
}       /* end parab_find_state_at_crossing */

void Incompress_Solver_Smooth_2D_Cartesian::
        computeDiffusionParab(void)
{
        static PARABOLIC_SOLVER parab_solver(*front);

        COMPONENT comp;
        int index;
        int i,j,l,icoords[MAXD];
        double source[MAXD];
        double **vel = field->vel;
        double **f_surf = field->f_surf;
        static double *nu;
	static boolean first = YES;

        if (debugging("trace"))
            (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionParab()\n");

        if (nu == NULL)
        {
            int size = 1;
            for (i = 0; i < dim; ++i)  size *= (top_gmax[i] + 1);
            FT_VectorMemoryAlloc((POINTER*)&nu,size,sizeof(double));
        }
        setIndexMap();
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index2d(i,j,top_gmax);
            nu[index] = field->mu[index]/field->rho[index];
        }
        FT_ParallelExchGridArrayBuffer(nu,front,NULL);

        parab_solver.soln_comp = LIQUID_COMP2;
        parab_solver.obst_comp = SOLID_COMP;
        parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
        parab_solver.dt = m_dt;
        parab_solver.order = 2;
        parab_solver.a = NULL;
        parab_solver.findStateAtCrossing = parab_find_state_at_crossing;
        parab_solver.first = first;
        parab_solver.set_solver_domain();
	first = NO;
        switch(dim)
        {
        case 2:
            parab_solver.ij_to_I = ij_to_I;
            break;
        case 3:
            parab_solver.ijk_to_I = ijk_to_I;
            break;
        }
        for (l = 0; l < dim; ++l)
        {
            parab_solver.var = vel[l];
            parab_solver.soln = vel[l];
            parab_solver.getStateVarFunc = getStateVel[l];
            parab_solver.solveIM();
        }

        if (debugging("trace"))
            (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
                        "computeDiffusionParab()\n");
}       /* end computeDiffusionParab */

void Incompress_Solver_Smooth_2D_Cartesian::
	computeDiffusionExplicit(void)
{
        COMPONENT comp;
        int index,index_nb[4],size;
        int I,I_nb[4];
	int i,j,l,nb,icoords[MAXD];
	double coords[MAXD], crx_coords[MAXD];
	double coeff[4],mu[4],mu0,rho,rhs,U_nb[4];
        double *x;
	GRID_DIRECTION dir[4] = {WEST,EAST,SOUTH,NORTH};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double aII;
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusionExplicit()\n");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ij_to_I[i][j];
            	if (I == -1) continue;

            	index  = d_index2d(i,j,top_gmax);
            	index_nb[0] = d_index2d(i-1,j,top_gmax);
            	index_nb[1] = d_index2d(i+1,j,top_gmax);
            	index_nb[2] = d_index2d(i,j-1,top_gmax);
            	index_nb[3] = d_index2d(i,j+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		comp = top_comp[index];

            	mu0   = field->mu[index];
            	rho   = field->rho[index];

            	for (nb = 0; nb < 4; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
                    	    U_nb[nb] = vel[l][index];
			else
			    U_nb[nb] = getStateVel[l](intfc_state);
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    neumann_type_bdry(wave_type(hs)))
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

            	coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		rhs = (-coeff[0]-coeff[1]-coeff[2]-coeff[3])*vel[l][index];

		int num_nb = 0;
		for(nb = 0; nb < 4; nb++)
		{
		    rhs += coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
                x[I-ilower] = vel[l][index] + rhs;
            }

            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ij_to_I[i][j];
                index = d_index2d(i,j,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
        FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_2D_Cartesian::"
			"computeDiffusionExplicit()\n");
}       /* end computeDiffusionExplicit */

void Incompress_Solver_Smooth_2D_Cartesian::vtk_plot_scalar(
        char *outname, const char* varname)
{
}	/* end vtk_plot_scalar */

void Incompress_Solver_Smooth_2D_Cartesian::extractFlowThroughVelocity()
{
	int index,index_nb,index_op;
	int i,j,k,l,nb,icoords[MAXD],icoords_nb[MAXD],icoords_op[MAXD];
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	POINTER intfc_state;
        HYPER_SURF *hs;
	double crx_coords[MAXD];
	double vel_save,**vel = field->vel;
	double div;
	int status;

	/* Extract velocity for zero interior divergence */
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
	{
	    if (i >= imin && i <= imax && j >= jmin && j <= jmax)
		continue;

	    icoords[0] = i;
	    icoords[1] = j;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index])) continue;
	    for (l = 0; l < dim; ++l)
		icoords_nb[l] = icoords_op[l] = icoords[l];
	    for (l = 0; l < dim; ++l)
	    for (nb = 0; nb < 2; ++nb)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l][nb],
                                top_comp[index],&intfc_state,&hs,crx_coords);
		if (status == NO_PDE_BOUNDARY) continue;
		icoords_nb[l] = (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;
		icoords_op[l] = (nb == 0) ? icoords[l] - 2 : icoords[l] + 2;
            	index_nb = d_index(icoords_nb,top_gmax,dim);
            	index_op = d_index(icoords_op,top_gmax,dim);
		if (!ifluid_comp(top_comp[index_nb]))
		    continue;
		vel[l][index] = 0.0;
                div = computeFieldPointDiv(icoords_nb,vel);
                vel[l][index] = (nb == 0) ? vel[l][index] - 2.0*top_h[l]*div
                                : vel[l][index] + 2.0*top_h[l]*div;
	    }
	}
}	/* end extractFlowThroughVelocity */

void Incompress_Solver_Smooth_2D_Cartesian::computeVarIncrement(
	double *var_old,
	double *var_new,
	boolean use_dual_grid)
{
	int i,j,k,index,size;
	double mag;
	double Lnorm[3];

	if (use_dual_grid)
	{
	    ;	// To add dual grid computation.
	}
	else
	{
	    Lnorm[0] = Lnorm[1] = Lnorm[2] = 0.0;
	    size = (imax - imin + 1)*(jmax - jmin + 1);
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index2d(i,j,top_gmax);
		mag = fabs(var_new[index] - var_old[index]);
		Lnorm[0] += mag;
		Lnorm[1] += sqr(mag);
		if (Lnorm[2] < mag)
		{
		    Lnorm[2] = mag;
		}
            }
	    Lnorm[0] /= (imax - imin + 1)*(jmax - jmin + 1);
	    Lnorm[1] = sqrt(Lnorm[1]/size);
	    (void) printf("L-1 norm = %20.14f  L-1/dt = %20.14f\n",
					Lnorm[0],Lnorm[0]/m_dt);
	    (void) printf("L-2 norm = %20.14f  L-2/dt = %20.14f\n",
					Lnorm[1],Lnorm[1]/m_dt);
	    (void) printf("L-I norm = %20.14f  L-I/dt = %20.14f\n",
					Lnorm[2],Lnorm[2]/m_dt);
	    Lnorm[0] = 0.0;
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		mag = fabs(var_new[index]);
		Lnorm[0] += mag;
	    }
	    Lnorm[0] /= (imax - imin + 1)*(jmax - jmin + 1);
	    (void) printf("L-1 norm of old variable = %20.14f\n",Lnorm[0]);
	}
}	/* end computeVarIncrement */	

void Incompress_Solver_Smooth_2D_Cartesian::computeVelDivergence()
{
	double *div_U = field->div_U;
	double **vel = field->vel;
	int i,j,index,icoords[MAXD];

	/* Compute velocity divergence */
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    index  = d_index(icoords,top_gmax,dim);
	    if (!ifluid_comp(top_comp[index]))
            div_U[index] = 0.0;
	    div_U[index] = computeFieldPointDiv(icoords,vel);
	}
}	/* end computeVelDivergence */
