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

#include "solver.h"

ELLIPTIC_SOLVER::ELLIPTIC_SOLVER(Front* fr)
    : front(fr)
{
    //TODO: remove after modifying solve2d()
    porosity = 0.0;
}

void ELLIPTIC_SOLVER::set_solver_domain(void)
{

	static boolean first = YES;
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
    struct Table *T = table_of_interface(front->grid_intfc);
    int *lbuf = front->rect_grid->lbuf;
    int *ubuf = front->rect_grid->ubuf;
	int i;

	dim = Dimension(front->interf);
    top_comp = T->components;
    top_gmax = rgr->gmax;
	top_h = rgr->h;
	top_L = rgr->L;
	
    if (first)
	{
	    first = NO;
	    switch(dim)
	    {
	    case 1:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
		FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
		break;
	    case 2:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		FT_VectorMemoryAlloc((POINTER*)&array,
                                (top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
		break;
	    case 3:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
		FT_VectorMemoryAlloc((POINTER*)&array,
                        (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
		break;
	    }
	    array_size = 1;
	    for (i = 0; i < dim; ++i)
	    	array_size *= (top_gmax[i] + 1);
	}
}	/* end set_solver_domain */

void ELLIPTIC_SOLVER::solve(double *soln)
{
	switch (dim)
	{
	case 1:
	    return solve1d(soln);
	case 2:
	    return solve2d(soln);
	case 3:
	    return solve3d(soln);
	}
}	/* end solve */

void ELLIPTIC_SOLVER::solve1d(double *soln)
{
	int index,index_nb[2],size;
	double k0,k_nb[2];
	double rhs,coeff[2];
	int I,I_nb[2];
	int i,l,icoords[MAXD];
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[2] = {WEST,EAST};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
        int status;
	POINTER intfc_state;

	PETSc solver;
	solver.Create(ilower, iupper-1, 3, 3);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = iupper - ilower;
	max_soln = -HUGE;
	min_soln = HUGE;

        for (i = imin; i <= imax; i++)
	{
	    index  = d_index1d(i,top_gmax);
	    comp = top_comp[index];
	    I = i_to_I[i];
	    if (I == -1) continue;

	    index_nb[0] = d_index1d(i-1,top_gmax);
	    index_nb[1] = d_index1d(i+1,top_gmax);
	    I_nb[0] = i_to_I[i-1];
	    I_nb[1] = i_to_I[i+1];
	    icoords[0] = i;
	
	    k0 = D[index];
	    num_nb = 0;
	    for (l = 0; l < 2; ++l)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
		if (status != CONST_V_PDE_BOUNDARY)
		    num_nb++;
                if (status == CONST_V_PDE_BOUNDARY ||
		    status == CONST_P_PDE_BOUNDARY)
		    index_nb[l] = index;
		k_nb[l] = 0.5*(k0 + D[index_nb[l]]);
	    	coeff[l] = k_nb[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 2; ++l)
	    {
		if (num_nb == 0) break;
		status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
                {
                    solver.Set_A(I,I_nb[l],coeff[l]);
                    aII += -coeff[l];
                }
                else if (status == CONST_P_PDE_BOUNDARY)
                {
                    aII += -coeff[l];
		    rhs += -coeff[l]*getStateVar(intfc_state);
		    use_neumann_solver = NO;
                }
	    }
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
	
		(void) printf("WARNING: isolated value!\n");
                solver.Set_A(I,I,1.0);
		rhs = soln[index];
            }
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    (void) printf("\nUsing Neumann Solver!\n");
	    if (size < 4)
	    {
	    	(void) printf("Isolated small region for solve1d()\n");
		stop_clock("Petsc Solver");
	    }
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	else
	{
	    (void) printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if(rel_residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",rel_residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetFinalRelativeResidualNorm(&rel_residual);
	    }

	}
	stop_clock("Petsc Solver");

	double *x;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
        for (i = imin; i <= imax; i++)
	{
	    index = d_index1d(i,top_gmax);
	    I = i_to_I[i];
	    if (I == -1) continue;
	    soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) max_soln = soln[index];
	    if (min_soln > soln[index]) min_soln = soln[index];
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if(debugging("step_size"))
	{
	    printf("\nThe max solution value is %.16g\n",max_soln);
	    printf("\nThe min solution value is %.16g\n",min_soln);
	}

	FT_FreeThese(1,x);
}	/* end solve1d */

void ELLIPTIC_SOLVER::solve2d(double *soln)
{
	int index,index_nb[4];
	double k0,k_nb[4];
	double rhs,coeff[4];
	int I,I_nb[4],I_oppnb[4];
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

    static PETSc solver;
	static double *x;
	static bool first = true;

	int size = iupper - ilower;
    if (first)
    {
	    solver.Create(ilower, iupper-1, 5, 5);
	    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        first = false;
    }
    else
    {
        solver.SetPrevSolnInitialGuess();
    }
	
    if (debugging("check_div"))
            printf("Enterng solve2d()\n");
	
    solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	max_soln = -HUGE;
	min_soln = HUGE;

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
	    I_oppnb[0] = ij_to_I[i+1][j];
	    I_oppnb[1] = ij_to_I[i-1][j];
	    I_oppnb[2] = ij_to_I[i][j+1];
	    I_oppnb[3] = ij_to_I[i][j-1];
	    icoords[0] = i;
	    icoords[1] = j;
	
	    k0 = D[index];
	    num_nb = 0;
	    for (l = 0; l < 4; ++l)
	    {
            status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                    &intfc_state,&hs,crx_coords);

            if (status != CONST_V_PDE_BOUNDARY)
                num_nb++;
                    
            if (status == CONST_V_PDE_BOUNDARY || status == CONST_P_PDE_BOUNDARY)
            {
                //TODO: for OUTLET use regular index_nb[l] from above??
                index_nb[l] = index;
            }

            k_nb[l] = 0.5*(k0 + D[index_nb[l]]);
            coeff[l] = k_nb[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 4; ++l)
        {
            refl_side[l] = NO;
            if (num_nb == 0) break;

            status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                    &intfc_state,&hs,crx_coords);
            
            if (status == NO_PDE_BOUNDARY)
            {
                solver.Set_A(I,I_nb[l],coeff[l]);
                aII += -coeff[l];
            }
            else if (status == CONST_P_PDE_BOUNDARY)
            {
                //TODO: should this include the inlet as well?
                //      inlet has constant pressure, but
                //      ifluid_find_state_at_crossing() returns
                //      CONST_V_PDE_BOUNDARY for it.
                
                //TODO: getStateVar() will return phi, not pressure.
                //      The pressure is updated at the flow through boundary,
                //      but phi is not.
                //
                //TODO: for OUTLET use regular index_nb[l] from above??
                rhs += -coeff[l]*getStateVar(intfc_state);
                aII += -coeff[l];
                use_neumann_solver = NO;
            }
            else if (status == CONST_V_PDE_BOUNDARY &&
                      (wave_type(hs) == NEUMANN_BOUNDARY ||
                       wave_type(hs) == MOVABLE_BODY_BOUNDARY))
            {
                //dp/dn = 0 (reflecting boundary for pressure)
                int icoords_ghost[MAXD];
                for (int m = 0; m < dim; ++m)
                    icoords_ghost[m] = icoords[m];
                
                int idir = l/2; int nb = l%2;
                icoords_ghost[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                
                double coords_ghost[MAXD];
                for (int m = 0; m < dim; ++m)
                    coords_ghost[m] = top_L[m] + icoords_ghost[m]*top_h[m];


                //TODO: solve() (this function) should be a method of the incompressible fluid solver
                //      so we can use functions such as the one below, and avoid hard coding
                //      variable types into the solver library.
                    //getRectangleCenter(index_nb[l],coords_ghost);
                    
                double intfc_crx_coords[MAXD];
                double coords_reflect[MAXD];
                double nor[MAXD];
                
                FT_ReflectPointThroughBdry(front,hs,coords_ghost,
                        comp,intfc_crx_coords,coords_reflect,nor);

                /*
                 
                double nor[MAXD];
                FT_NormalAtGridCrossing(front,icoords,
                        dir[l],comp,nor,&hs,crx_coords);
                        
                //Reflect the ghost point through intfc-mirror at crossing.
                //first reflect across the grid line containing intfc crossing,
                //which is just the coords at the current index.
                double coords_reflect[MAXD];
                for (int m = 0; m < dim; ++m)
                    coords_reflect[m] = top_L[m] + top_h[m]*icoords[m];

                //Reflect the displacement vector across the line
                //containing the intfc normal vector
                double v[MAXD];
                double vn = 0.0;

                for (int m = 0; m < dim; ++m)
                {
                    v[m] =  coords_reflect[m] - crx_coords[m];
                    vn += v[m]*nor[m];
                }

                for (int m = 0; m < dim; ++m)
                    v[m] = 2.0*vn*nor[m] - v[m];

                //The desired reflected point
                for (int m = 0; m < dim; ++m)
                    coords_reflect[m] = crx_coords[m] + v[m];

                */

                //Interpolate the pressure at the reflected point,
                //which will serve as the ghost point pressure.
                double pres_reflect;
                FT_IntrpStateVarAtCoords(front,comp,
                        coords_reflect,soln,getStateVar,
                        &pres_reflect,&soln[index]);
                //TODO: getStateVar() returns phi which is what we are solving for.
                //      More correct method would place the weights of the points
                //      used to interpolate at the reflected point into the matrix.
                //
                //      try using FrontGetRectCellIntrpCoeffs() to modify matrix

                aII += -coeff[l];
                rhs -= coeff[l]*pres_reflect; 
                use_neumann_solver = NO;

          
                /*  
                if (porosity != 0.0)
                {
                    solver.Set_A(I,I_nb[l],porosity*coeff[l]);
                    aII += -porosity*coeff[l];
                }
                else
                {
                    refl_side[l] = YES;
                }
                */
            }

        }

        /*
	    for (l = 0; l < 4; ++l)
	    {
            break;
            if (refl_side[l] == YES)
            {
                double alpha = 1.0;
                        solver.Set_A(I,I_oppnb[l],(1.0-alpha)*coeff[l]);
                        aII -= -alpha*coeff[l];
            }
	    }
        */

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
            rhs = soln[index];//TODO: or use pressure (instead of phi)?
        }

        solver.Set_b(I,rhs);
	}

    use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

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
	    if (size < 20)
	    {
	    	if (debugging("linear_solver"))
	    	    (void) printf("Isolated small region for solve2d()\n");
		stop_clock("Petsc Solver");
		return;
	    }
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetResidualNorm(&residual);
	    if(residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",residual);
		solver.Reset_x();
		solver.Solve_withPureNeumann_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetResidualNorm(&residual);
	    }

	}
	else
	{
	    if (debugging("linear_solver"))
	        (void) printf("\nUsing non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetResidualNorm(&residual);

	    if(residual > 1)
	    {
		(void) printf("\n The solution diverges! The residual "
		       "is %g. Solve again using GMRES!\n",residual);
		solver.Reset_x();
		solver.Solve_GMRES();
		solver.GetNumIterations(&num_iter);
		solver.GetResidualNorm(&residual);
	    }

	}
	stop_clock("Petsc Solver");

	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, residual = %g \n", 
			num_iter, residual);
	
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];
	    if (I == -1) continue;
	    else soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) 
	    {
		icrds_max[0] = i;
		icrds_max[1] = j;
		max_soln = soln[index];
	    }
	    if (min_soln > soln[index]) 
	    {
		icrds_min[0] = i;
		icrds_min[1] = j;
		min_soln = soln[index];
	    }
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
	    (void) printf("Max solution = %20.14f occuring at: %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1]);
            checkSolver(icrds_max,YES);
            (void) printf("Min solution = %20.14f occuring at: %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1]);
            checkSolver(icrds_min,YES);
	}
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
            (void) printf("In elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d)\n",icrds_max[0],icrds_max[1]);
            error = checkSolver(icrds_max,YES);
	}
	if (debugging("check_div"))
            printf("Leaving solve2d()\n");

	//FT_FreeThese(1,x);
}	/* end solve2d */

void ELLIPTIC_SOLVER::printIsolatedCells()
{
        int i,j,k,I,index,num_cells;

        num_cells = 0;
        switch (dim)
        {
        case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
	        I = ijk_to_I[i][j][k];
	        if (I == -1) continue;
                index = d_index3d(i,j,k,top_gmax);
                printf("Isolated cell: (%d,%d,%d) soln = %f\n",
                                i,j,k,soln[index]);
                num_cells++;
            }
            break;
        }
}       /* end printIsolatedCells */

/*
void ELLIPTIC_SOLVER::solve3d(double *soln)
{
	const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };
	
	POINTER intfc_state;
    HYPER_SURF *hs;

	COMPONENT comp;
	int index,index_nb;
    int I,I_nb;
	int crx_status;
    boolean status;

    int icoords[MAXD], icnb[MAXD];
	double crx_coords[MAXD];
	double nor[MAXD];

    //TODO: Use last time step soln as initial guess for x.
	PETSc solver;
	solver.Create(ilower, iupper-1, 7, 7);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
    boolean use_neumann_solver = YES;

	for (int k = kmin; k <= kmax; k++)
	for (int j = jmin; j <= jmax; j++)
    for (int i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;

        double aII = 0.0;
        double RHS = source[index];

        //double D_index = D[index];
        double rho_index = rho[index];
            
        int skip_count = 0;

        for (int idir = 0; idir < 3; ++idir)
        {
            for (int m = 0; m < 3; ++m)
                icnb[m] = icoords[m];

            double lambda = 1.0/sqr(top_h[idir]);

            for (int nb = 0; nb < 2; ++nb)
            {
                icnb[idir] = (nb == 0) ?
                    icoords[idir] - 1 : icoords[idir] + 1;

                index_nb  = d_index(icnb,top_gmax,3);
                I_nb = ijk_to_I[icnb[0]][icnb[1]][icnb[2]];

                double coeff_nb = 0.0;
                //double D_halfidx = 0.5*D_index;
                double rho_halfidx = 0.5*rho_index;

                crx_status = (*findStateAtCrossing)(front,icoords,
                        dir[idir][nb],comp,&intfc_state,&hs,crx_coords);

                if (crx_status)
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //INFLOW/OUTFLOW BOUNDARY ONLY, NOT NOSLIP WALL BOUNDARY!
                        if (crx_status == CONST_P_PDE_BOUNDARY)
                        {
                             //coeff_nb = 1.0*lambda/rho_index;
                            //D_halfidx += 0.5*D[index_nb];
                            //coeff_nb = lambda*D_halfidx;
                            rho_halfidx += 0.5*rho[index_nb];
                            coeff_nb = lambda/rho_halfidx;
                            aII -= coeff_nb;
                            RHS -= coeff_nb*getStateVar(intfc_state);
                            use_neumann_solver = NO;
                        }
                        else
                        {
                            //skip for constant velocity inlet boundary since
                            //we require the Bdry(grad_phi) = 0.
                            skip_count++;
                        }
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                             wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                    {
                        //dp/dn = 0 (reflecting boundary for pressure)
                        status = FT_NormalAtGridCrossing(front,icoords,
                                dir[idir][nb],comp,nor,&hs,crx_coords);

                        //Reflect the ghost point through intfc-mirror at crossing.
                        //first reflect across the grid line containing intfc crossing,
                        //which is just the coords at the current index.
                        double coords_reflect[MAXD];
                        for (int m = 0; m < 3; ++m)
                            coords_reflect[m] = top_L[m] + top_h[m]*icoords[m];

                        //Reflect the displacement vector across the line
                        //containing the intfc normal vector
                        double v[MAXD];
                        double vn = 0.0;

                        for (int m = 0; m < 3; ++m)
                        {
                            v[m] =  coords_reflect[m] - crx_coords[m];
                            vn += v[m]*nor[m];
                        }

                        for (int m = 0; m < 3; ++m)
                            v[m] = 2.0*vn*nor[m] - v[m];

                        //The desired reflected point
                        for (int m = 0; m < 3; ++m)
                            coords_reflect[m] = crx_coords[m] + v[m];

                        //Interpolate the pressure at the reflected point,
                        //which will serve as the ghost point pressure.
                        double pres_reflect;
                        FT_IntrpStateVarAtCoords(front,comp,
                                coords_reflect,soln,getStateVar,
                                &pres_reflect,&soln[index]);

                        //coeff_nb = 1.0*lambda/rho_index;
                        rho_halfidx += 0.5*rho_index;
                        coeff_nb = lambda/rho_halfidx;
                        aII -= coeff_nb;
                        RHS -= coeff_nb*pres_reflect;
                    }
                }
                else
                {
                    //NO_PDE_BOUNDARY
                    //coeff_nb = 1.0*lambda/rho_index;
                    rho_halfidx += 0.5*rho[index_nb];
                    coeff_nb = lambda/rho_halfidx;
                    solver.Set_A(I,I_nb,coeff_nb);
                    aII -= coeff_nb;
                }
            }
        }

        if (skip_count == 6)
        {
            printf("\n\nELLIPTIC_SOLVER::solve3d(): WARNING, \
                    skip_count = %d\n\n",skip_count);
            solver.Set_A(I,I,1.0);
            RHS = soln[index];
        }
        else
        {
            solver.Set_A(I,I,aII);
        }

        solver.Set_b(I,RHS);
    }
	   
        
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	use_neumann_solver = pp_min_status(use_neumann_solver);
    bool Try_GMRES = false;

	PetscInt num_iter;
    double rel_residual;

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    if (skip_neumann_solver)
	    {
	        if (debugging("PETSc"))
            {
                printf("Skip isolated small region for solve3d()\n");
                printIsolatedCells();
            }
            stop_clock("Petsc Solver");
            return;
            //TODO: go on to gmres??
	    }

        printf("\nELLIPTIC_SOLVER: Using Neumann Solver!\n");
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);
	    if (rel_residual > 1)
	    {
		    printf("\n The solution diverges! The residual "
                    "is %g. Solve again using GMRES!\n",rel_residual);
            Try_GMRES = true;
	    }
	}
	else
	{
	    printf("\nELLIPTIC_SOLVER: Using non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetFinalRelativeResidualNorm(&rel_residual);

	    if (rel_residual > 1)
	    {
            printf("\n The solution diverges! The residual "
                    "is %g. Solve again using GMRES!\n",rel_residual);
            Try_GMRES = true;
	    }
	}

    if (Try_GMRES)
    {
        solver.Reset_x();
        solver.Solve_GMRES();
        solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);
	    
        if (rel_residual > 1)
	    {
            printf("\n The solution diverges using GMRES. \
                    The residual is %g. Exiting ...\n",rel_residual);
            clean_up(EXIT_FAILURE);
        }
    }

	stop_clock("Petsc Solver");

	double *x;
	int size = iupper - ilower;
	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
    {
        printf("In poisson_solver(): \
                num_iter = %d, rel_residual = %g\n", 
                num_iter,rel_residual);
    }

	for (int k = kmin; k <= kmax; k++)
	for (int j = jmin; j <= jmax; j++)
    for (int i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;
	    
        soln[index] = x[I-ilower];

//        //
//	    if (max_soln < soln[index]) 
//	    {
//            max_soln = soln[index];
//            icrds_max[0] = i;
//            icrds_max[1] = j;
//            icrds_max[2] = k;
//	    }
//	    if (min_soln > soln[index]) 
//	    {
//            min_soln = soln[index];
//            icrds_min[0] = i;
//            icrds_min[1] = j;
//            icrds_min[2] = k;
//	    }
//        //
	}
//	//
//    pp_global_max(&max_soln,1);
//	pp_global_min(&min_soln,1);
//
//	if (debugging("step_size"))
//	{
//	    printf("Max solution = %20.14f occuring at: %d %d %d\n",
//                        max_soln,icrds_max[0],icrds_max[1],icrds_max[2]);
//            checkSolver(icrds_max,YES);
//            printf("Min solution = %20.14f occuring at: %d %d %d\n",
//                        min_soln,icrds_min[0],icrds_min[1],icrds_min[2]);
//            checkSolver(icrds_min,YES);
//	}
//	if (debugging("elliptic_error"))
//        {
//            double error,max_error = 0.0;
//            for (k = kmin; k <= kmax; k++)
//            for (j = jmin; j <= jmax; j++)
//            for (i = imin; i <= imax; i++)
//            {
//                icoords[0] = i;
//                icoords[1] = j;
//                icoords[2] = k;
//                if (ijk_to_I[i][j][k] == -1) continue;
//                error = checkSolver(icoords,NO);
//                if (error > max_error)
//                {
//                    max_error = error;
//                    icrds_max[0] = i;
//                    icrds_max[1] = j;
//                    icrds_max[2] = k;
//                }
//            }
//            (void) printf("In elliptic solver:\n");
//            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
//            (void) printf("Occuring at (%d %d %d)\n",icrds_max[0],
//                                icrds_max[1],icrds_max[2]);
//            error = checkSolver(icrds_max,YES);
//        }
//    //
	FT_FreeThese(1,x);
}*/	/* end solve3d */

void ELLIPTIC_SOLVER::solve3d(double *soln)
{
	int index,index_nb[6];
	double k0,k_nb[6];
	double rhs,coeff[6];
	int I,I_nb[6];
	int i,j,k,l,icoords[MAXD],icrds_max[MAXD],icrds_min[MAXD];
	COMPONENT comp;
	double aII;
	int num_nb;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	int status;
	POINTER intfc_state;

    static PETSc solver;
	static double *x;
	static bool first = true;
	
    int size = iupper - ilower;
    
    if (first)
    {
	    solver.Create(ilower, iupper-1, 7, 7);
	    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
        first = false;
    }
    else
    {
        solver.SetPrevSolnInitialGuess();
    }
	
    solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	max_soln = -HUGE;
	min_soln = HUGE;

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;

	    index_nb[0] = d_index3d(i-1,j,k,top_gmax);
	    index_nb[1] = d_index3d(i+1,j,k,top_gmax);
	    index_nb[2] = d_index3d(i,j-1,k,top_gmax);
	    index_nb[3] = d_index3d(i,j+1,k,top_gmax);
	    index_nb[4] = d_index3d(i,j,k-1,top_gmax);
	    index_nb[5] = d_index3d(i,j,k+1,top_gmax);
	    I_nb[0] = ijk_to_I[i-1][j][k];
	    I_nb[1] = ijk_to_I[i+1][j][k];
	    I_nb[2] = ijk_to_I[i][j-1][k];
	    I_nb[3] = ijk_to_I[i][j+1][k];
	    I_nb[4] = ijk_to_I[i][j][k-1];
	    I_nb[5] = ijk_to_I[i][j][k+1];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	
	    k0 = D[index];
	    num_nb = 0;
	    for (l = 0; l < 6; ++l)
	    {
            status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                    &intfc_state,&hs,crx_coords);
    
            if (status != CONST_V_PDE_BOUNDARY)
                num_nb++;

            if (status == CONST_V_PDE_BOUNDARY ||
                    status == CONST_P_PDE_BOUNDARY)
            {
                //TODO: for OUTLET use regular index_nb[l] from above??
                index_nb[l] = index;
            }
    
            k_nb[l] = 0.5*(k0 + D[index_nb[l]]);
	    	coeff[l] = k_nb[l]/(top_h[l/2]*top_h[l/2]); 
	    }

	    rhs = source[index];

	    aII = 0.0;
	    for (l = 0; l < 6; ++l)
	    {
		    if (num_nb == 0) break;
		
            status = (*findStateAtCrossing)(front,icoords,dir[l],comp,
                                &intfc_state,&hs,crx_coords);
            
            if (status == NO_PDE_BOUNDARY)
            {
                solver.Set_A(I,I_nb[l],coeff[l]);
                aII += -coeff[l];
            }
            else if (status == CONST_P_PDE_BOUNDARY)
            {
                //TODO: should this include the inlet as well?
                //      inlet has constant pressure, but
                //      ifluid_find_state_at_crossing() returns
                //      CONST_V_PDE_BOUNDARY for it.
                
                //TODO: getStateVar() will return phi, not pressure.
                //      The pressure is updated at the flow through boundary,
                //      but phi is not.
                
                //TODO: for OUTLET use regular index_nb[l] from above??
                aII += -coeff[l];
                rhs += -coeff[l]*getStateVar(intfc_state);
                use_neumann_solver = NO;
            }
            else if (status == CONST_V_PDE_BOUNDARY &&
                      (wave_type(hs) == NEUMANN_BOUNDARY ||
                       wave_type(hs) == MOVABLE_BODY_BOUNDARY))
            {
                //dp/dn = 0 (reflecting boundary for pressure)
                int icoords_ghost[MAXD];
                for (int m = 0; m < dim; ++m)
                    icoords_ghost[m] = icoords[m];
                
                int idir = l/2; int nb = l%2;
                icoords_ghost[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                
                double coords_ghost[MAXD];
                for (int m = 0; m < dim; ++m)
                    coords_ghost[m] = top_L[m] + icoords_ghost[m]*top_h[m];


                //TODO: solve() (this function) should be a method of the incompressible fluid solver
                //      so we can use functions such as the one below, and avoid hard coding
                //      variable types into the solver library.
                    //getRectangleCenter(index_nb[l],coords_ghost);
                    
                double intfc_crx_coords[MAXD];
                double coords_reflect[MAXD];
                double nor[MAXD];
                
                FT_ReflectPointThroughBdry(front,hs,coords_ghost,
                        comp,intfc_crx_coords,coords_reflect,nor);

                /*
                
                double nor[MAXD];
                FT_NormalAtGridCrossing(front,icoords,
                        dir[l],comp,nor,&hs,crx_coords);
                        
                //Reflect the ghost point through intfc-mirror at crossing.
                //first reflect across the grid line containing intfc crossing,
                //which is just the coords at the current index.
                double coords_reflect[MAXD];
                for (int m = 0; m < dim; ++m)
                    coords_reflect[m] = top_L[m] + top_h[m]*icoords[m];

                //Reflect the displacement vector across the line
                //containing the intfc normal vector
                double v[MAXD];
                double vn = 0.0;

                for (int m = 0; m < dim; ++m)
                {
                    v[m] =  coords_reflect[m] - crx_coords[m];
                    vn += v[m]*nor[m];
                }

                for (int m = 0; m < dim; ++m)
                    v[m] = 2.0*vn*nor[m] - v[m];

                //The desired reflected point
                for (int m = 0; m < dim; ++m)
                    coords_reflect[m] = crx_coords[m] + v[m];
                
                */

                //Interpolate the pressure at the reflected point,
                //which will serve as the ghost point pressure.
                double pres_reflect;
                FT_IntrpStateVarAtCoords(front,comp,coords_reflect,soln,
                        getStateVar,&pres_reflect,&soln[index]);
                //TODO: getStateVar() returns phi which is what we are solving for.
                //      More correct method would place the weights of the points
                //      used to interpolate at the reflected point into the matrix.
                //
                //      try using FrontGetRectCellIntrpCoeffs() to modify matrix

                aII += -coeff[l];
                rhs -= coeff[l]*pres_reflect; 
                use_neumann_solver = NO;
            }
        }

	    ///
	    //  This change reflects the need to treat point with only one
	    //  interior neighbor (a convex point). Not sure why PETSc cannot
	    //  handle such case. If we have better understanding, this should
	    //  be changed back.
	    // /
	    if(num_nb > 0)
	    {
            solver.Set_A(I,I,aII);
	    }
        else
        {
            (void) printf("WARNING: isolated value!\n");
            solver.Set_A(I,I,1.0);
            rhs = soln[index];
        }

        solver.Set_b(I,rhs);
	}

	solver.SetMaxIter(40000);
	solver.SetTol(1.0e-10);

	use_neumann_solver = pp_min_status(use_neumann_solver);
    bool Try_GMRES = false;

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    if (skip_neumann_solver)
	    {
	        if (debugging("PETSc"))
            {
            (void) printf("Skip isolated small region for solve3d()\n");
                printIsolatedCells();
            }
            stop_clock("Petsc Solver");
            return;
            //TODO: go on to gmres??
	    }
	    printf("\nELLIPTIC_SOLVER: Using Neumann Solver!\n");
	    solver.Solve_withPureNeumann();
	    solver.GetNumIterations(&num_iter);
	    solver.GetResidualNorm(&residual);
	    if(residual > 1)
	    {
		    printf("\n The solution diverges! The residual "
                    "is %g after %d iterations. Solve again using GMRES!\n",
                    residual,num_iter);
            //clean_up(ERROR);
            //TODO: go on to gmres??
            Try_GMRES = true;
	    }

	}
	else
	{
	    printf("\nELLIPTIC_SOLVER: Using non-Neumann Solver!\n");
	    solver.Solve();
	    solver.GetNumIterations(&num_iter);
	    solver.GetResidualNorm(&residual);

	    if(residual > 1)
	    {
		    printf("\n The solution diverges! The residual "
                    "is %g after %d iterations. Solve again using GMRES!\n",
                    residual,num_iter);
            ///
            //solver.Reset_x();
            //solver.Solve_GMRES();
            //solver.GetNumIterations(&num_iter);
            //solver.GetResidualNorm(&residual);
            ///
            Try_GMRES = true;
	    }
	}

    if (Try_GMRES)
    {
        solver.Reset_x();
        solver.Solve_GMRES();
        solver.GetNumIterations(&num_iter);
        solver.GetResidualNorm(&residual);
	    
        if(residual > 1)
	    {
            printf("\n The solution diverges using GMRES. \
                    The residual is %g after %d iterations. Exiting ...\n",
                    residual,num_iter);
            clean_up(EXIT_FAILURE);
        }
    }

	stop_clock("Petsc Solver");

	solver.Get_x(x);

	if (debugging("PETSc"))
    {
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, residual = %g \n", 
			num_iter, residual);
    }

	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = ijk_to_I[i][j][k];
	    if (I == -1) continue;
	    soln[index] = x[I-ilower];
	    if (max_soln < soln[index]) 
	    {
		max_soln = soln[index];
		icrds_max[0] = i;
		icrds_max[1] = j;
		icrds_max[2] = k;
	    }
	    if (min_soln > soln[index]) 
	    {
		min_soln = soln[index];
		icrds_min[0] = i;
		icrds_min[1] = j;
		icrds_min[2] = k;
	    }
	}
	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
	    printf("Max solution = %20.14f occuring at: %d %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1],icrds_max[2]);
            checkSolver(icrds_max,YES);
            printf("Min solution = %20.14f occuring at: %d %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1],icrds_min[2]);
            checkSolver(icrds_min,YES);
	}
	if (debugging("elliptic_error"))
        {
            double error,max_error = 0.0;
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                if (ijk_to_I[i][j][k] == -1) continue;
                error = checkSolver(icoords,NO);
                if (error > max_error)
                {
                    max_error = error;
                    icrds_max[0] = i;
                    icrds_max[1] = j;
                    icrds_max[2] = k;
                }
            }
            (void) printf("In elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d %d)\n",icrds_max[0],
                                icrds_max[1],icrds_max[2]);
            error = checkSolver(icrds_max,YES);
        }

    //NOTE: Was causing a double free() call somehow.
    //      Made x static and only allocate on first call
    //      to circumvent it.
	
    //FT_FreeThese(1,x);
}   /* end solve3d */

double ELLIPTIC_SOLVER::checkSolver(
	int *icoords,
	boolean print_details)
{
	int i,l,m;
	int comp;
	double w[2];
	int id0,index_nb;
	double dw[2],coefs[2],lhs,rhs;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
        double crx_coords[MAXD];
        int status;
        POINTER intfc_state;
	int icnb[MAXD];
	double denom = 0.0;

	if (print_details)
	{
	    (void) printf("\nEntering checkSolver()\n");
	    (void) printf("icoords = ");
	    for (i = 0; i < dim; ++i)
	    	(void) printf("%d ",icoords[i]);
	    (void) printf("\n");
	}

	id0 = d_index(icoords,top_gmax,dim);
	comp = top_comp[id0];
	lhs = 0;
	for (l = 0; l < dim; ++l)
	{
	    if (print_details)
	    	(void) printf("Direction %d:\n",l);
	    for (i = 0; i < dim; ++i)
		icnb[i] = icoords[i];
	    for (m = 0; m < 2; ++m)
	    {
		status = (*findStateAtCrossing)(front,icoords,dir[l][m],comp,
                                &intfc_state,&hs,crx_coords);
		icnb[l] = (m == 0) ? icoords[l] - 1 : icoords[l] + 1;
		index_nb = d_index(icnb,top_gmax,dim);
		switch(status)
		{
		case NO_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d NO_PDE_BOUNDARY\n",m);
		    coefs[m] = 0.5*(D[id0] + D[index_nb]);
		    w[0] = soln[id0];
		    w[1] = soln[index_nb];
		    break;
		case CONST_V_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d CONST_V_PDE_BOUNDARY\n",m);
		    coefs[m] = D[id0];
		    w[0] = soln[id0];
		    w[1] = soln[id0];
		    break;
		case CONST_P_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d CONST_P_PDE_BOUNDARY\n",m);
		    coefs[m] = D[id0];
		    w[0] = soln[id0];
		    w[1] = getStateVar(intfc_state);
		    break;
		case NEUMANN_PDE_BOUNDARY:
	    	    if (print_details)
		    	(void) printf("Side %d NEUMANN_PDE_BOUNDARY\n",m);
		    w[0] = soln[id0];
		    w[1] = soln[id0];
		    break;
		default:
		    (void) printf("Side %d Unknown BOUNDARY\n",m);
		    clean_up(ERROR);
		}
		dw[m] = (w[1] - w[0])/top_h[l];
		if (denom < fabs(coefs[m]*dw[m]/2.0/top_h[l]))
		    denom = fabs(coefs[m]*dw[m]/2.0/top_h[l]);
	    }
	    if (print_details)
	    {
	    	(void) printf("Coefs: %f %f\n",coefs[0],coefs[1]);
	    	(void) printf("C*dw: %f %f\n",coefs[0]*dw[0],coefs[1]*dw[1]);
	    }
	    lhs += (coefs[1]*dw[1] + coefs[0]*dw[0])/top_h[l];
	}
	rhs = source[id0];
	if (print_details)
	{
	    (void) printf("Solution = %20.14f\n",soln[id0]);
	    (void) printf("LHS = %20.14f  RHS = %20.14f\n",lhs,rhs);
	    (void) printf("LHS - RHS = %20.14f\n",lhs-rhs);
	    if (denom != 0.0)
	    	(void) printf("Relative error = %20.14g\n",fabs(lhs-rhs)/denom);
	    else
	    	(void) printf("Relative error = %20.14g\n",0.0);
	    (void) printf("Leaving checkSolver()\n\n");
	}
	return fabs(lhs-rhs)/denom;
}	/* end checkSolver */
