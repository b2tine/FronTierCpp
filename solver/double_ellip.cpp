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

DOUBLE_ELLIPTIC_SOLVER::DOUBLE_ELLIPTIC_SOLVER(Front &front):front(&front)
{
        porosity = 0.0;
}

DOUBLE_ELLIPTIC_SOLVER::~DOUBLE_ELLIPTIC_SOLVER()
{
        FT_FreeThese(3,ext_source,ext_D,ext_array);
}

void DOUBLE_ELLIPTIC_SOLVER::set_solver_domain(void)
{

	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
	int i;

	dim = Dimension(front->interf);
        top_comp = T->components;
        top_gmax = rgr->gmax;
	top_h = rgr->h;
	array_size = 1;
	for (i = 0; i < dim; ++i) array_size *= (ext_gmax[i] + 1);
	FT_VectorMemoryAlloc((POINTER*)&ext_source,array_size,FLOAT);
	FT_VectorMemoryAlloc((POINTER*)&ext_D,array_size,FLOAT);
	FT_VectorMemoryAlloc((POINTER*)&ext_array,array_size,FLOAT);

	switch(dim)
	{
	case 2:
	    imin[0] = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imin[1] = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax[0] = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    imax[1] = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    break;
	case 3:
	    imin[0] = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imin[1] = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imin[2] = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax[0] = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    imax[1] = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    imax[2] = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    break;
	}

}	/* end set_solver_domain */

void DOUBLE_ELLIPTIC_SOLVER::dsolve(double *soln)
{
        switch (dim)
        {
        case 2:
            return dsolve2d(soln);
        case 3:
            return dsolve3d(soln);
        }
}       /* end solve */

void DOUBLE_ELLIPTIC_SOLVER::set_extension()
{
        int i,j,k,l,ic,id;
        int icoords[MAXD];
        double copy_D;
        static boolean first = YES;

        printf("Enter set_extension()\n");
        /*
        printf("Check H-source:\n");
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            ic = d_index2d(i,5,top_gmax);
            printf("source[%d][%d] = %f\n",i,5,source[ic]);
        }
        printf("Check V-source:\n");
        for (j = 0; j <= top_gmax[1]; ++j)
        {
            ic = d_index2d(5,j,top_gmax);
            printf("source[%d][%d] = %f\n",5,j,source[ic]);
        }
        */
        printf(" imin = %d %d\n",imin[0],imin[1]);
        printf(" imax = %d %d\n",imax[0],imax[1]);
        printf("ext_imin = %d %d\n",ext_imin[0],ext_imin[1]);
        printf("ext_imax = %d %d\n",ext_imax[0],ext_imax[1]);
        printf("ext_l = %d %d\n",ext_l[0],ext_l[1]);
        printf("ext_u = %d %d\n",ext_u[0],ext_u[1]);
        switch(dim)
        {
        case 2:
            for (i = 0; i <= ext_gmax[0]; ++i)
            for (j = 0; j <= ext_gmax[1]; ++j)
            {
                ic = d_index2d(i,j,ext_gmax);
                ext_D[ic] = ext_source[ic] = 0.0;
            }
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
                ic = d_index2d(i,j,top_gmax);
                id = d_index2d(i+ext_l[0],j+ext_l[1],ext_gmax);
                ext_D[id] = D[ic];
                ext_source[id] = source[ic];
            }
            /*
            printf("After assignment V-dir:\n");
            for (j = 0; j <= ext_gmax[1]; ++j)
            {
                ic = d_index2d(7,j,ext_gmax);
                printf("(%2d,%2d): ext_D = %f  ext_source = %f\n",
                            7,j,ext_D[ic],ext_source[ic]);
            }
            printf("After assignment H-dir:\n");
            for (i = 0; i <= ext_gmax[0]; ++i)
            {
                ic = d_index2d(i,7,ext_gmax);
                printf("(%2d,%2d): ext_D = %f  ext_source = %f\n",
                            i,7,ext_D[ic],ext_source[ic]);
            }
            */
            for (j = ext_imin[1]; j <= ext_imax[1]; ++j)
            {
                icoords[1] = j;
                if (ext_l[0] != 0)
                {
                    icoords[0] = imin[0] + ext_l[0];
                    id = d_index(icoords,ext_gmax,2);
                    copy_D = ext_D[id];
                    for (i = 0; i <= ext_l[0]; ++i)
                    {
                        icoords[0] = i;
                        id = d_index(icoords,ext_gmax,2);
                        ext_D[id] = copy_D;
                    }
                }
                if (ext_u[0] != 0)
                {
                    icoords[0] = ext_gmax[0] - ext_u[0] - 1;
                    id = d_index(icoords,ext_gmax,2);
                    copy_D = ext_D[id];
                    for (i = 0; i <= ext_u[0]; ++i)
                    {
                        icoords[0] = ext_gmax[0] - i;
                        id = d_index(icoords,ext_gmax,2);
                        ext_D[id] = copy_D;
                    }
                }
            }
            for (i = ext_imin[0]; i <= ext_imax[0]; ++i)
            {
                icoords[0] = i;
                if (ext_l[1] != 0)
                {
                    icoords[1] = imin[1] + ext_l[1];
                    id = d_index(icoords,ext_gmax,2);
                    copy_D = ext_D[id];
                    for (j = 0; j <= ext_l[1]; ++j)
                    {
                        icoords[1] = j;
                        id = d_index(icoords,ext_gmax,2);
                        ext_D[id] = copy_D;
                    }
                }
                if (ext_u[1] != 0)
                {
                    icoords[1] = ext_gmax[1] - ext_u[1] - 1;
                    id = d_index(icoords,ext_gmax,2);
                    copy_D = ext_D[id];
                    for (j = 0; j <= ext_u[1]; ++j)
                    {
                        icoords[1] = ext_gmax[1] - j;
                        id = d_index(icoords,ext_gmax,2);
                        ext_D[id] = copy_D;
                    }
                }
            }
            /*
            printf("After extension V-dir:\n");
            for (j = 0; j <= ext_gmax[1]; ++j)
            {
                ic = d_index2d(7,j,ext_gmax);
                printf("(%2d,%2d): ext_D = %f\n",7,j,ext_D[ic]);
            }
            printf("After extension H-dir:\n");
            for (i = 0; i <= ext_gmax[0]; ++i)
            {
                ic = d_index2d(i,7,ext_gmax);
                printf("(%2d,%2d): ext_D = %f\n",i,7,ext_D[ic]);
            }
            */
            break;
        case 3:
            break;
        }
        
}       /* end set_extension */

void DOUBLE_ELLIPTIC_SOLVER::dsolve2d(double *soln)
{
	int index,index_nb,num_nb,size,ic;
	int I,I_nb;
	int i,j,l,icoords[MAXD],icn[MAXD],icnb[MAXD],iknb[MAXD];
	int icrds_max[MAXD],icrds_min[MAXD];
	int idir,nb;
	int status;
	double rhs;
        double k0,k_nb,coeff_nb;
	double crx_coords[MAXD];
	double aII;
	double rel_residual = 0.0;
	double h2[MAXD];
	double *x;
	COMPONENT comp;
	GRID_DIRECTION dir[2][2] = {{WEST,EAST},{SOUTH,NORTH}};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	HYPER_SURF *hs;
	POINTER intfc_state;

	PETSc solver;

	if (debugging("check_div"))
            printf("Entering dsolve2d()\n");
        printf("eilower = %d  eiupper = %d\n",eilower,eiupper);

	solver.Create(eilower, eiupper-1, 5, 5);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();

	size = eiupper - eilower;
	max_soln = -HUGE;
	min_soln = HUGE;

	for (i = 0; i < dim; ++i)
	    h2[i] = 4.0*sqr(top_h[i]);

        double soln_set;
        /* Set original layer coefficients */
        max_rhs = 0.0;
	for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
	{
            index  = d_index2d(i,j,top_gmax);
            comp = top_comp[index];
            I = dij_to_I[i+ext_l[0]][j+ext_l[1]];
            icoords[0] = i;
            icoords[1] = j;

            if (I == -1) continue;

            rhs = source[index];
            if (max_rhs < fabs(rhs)) max_rhs = fabs(rhs);

            aII = 0.0;
            for (idir = 0; idir < dim; ++idir)
            {
                for (nb = 0; nb < 2; ++nb)
                {
                    icnb[0] = iknb[0] = icoords[0];
                    icnb[1] = iknb[1] = icoords[1];
                    icnb[idir] = (nb == 0) ? icoords[idir]-2 : icoords[idir]+2;
                    iknb[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;

                    status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                                    comp,&intfc_state,&hs,crx_coords);
                    if (status == CONST_V_PDE_BOUNDARY)
                    {
                        if (wave_type(hs) == NEUMANN_BOUNDARY)
                        {
                            icnb[idir] = (nb == 0) ? icoords[idir]+1 : 
                                    icoords[idir]-1;
                            iknb[idir] = icoords[idir];
                        }
                    }
                    else if (status == NO_PDE_BOUNDARY)
                    {
                        icn[idir] = (nb == 0) ? icoords[idir]-1 :
                                    icoords[idir]+1;
                        status = (*findStateAtCrossing)(front,icn,
                                        dir[idir][nb],comp,&intfc_state,&hs,
                                        crx_coords);
                        if (status == CONST_V_PDE_BOUNDARY)
                        {
                            if (wave_type(hs) == NEUMANN_BOUNDARY)
                            {
                                icnb[idir] = (nb == 0) ? icoords[idir]-1 : 
                                        icoords[idir]+1;
                            }
                        }
                    }

                    I_nb = dij_to_I[icnb[0]+ext_l[0]][icnb[1]+ext_l[1]];
                    iknb[0] += ext_l[0];
                    iknb[1] += ext_l[1];
                    index_nb = d_index(iknb,ext_gmax,dim);
                    k_nb = ext_D[index_nb];
                    coeff_nb = k_nb/h2[idir];

                    /* Set neighbor at boundary */

                    solver.Set_A(I,I_nb,coeff_nb);
                    aII += -coeff_nb;
                }
            }
            
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
	}

        /* Set extended layer coefficients */
        boolean extended_domain = NO;
	for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	{
            if (i >= ext_imin[0]+ext_l[0] && i <= ext_imax[0]-ext_u[0] &&
                j >= ext_imin[1]+ext_l[1] && j <= ext_imax[1]-ext_u[1])
                continue;
            extended_domain = YES;
            index  = d_index2d(i,j,ext_gmax);
            icn[0] = icoords[0] = i;
            icn[1] = icoords[1] = j;
            I = dij_to_I[i][j];
            boolean fixed_point = NO;
            for (idir = 0; idir < dim; ++idir)
            {
                if (icoords[idir] == ext_imin[idir] && ext_l[idir] != 0)
                {
                    ic = d_index(icoords,ext_gmax,dim);
                    comp = ext_comp[ic];
                    status = (*findStateAtCrossing)(front,icn,dir[idir][0],
                                    comp,&intfc_state,&hs,crx_coords);
                    if (status == CONST_V_PDE_BOUNDARY &&
                        wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        soln_set = getStateVar(intfc_state);
                        solver.Set_A(I,I,1.0);
                        solver.Set_b(I,soln_set);
                        fixed_point = YES;
                    }
                }
                else if (icoords[idir] == ext_imax[idir] && ext_u[idir] != 0)
                {
                    ic = d_index(icoords,ext_gmax,dim);
                    comp = ext_comp[ic];
                    icn[idir] -= (ext_l[idir] + ext_u[idir]);
                    status = (*findStateAtCrossing)(front,icn,dir[idir][1],
                                    comp,&intfc_state,&hs,crx_coords);
                    if (status == CONST_V_PDE_BOUNDARY &&
                        wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        soln_set = getStateVar(intfc_state);
                        solver.Set_A(I,I,1.0);
                        solver.Set_b(I,soln_set);
                        fixed_point = YES;
                    }
                }
            }
            if (fixed_point == YES) continue;

            k0 = ext_D[index];
            aII = 0.0;
            rhs = ext_source[index];
            for (idir = 0; idir < dim; ++idir)
            {
                for (nb = 0; nb < 2; ++nb)
                {
                    icn[0] = icoords[0];
                    icn[1] = icoords[1];
                    icn[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;
                    if ((ext_l[idir] != 0 && icn[idir] < ext_imin[idir]) || 
                        (ext_u[idir] != 0 && icn[idir] > ext_imax[idir]))
                        continue;       // treat as Neumann boundary

                    /* Single spacing discretizatin */
                    I_nb = dij_to_I[icn[0]][icn[1]];
                    coeff_nb = ext_D[index]/sqr(top_h[idir]);
                    solver.Set_A(I,I_nb,coeff_nb);
                    aII += -coeff_nb;
                    if (ext_u[idir] != 0 && icn[idir] > ext_imax[idir])
                        rhs -= coeff_nb*soln_set;
                }
            }
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
        }

	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
        solver.Solve_PetscDecide();
	stop_clock("Petsc Solver");

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

        /* Move to solver */
	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
        /* Extracting ext_array from Petsc solver */
	for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	{
	    index = d_index2d(i,j,ext_gmax);
	    I = dij_to_I[i][j];
	    if (I == -1) continue;
	    ext_array[index] = x[I-eilower];
	}
        FT_ParallelExchExtendedGridArrayBuffer(ext_array,front,ext_gmax,NULL);

        /* Extracting solution from ext_array */
	for (j = imin[1]-1; j <= imax[1]+1; j++)
        for (i = imin[0]-1; i <= imax[0]+1; i++)
	{
            ic = d_index2d(i,j,top_gmax); 
            index = d_index2d(i+ext_l[0],j+ext_l[1],ext_gmax);
            soln[ic] = ext_array[index];
	    if (max_soln < soln[ic]) 
	    {
		max_soln = soln[ic];
		icrds_max[0] = i;
		icrds_max[1] = j;
	    }
	    if (min_soln > soln[ic]) 
	    {
		min_soln = soln[ic];
		icrds_min[0] = i;
		icrds_min[1] = j;
	    }
        }

        FILE *xfile = fopen("test1.xg","w");
        double norm = 0.0;
	for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        {
            i = (ext_imin[0] + ext_imax[0])/2;
	    I = dij_to_I[i][j];
            if (fabs(x[I-eilower]) > norm) norm = fabs(x[I-eilower]);
        }
	for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        {
            i = (ext_imin[0] + ext_imax[0])/2;
	    I = dij_to_I[i][j];
            fprintf(xfile,"%f %f\n",(double)j,x[I-eilower]/norm);
        }
        norm = 0.0;
	for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        {
            i = (ext_imin[0] + ext_imax[0])/2;
            ic = d_index2d(i,j,ext_gmax);
            if (fabs(ext_source[ic]) > norm) norm = fabs(ext_source[ic]);
        }
        fprintf(xfile,"Next\n");
        fprintf(xfile,"color=red\n");
	for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        {
            i = (ext_imin[0] + ext_imax[0])/2;
            ic = d_index2d(i,j,ext_gmax);
            fprintf(xfile,"%f %f\n",(double)j,ext_source[ic]/norm);
        }
        fclose(xfile);

	pp_global_max(&max_soln,1);
	pp_global_min(&min_soln,1);

	if (debugging("step_size"))
	{
	    (void) printf("Max solution = %20.14f occuring at: %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1]);
            dcheckSolver(icrds_max,YES);
            (void) printf("Min solution = %20.14f occuring at: %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1]);
            dcheckSolver(icrds_min,YES);
	}
	if (debugging("elliptic_error"))
        {
            double error,max_error = 0.0;
            for (j = ext_imin[1]; j <= ext_imax[1]; j++)
            for (i = ext_imin[0]; i <= ext_imax[0]; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                if (dij_to_I[i][j] == -1) continue;
                error = dcheckSolver(icoords,NO);
                if (error > max_error)
                {
                    max_error = error;
                    icrds_max[0] = i;
                    icrds_max[1] = j;
                }
            }
            (void) printf("In double elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14g\n",max_error);
            (void) printf("Occuring at (%d %d)\n",icrds_max[0],icrds_max[1]);
            error = dcheckSolver(icrds_max,YES);
        }
        //if (extended_domain) exit(0);
	if (debugging("check_div"))
            printf("Leaving dsolve2d()\n");

	FT_FreeThese(1,x);
}	/* end dsolve2d */

void DOUBLE_ELLIPTIC_SOLVER::dsolve3d(double *soln)
{
	int index,index_nb,size;
	double rhs,coeff[3][2];
	int I,I_nb;
	int i,j,k,l,icoords[MAXD],icnb[MAXD];
	int icrds_max[MAXD],icrds_min[MAXD];
	COMPONENT comp;
	double aII;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	boolean use_neumann_solver = YES;
	PetscInt num_iter = 0;
	double rel_residual = 0.0;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	int status;
	POINTER intfc_state;
	int idir,nb;
	double h2[MAXD];
	double *x;

	PETSc solver;
	solver.Create(eilower, eiupper-1, 13, 13);
	solver.Reset_A();
	solver.Reset_b();
	solver.Reset_x();
	size = eiupper - eilower;
	max_soln = -HUGE;
	min_soln = HUGE;

	for (i = 0; i < dim; ++i)
	    h2[i] = 4.0*sqr(top_h[i]);

	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    I = dijk_to_I[i][j][k];

	    index  = d_index(icoords,top_gmax,dim);
	    comp = top_comp[index];
	    if (I == -1) continue;
	
	    rhs = source[index];
	    aII = 0.0;

	    for (idir = 0; idir < dim; ++idir)
	    for (nb = 0; nb < 2; ++nb)
	    {
	    	for (l = 0; l < dim; ++l)
	    	    icnb[l] = icoords[l];
		status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
				comp,&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
		{
		    icnb[idir] = (nb == 0) ? icoords[idir] - 1 : 
					icoords[idir] + 1;
		    index_nb = d_index(icnb,top_gmax,dim);
		    status = (*findStateAtCrossing)(front,icnb,dir[idir][nb],
				    comp,&intfc_state,&hs,crx_coords);
                    if (status == NO_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = D[index_nb]/h2[idir];
			icnb[idir] = (nb == 0) ? icoords[idir] - 2:
                                        icoords[idir] + 2;
			I_nb = dijk_to_I[icnb[0]][icnb[1]][icnb[2]];
		    	solver.Set_A(I,I_nb,coeff[idir][nb]);
                    	aII += -coeff[idir][nb];
		    }
		    else if (status == CONST_P_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = 0.5*(D[index] + D[index_nb])/h2[idir];
                    	aII += -coeff[idir][nb];
			rhs += -coeff[idir][nb]*getStateVar(intfc_state);
			use_neumann_solver = NO;
		    }
		    else if (status == CONST_V_PDE_BOUNDARY)
		    {
		    	coeff[idir][nb] = 0.5*(D[index] + D[index_nb])/h2[idir];
			I_nb = dijk_to_I[icnb[0]][icnb[1]][icnb[2]];
		    	solver.Set_A(I,I_nb,coeff[idir][nb]);
                    	aII += -coeff[idir][nb];
		    }
		}
		else if (status == CONST_P_PDE_BOUNDARY)
		{
		    coeff[idir][nb] = D[index]/h2[idir];
		    icnb[idir] = (nb == 0) ? icoords[idir] + 1:
                                        icoords[idir] - 1;
		    I_nb = dijk_to_I[icnb[0]][icnb[1]][icnb[2]];
		    solver.Set_A(I,I_nb,-coeff[idir][nb]);
		    rhs += -coeff[idir][nb]*getStateVar(intfc_state);
		    use_neumann_solver = NO;
		}
	    }
	    /*
	     * This change reflects the need to treat point with only one
	     * interior neighbor (a convex point). Not sure why PETSc cannot
	     * handle such case. If we have better understanding, this should
	     * be changed back.
	     */
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
	}
	use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1e-10);

	start_clock("Petsc Solver");
	if (use_neumann_solver)
	{
	    (void) printf("\nUsing Neumann Solver!\n");
	    if (size < 8)
	    {
	    	(void) printf("Isolated small region for dsolve3d()\n");
		stop_clock("Petsc Solver");
		return;
	    }
	    (void) printf("Neumann solver not working for dsolve3d()\n");
	    clean_up(ERROR);
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

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	
	for (k = imin[2]; k <= imax[2]; k++)
	for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    I = dijk_to_I[i][j][k];
	    if (I == -1) continue;
	    soln[index] = x[I-eilower];
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
	    (void) printf("Max solution = %20.14f occuring at: %d %d %d\n",
                        max_soln,icrds_max[0],icrds_max[1],icrds_max[2]);
            dcheckSolver(icrds_max,YES);
            (void) printf("Min solution = %20.14f occuring at: %d %d %d\n",
                        min_soln,icrds_min[0],icrds_min[1],icrds_min[2]);
            dcheckSolver(icrds_min,YES);
	}
	if (debugging("elliptic_error"))
        {
            double error,max_error = 0.0;
            for (k = imin[2]; k <= imax[2]; k++)
            for (j = imin[1]; j <= imax[1]; j++)
            for (i = imin[0]; i <= imax[0]; i++)
            {
                icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
                if (dijk_to_I[i][j][k] == -1) continue;
                error = dcheckSolver(icoords,NO);
                if (error > max_error)
                {
                    max_error = error;
                    icrds_max[0] = i;
                    icrds_max[1] = j;
                    icrds_max[2] = k;
                }
            }
            (void) printf("In double elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d %d)\n",icrds_max[0],
                                icrds_max[1],icrds_max[2]);
            error = dcheckSolver(icrds_max,YES);
        }

	FT_FreeThese(1,x);
}	/* end dsolve3d */

double DOUBLE_ELLIPTIC_SOLVER::dcheckSolver(
	int *icoords,
	boolean print_details)
{
        if (icoordsInterior(icoords))
            return dcheckSolverInterior(icoords,print_details);
        else
            return dcheckSolverExtended(icoords,print_details);
}

bool DOUBLE_ELLIPTIC_SOLVER::icoordsInterior(int *icoords)
{
        for (int i = 0; i < dim; ++i)
        {
            if (ext_l[i] != 0 & ext_u[i] != 0)
            {
                if (icoords[i] < ext_imin[i] + ext_l[i] ||
                    icoords[i] > ext_imax[i] - ext_u[i])
                {
                    return false;
                }
            }
        }
        return true;
}


bool DOUBLE_ELLIPTIC_SOLVER::icoordsBoundary(int *icoords)
{
        for (int i = 0; i < dim; ++i)
        { 
            if ( ext_l[i] != 0 & ext_u[i] != 0)
            {
                if (icoords[i] <= ext_imin[i] || icoords[i] >= ext_imax[i])
                    return true;
            }
            else if (icoords[i] < ext_imin[i] || icoords[i] > ext_imax[i] )
            {
                return true;
            }
        }
        return false;
}

double DOUBLE_ELLIPTIC_SOLVER::dcheckSolverInterior(
	int *icoords,
	boolean print_details)
{
	int i,j,l;
	int comp;
	double w[2];
	int id0,index_nb,index_knb;
	double dw[2],coefs[2],lhs,rhs;
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
        double crx_coords[MAXD],w_crx;
        int status;
        POINTER intfc_state;
	int icn[MAXD];
	int icnb[MAXD];
	int iknb[MAXD];
	double denom = 0.0;

	if (print_details)
	{
	    (void) printf("\nEntering dcheckSolverInterior()\n");
	    (void) printf("icoords = ");
	    for (i = 0; i < dim; ++i)
	    	(void) printf("%d ",icoords[i]);
	    (void) printf("\n");
	}

	id0 = d_index(icoords,ext_gmax,dim);
	comp = ext_comp[id0];
	lhs = 0.0;

	for (int idir = 0; idir < dim; ++idir)
	{
	    if (print_details)
	    	printf("Direction %d:\n",idir);
	    
            for (i = 0; i < dim; ++i)
            {
	    	icnb[i] = icoords[i];
	    	iknb[i] = icoords[i];
            }

	    for (int nb = 0; nb < 2; ++nb)
	    {
                icnb[idir] = (nb == 0) ? icoords[idir]-2 : icoords[idir]+2;
                iknb[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;
                index_nb = d_index(icnb,ext_gmax,dim);
	        index_knb = d_index(iknb,ext_gmax,dim);
                coefs[nb] = ext_D[index_knb];
                status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                            comp,&intfc_state,&hs,crx_coords);
                if (status == CONST_V_PDE_BOUNDARY)
                {
                    if (wave_type(hs) == NEUMANN_BOUNDARY)
                    {
                        if (print_details)
                            (void) printf("Side %d CONST_V_PDE_BOUNDARY -- NEUMANN\n",nb);
                        icnb[idir] = (nb == 0) ? icoords[idir] + 1 : 
                                    icoords[idir] - 1;
                        iknb[idir] = icoords[idir];
                        index_knb = d_index(iknb,ext_gmax,dim);
                        coefs[nb] = ext_D[id0];
                        index_nb = d_index(icnb,ext_gmax,dim);
                    }
                }
                else if (status == NO_PDE_BOUNDARY)
                {
                    if (print_details)
                        (void) printf("Side %d NO_PDE_BOUNDARY\n",nb);
                    icn[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                    status = (*findStateAtCrossing)(front,icn,dir[idir][nb],
                                comp,&intfc_state,&hs,crx_coords);

                    if (status == CONST_V_PDE_BOUNDARY)
                    {
                        if (wave_type(hs) == NEUMANN_BOUNDARY)
                        {
                            if (print_details)
                                (void) printf("Side %d CONST_V_PDE_BOUNDARY -- NEUMANN\n",nb);
                    
                            icnb[idir] = (nb == 0) ? icoords[idir] - 1 : 
                                            icoords[idir] + 1;
                            index_nb = d_index(icnb,ext_gmax,dim);
                        }
                    }
                }
                w[0] = ext_array[id0];
                w[1] = ext_array[index_nb];
                dw[nb] = (w[1] - w[0])/2.0/top_h[idir];
            }
            if (print_details)
            {
                (void) printf("Coefs: %f %f\n",coefs[0],coefs[1]);
                (void) printf("C*dw: %f %f\n",coefs[0]*dw[0],coefs[1]*dw[1]);
            }
	    lhs += (coefs[1]*dw[1] + coefs[0]*dw[0])/2.0/top_h[idir];
        }
        rhs = ext_source[id0];
	if (print_details)
        {
	    (void) printf("Solution = %20.14g\n",ext_array[id0]);
	    (void) printf("LHS = %20.14g  RHS = %20.14f\n",lhs,rhs);
	    (void) printf("LHS - RHS = %20.14g\n",lhs-rhs);
	    (void) printf("Relative error = %20.14g\n",fabs(lhs-rhs)/max_rhs);
	    (void) printf("Leaving dcheckSolverInterior()\n\n");
	}
        if (fabs(lhs-rhs)/max_rhs > 0.001)
            printf("icoords = %d %d  rel_error = %20.14g\n",icoords[0],
                                icoords[1],fabs(lhs-rhs)/max_rhs);
        return fabs(lhs-rhs)/max_rhs;
}       /* end dcheckSolverInterior */


bool absmax(double a, double b)
{
        a = fabs(a);
        b = fabs(b);
        return a < b;
}

double DOUBLE_ELLIPTIC_SOLVER::dcheckSolverExtended(
	int *icoords,
	boolean print_details)
{
	int comp;
	int i,idir,id0,index_nb,index_knb;
	int icn[MAXD],icnb[MAXD],iknb[MAXD];
        int status;
	double w[2];
	double dw[2],coefs[2],lhs,rhs;
        double crx_coords[MAXD],w_crx,soln_set;
        bool free_end[2];
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	HYPER_SURF *hs;
        POINTER intfc_state;

	if (print_details)
	{
	    (void) printf("\nEntering dcheckSolverExtended()\n");
	    (void) printf("icoords = ");
	    for (int i = 0; i < dim; ++i)
	    	(void) printf("%d ",icoords[i]);
	    (void) printf("\n");
	}

	id0 = d_index(icoords,ext_gmax,dim);
	comp = ext_comp[id0];
	lhs = 0.0;
	rhs = ext_source[id0];

	for (idir = 0; idir < dim; ++idir)
	{
	    if (print_details)
	    	printf("Direction %d:\n",idir);
	    
            free_end[0] = free_end[1] = false;
            for (i = 0; i < dim; ++i)
            {
	    	icn[i] = icoords[i];
	    	icnb[i] = icoords[i];
            }
            if (icoords[idir] == ext_imin[idir] && ext_l[idir] != 0)
            {
                status = (*findStateAtCrossing)(front,icn,dir[idir][0],
                                comp,&intfc_state,&hs,crx_coords);
                if (status == CONST_V_PDE_BOUNDARY &&
                    wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    soln_set = getStateVar(intfc_state);
                    lhs = ext_array[id0];
                    rhs = soln_set;
                    break;
                }
                else
                    free_end[0] = true;
            }
            else if (icoords[idir] == ext_imax[idir] && ext_u[idir] != 0)
            {
                icn[idir] -= (ext_l[idir] + ext_u[idir]);
                status = (*findStateAtCrossing)(front,icn,dir[idir][1],
                                    comp,&intfc_state,&hs,crx_coords);
                if (status == CONST_V_PDE_BOUNDARY &&
                    wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    soln_set = getStateVar(intfc_state);
                    lhs = ext_array[id0];
                    rhs = soln_set;
                    break;
                }
                else
                    free_end[1] = true;
            }


	    for (int nb = 0; nb < 2; ++nb)
	    {
                if (free_end[nb] == true)
                {
                    dw[nb] = 0.0;
                    continue;
                }
                icnb[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;
                index_nb = d_index(icnb,ext_gmax,dim);
	        
                coefs[nb] = ext_D[id0];
                w[0] = ext_array[id0];
                w[1] = ext_array[index_nb];
                dw[nb] = (w[1] - w[0])/top_h[idir];
            }

            if (print_details)
            {
                (void) printf("Coefs: %f %f\n",coefs[0],coefs[1]);
                (void) printf("C*dw: %f %f\n",coefs[0]*dw[0],coefs[1]*dw[1]);
            }
	    lhs += (coefs[1]*dw[1] + coefs[0]*dw[0])/top_h[idir];
        }

	if (print_details)
        {
	    (void) printf("Solution = %20.14g\n",ext_array[id0]);
	    (void) printf("LHS = %20.14g  RHS = %20.14g\n",lhs,rhs);
	    (void) printf("LHS - RHS = %20.14g\n",lhs-rhs);
	    (void) printf("Relative error = %20.14g\n",fabs(lhs-rhs)/max_rhs);
	    (void) printf("Leaving dcheckSolverExtended()\n\n");
	}
        return fabs(lhs-rhs)/max_rhs;
}
