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

DOUBLE_ELLIPTIC_SOLVER::DOUBLE_ELLIPTIC_SOLVER(Front &front)
    : front(&front)
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

            //TODO: is this off by an index?
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

void DOUBLE_ELLIPTIC_SOLVER::dsolve(double *soln)
{
    switch (dim)
    {
        case 2:
            return dsolve2d(soln);
        case 3:
            return dsolve3d(soln);
    }
}       /* end dsolve */

void DOUBLE_ELLIPTIC_SOLVER::dsolve2d(double *soln)
{
	int index,index_nb,num_nb,size,ic;
	int I,I_nb;
	int i,j,l,icoords[MAXD],icn[MAXD],icnb[MAXD],iknb[MAXD];
	int icrds_max[MAXD],icrds_min[MAXD];
	int idir,nb;
	int status;
	double rhs;
    double k_nb,coeff_nb;
	double crx_coords[MAXD];
	double aII;
	double rel_residual = 0.0;
	double h2[MAXD];
	double *x;

	COMPONENT comp;
	GRID_DIRECTION dir[2][2] = {{WEST,EAST},{SOUTH,NORTH}};
	PetscInt num_iter = 0;
	HYPER_SURF *hs;
	POINTER intfc_state;

	PETSc solver;

	if (debugging("check_div"))
            printf("Entering dsolve2d()\n");
        printf("eilower = %d  eiupper = %d\n",eilower,eiupper);

	solver.Create(eilower, eiupper-1, 5, 5);

	size = eiupper - eilower;

	for (i = 0; i < dim; ++i)
	    h2[i] = 4.0*sqr(top_h[i]);

    /*
    ///////////////////////////////////////////////
        printf("J RANGE\n\n");
        printf("\next_gmax[1] = %d\n",ext_gmax[1]);
        for( int j = 0; j <= ext_gmax[1]; ++j )
        {
            printf("j = %d",j);
            
            if( j == 0 || j == ext_gmax[1] )
                printf("\t\tpast boundary");
            if( j == ext_imin[1] )
                printf("\t\text_imin[1]");
            if( j == ext_imin[1]+ext_l[1] )
                printf("\t\text_imin[1]+ext_l[1]");
            if( j == imin[1] )
                printf("\t\timin[1]");
            if( j == imax[1] )
                printf("\t\timax[1]");
            if( j == ext_imax[1]-ext_u[1] )
                printf("\t\text_imax[1]-ext_u[1]");
            if( j == ext_imax[1] )
                printf("\t\text_imax[1]");
            printf("\n");
        }
        printf("\n");
    ///////////////////////////////////////////////
    */

    /* Set original layer coefficients */
	for (j = imin[1]; j <= imax[1]; j++)
    {
        for (i = imin[0]; i <= imax[0]; i++)
	    {
            icoords[0] = i;
            icoords[1] = j;

            index  = d_index2d(i,j,top_gmax);
            comp = top_comp[index];

            I = dij_to_I[i+ext_l[0]][j+ext_l[1]];
            if (I == -1) continue;

            //k0 = 1.0;
            //double k0 = D[index];
            
            rhs = 1.0;
            //rhs = source[index];

            if (i == 5)
                printf("j = %d rhs = %g\n",j+ext_l[1],rhs);


            aII = 0.0;

            for (idir = 0; idir < dim; ++idir)
            {
                for (nb = 0; nb < 2; ++nb)
                {
                    icnb[0] = icoords[0];   iknb[0] = icoords[0];
                    icnb[1] = icoords[1];   iknb[1] = icoords[1];

                    icnb[idir] = (nb == 0) ? icoords[idir]-2 : icoords[idir]+2;
                    iknb[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;

                    status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                                    comp,&intfc_state,&hs,crx_coords);

                    if (status == CONST_V_PDE_BOUNDARY)
                    {
                        if (wave_type(hs) == NEUMANN_BOUNDARY)
                        {
                            icnb[idir] = (nb == 0) ? icoords[idir]+1 : icoords[idir]-1;
                            iknb[idir] = icoords[idir];
                        }
                    }
                    else if (status == NO_PDE_BOUNDARY)
                    {
                        icn[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;

                        status = (*findStateAtCrossing)(front,icn,dir[idir][nb],
                                comp,&intfc_state,&hs,crx_coords);

                        if (status == CONST_V_PDE_BOUNDARY)
                        {
                            if (wave_type(hs) == NEUMANN_BOUNDARY)
                            {
                                icnb[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;
                            }
                        }
                    }

                    iknb[0] += ext_l[0];
                    iknb[1] += ext_l[1];
                    index_nb = d_index(iknb,ext_gmax,dim);

                    //k_nb = 1.0;
                    k_nb = ext_D[index_nb];
                    //coeff_nb = 1.0/k_nb/h2[idir];
                    coeff_nb = k_nb/h2[idir];

                    //Set neighbor
                    I_nb = dij_to_I[icnb[0]+ext_l[0]][icnb[1]+ext_l[1]];
                    solver.Set_A(I,I_nb,coeff_nb);

                    aII += -coeff_nb;

                    //NOTE: Using sparse laplacian discretization,
                    //      no cancellation occurs when enforcing
                    //      Neumann boundary.
                    //      i.e. Always set neighbor, and add -coeff_nb to aII.
                    
                }
            }
            
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
	    }
    }

    /* Set extended layer coefficients */
    boolean extended_domain = NO;
	for (j = ext_imin[1]; j <= ext_imax[1]; j++)
    {
        for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	    {
            if (i >= ext_imin[0]+ext_l[0] && i <= ext_imax[0]-ext_u[0] &&
                j >= ext_imin[1]+ext_l[1] && j <= ext_imax[1]-ext_u[1])
                continue;

            extended_domain = YES;
            index = d_index2d(i,j,ext_gmax);
            
            icoords[0] = i;
            icoords[1] = j;
            
            I = dij_to_I[i][j];
            if (I == -1) exit(0);
            
            //double k0 = 1.0;
            double k0 = ext_D[index];
            //rhs = 0.0;
            rhs = ext_source[index];

            if (i == 5)
                printf("j = %d rhs = %g\n",j,rhs);
            //printf("(%2d %2d) I: %d  I_nb: ",i,j,I);
            

            aII = 0.0;
            for (idir = 0; idir < dim; ++idir)
            {
                for (nb = 0; nb < 2; ++nb)
                {
                    icnb[0] = i;    iknb[0] = i;
                    icnb[1] = j;    iknb[1] = j;

                    /* Single spacing discretizatin */
                    icnb[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;
                    iknb[idir] = (nb == 0) ? icoords[idir]-1 : icoords[idir]+1;

                    //k_nb = 1.0;
                    //index_nb = d_index(iknb,ext_gmax,dim);
                    //k_nb = ext_D[index_nb];
                    
                    coeff_nb = k0/sqr(top_h[idir]);
                    //coeff_nb = 1.0/sqr(top_h[idir]);

                    //Need rho itself to evaluate at half index,
                    //      can't just adverage D = 1/rho.

                    //coeff_nb = 1.0/k_half/sqr(top_h[idir]);
                    //double k_half = 0.5*(k0 + k_nb);

                    if ((ext_l[idir] != 0 && icnb[idir] < ext_imin[idir]) || 
                        (ext_u[idir] != 0 && icnb[idir] > ext_imax[idir]))
                    {
                        int bdryface = 2*idir + nb;
                        if (bdryface == ConstantBdryPosition)
                        {
                            //Dirichlet Boundary (inlet) equal 0
                            double dirichletVal = 0.0;
                            rhs -= coeff_nb*dirichletVal;
                            aII += -coeff_nb;
                        }
                        else
                        {
                            //Neumann boundary (outlet) equal 0
                            //rhs -= 0.0;
                            continue;
                        }
                    }
                    else
                    {
                        I_nb = dij_to_I[icnb[0]][icnb[1]];
                        solver.Set_A(I,I_nb,coeff_nb);
                        aII += -coeff_nb;
                    }

                    //coeff_nb = k0/sqr(top_h[idir]);
                    //coeff_nb = ext_D[index]/sqr(top_h[idir]);
                    
                    //solver.Set_A(I,I_nb,coeff_nb);
                    //aII += -coeff_nb;
                    //printf("%d ",I_nb);
                }
            }
            
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
        }
    }

	//use_neumann_solver = pp_min_status(use_neumann_solver);
	
	solver.SetMaxIter(40000);
	solver.SetTol(1.0e-10);

	start_clock("Petsc Solver");
        solver.Solve_PetscDecide();
	stop_clock("Petsc Solver");

	FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));
	solver.Get_x(x);

	max_soln = -HUGE;
	min_soln = HUGE;

        /* Move to solver */
	if (debugging("PETSc"))
	    (void) printf("In poisson_solver(): "
	       		"num_iter = %d, rel_residual = %g \n",num_iter,rel_residual);
	
    
        /* Extracting ext_array from Petsc solver */
    
    //Need to bring the periodic buffer with ext_array (extended solution vector),
    //in order to check the the error in dcheckSolver().
    int* lbuf = front->rect_grid->lbuf;
    int* ubuf = front->rect_grid->ubuf;

    for (j = ext_imin[1] - lbuf[1]; j <= ext_imax[1] + ubuf[1]; j++)
        for (i = ext_imin[0] - lbuf[0]; i <= ext_imax[0] + ubuf[0]; i++)
	{
	    index = d_index2d(i,j,ext_gmax);
	    I = dij_to_I[i][j];
	    if (I == -1) continue;
	    ext_array[index] = x[I-eilower];
	}
    

        /* Extracting solution from ext_array */

    for (j = imin[1]; j <= imax[1]; j++)
        for (i = imin[0]; i <= imax[0]; i++)
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
            if (i == (imin[0] + imax[0])/2+1)
            {
                printf("soln[%d][%d] = %f\n",i,j,soln[ic]);
            }
        }

    if (extended_domain)
    {
        FILE *xfile = fopen("test1.xg","w");
	    for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        {
            i = (ext_imin[0] + ext_imax[0])/2;
	        I = dij_to_I[i][j];
            fprintf(xfile,"%f %f\n",j*top_h[1],x[I-eilower]);
        }
        fclose(xfile);
    }   

	if (debugging("elliptic_error"))
    {
    //This is prototype for dcheckSolver()
            
        int imid = (ext_imin[0]+ext_imax[0])/2; 
        printf("line: i = %d\n",imid);
        for (j = ext_imin[1]; j <= ext_imax[1]; j++)
        {
            for (i = ext_imin[0]; i <= ext_imax[0]; i++)
            {
                double LHS = 0.0;
                double RHS = 0.0;

                index = d_index2d(i,j,ext_gmax);

                if (j > ext_imin[1]+1 && j < ext_imax[1]-1 && i == imid)
                {
                    if( j >= ext_imin[1]+ext_l[1] &&
                            j <= ext_imax[1]-ext_u[1] )
                    {
                        int index_nb1 = d_index2d(i+2,j,ext_gmax);
                        int index_nb2 = d_index2d(i-2,j,ext_gmax);
                        int index_nb3 = d_index2d(i,j+2,ext_gmax);
                        int index_nb4 = d_index2d(i,j-2,ext_gmax);
                        
                        int index_knb1 = d_index2d(i+1,j,ext_gmax);
                        int index_knb2 = d_index2d(i-1,j,ext_gmax);
                        int index_knb3 = d_index2d(i,j+1,ext_gmax);
                        int index_knb4 = d_index2d(i,j-1,ext_gmax);

                        LHS = ( (ext_array[index_nb1] - ext_array[index])*ext_D[index_knb1]
                            +  (ext_array[index_nb2] - ext_array[index])*ext_D[index_knb2] )/h2[0]
                            +  ( (ext_array[index_nb3] - ext_array[index])*ext_D[index_knb3]
                            +  (ext_array[index_nb4] - ext_array[index])*ext_D[index_knb4] )/h2[1];
                        

                        /*
                        LHS = (ext_array[index_nb1] + ext_array[index_nb2]
                            + ext_array[index_nb3] + ext_array[index_nb4] 
                            - 4.0*ext_array[index])/4.0/top_h[1]/top_h[1];
                        */

                        RHS = 1.0;
                        //RHS = ext_source[index];
                    }
                    else
                    {
                        int index_nb1 = d_index2d(i+1,j,ext_gmax);
                        int index_nb2 = d_index2d(i-1,j,ext_gmax);
                        int index_nb3 = d_index2d(i,j+1,ext_gmax);
                        int index_nb4 = d_index2d(i,j-1,ext_gmax);

                        //TODO: Evaluate rho at 1/2 index in the single scheme?
                        //      If so, need rho itself not D = 1/rho
                        /*
                        double k0 = ext_D[index];
                        double knb1 = 0.5*(k0 + ext_D[index_nb1]);
                        double knb2 = 0.5*(k0 + ext_D[index_nb2]);
                        double knb3 = 0.5*(k0 + ext_D[index_nb3]);
                        double knb4 = 0.5*(k0 + ext_D[index_nb4]);

                        LHS = (ext_array[index_nb1] - ext_array[index])/knb1
                            +  (ext_array[index_nb2] - ext_array[index])/knb2
                            +  (ext_array[index_nb3] - ext_array[index])/knb3
                            +  (ext_array[index_nb4] - ext_array[index])/knb4;
                        */
                        
                        //double k0 = 1.0;
                        double k0 = ext_D[index];

                        LHS = ( (ext_array[index_nb1] - ext_array[index])
                            +  (ext_array[index_nb2] - ext_array[index]) )*k0/sqr(top_h[0])
                            +  ( (ext_array[index_nb3] - ext_array[index])
                            +  (ext_array[index_nb4] - ext_array[index]) )*k0/sqr(top_h[1]);


                        /*
                        LHS = (ext_array[index_nb1] + ext_array[index_nb2]
                            + ext_array[index_nb3] + ext_array[index_nb4] 
                            - 4.0*ext_array[index])/top_h[1]/top_h[1];
                        */

                        RHS = 0.0;
                        //RHS = ext_source[index];
                    }

                    double residual = LHS - RHS;

                    printf("j = %d  LHS = %g  RHS = %g  res = %8.6g\n",
                                    j,LHS,RHS,residual);
                }
            //End prototype for dchecksolve()


            /*
            if (j == ext_imin[1] && i == 10)
            {
                index = d_index2d(i,j,top_gmax);
                int index_nb1 = d_index2d(i+2,j,ext_gmax);
                int index_nb2 = d_index2d(i-2,j,ext_gmax);
                int index_nb3 = d_index2d(i,j+2,ext_gmax);
                int index_nb4 = d_index2d(i,j,top_gmax);
                double LHS = D[index]*(soln[index_nb1] + soln[index_nb2]
                        + soln[index_nb3] + soln[index_nb4] 
                        - 4.0*soln[index])/4.0/top_h[1]/top_h[1];
                double RHS = source[index];
                printf("j = %d  LHS = %f  RHS = %f  res = %5.2g\n",
                                j,LHS,RHS,LHS-RHS);
            } 
            if (j == ext_imin[1]+1 && i == 10)
            {
                index = d_index2d(i,j,top_gmax);
                int index_nb1 = d_index2d(i+2,j,ext_gmax);
                int index_nb2 = d_index2d(i-2,j,ext_gmax);
                int index_nb3 = d_index2d(i,j+2,ext_gmax);
                int index_nb4 = d_index2d(i,j-1,ext_gmax);
                double LHS = D[index]*(soln[index_nb1] + soln[index_nb2]
                        + soln[index_nb3] + soln[index_nb4] 
                        - 4.0*soln[index])/h2[1];
                double RHS = source[index];
                printf("j = %d  LHS = %f  RHS = %f  res = %5.2g\n",
                                j,LHS,RHS,LHS-RHS);
            } 
            if (j == ext_imax[1] && i == 10)
            {
                index = d_index2d(i,j,top_gmax);
                int index_nb1 = d_index2d(i+2,j,ext_gmax);
                int index_nb2 = d_index2d(i-2,j,ext_gmax);
                int index_nb3 = d_index2d(i,j+2,ext_gmax);
                int index_nb4 = d_index2d(i,j-2,ext_gmax);
                double LHS = D[index]*(soln[index_nb1] + soln[index_nb2]
                        + 1.0 + soln[index_nb4] 
                        - 4.0*soln[index])/h2[1];
                double RHS = source[index];
                printf("j = %d  LHS = %f  RHS = %f  res = %5.2g\n",
                                j,LHS,RHS,LHS-RHS);
            } 
            if (j == ext_imax[1]-1 && i == 10)
            {
                index = d_index2d(i,j,top_gmax);
                int index_nb1 = d_index2d(i+2,j,ext_gmax);
                int index_nb2 = d_index2d(i-2,j,ext_gmax);
                int index_nb3 = d_index2d(i,j+2,ext_gmax);
                int index_nb4 = d_index2d(i,j-2,ext_gmax);
                double LHS = D[index]*(soln[index_nb1] + soln[index_nb2]
                        + 1.0 + soln[index_nb4] 
                        - 4.0*soln[index])/4.0/top_h[1]/top_h[1];
                double RHS = source[index];
                printf("j = %d  LHS = %f  RHS = %f  res = %5.2g\n",
                                j,LHS,RHS,LHS-RHS);
            } 
            */
        
            }
        }

    }


	if (debugging("elliptic_error"))
    {
    //This is prototype for dcheckSolver()
                    
        int jline = 6;
        //int jmid = (ext_imin[1]+ext_imax[1])/2; 
        printf("\n line: j = %d\n",jline);
	    for (i = ext_imin[0]; i <= ext_imax[0]; i++)
        {
            for (j = ext_imin[1]; j <= ext_imax[1]; j++)
            {
                double LHS = 0.0;
                double RHS = 0.0;

                index = d_index2d(i,j,ext_gmax);

                //if (i > ext_imin[0]+1 && i < ext_imax[0]-1 && j == jmid)
                if (i >= ext_imin[0] && i <= ext_imax[0] && j == jline)
                {
                    if( i >= ext_imin[0]+ext_l[0] &&
                            i <= ext_imax[0]-ext_u[0] )
                    {
                        int index_nb1 = d_index2d(i+2,j,ext_gmax);
                        int index_nb2 = d_index2d(i-2,j,ext_gmax);
                        int index_nb3 = d_index2d(i,j+2,ext_gmax);
                        int index_nb4 = d_index2d(i,j-2,ext_gmax);
                        
                        int index_knb1 = d_index2d(i+1,j,ext_gmax);
                        int index_knb2 = d_index2d(i-1,j,ext_gmax);
                        int index_knb3 = d_index2d(i,j+1,ext_gmax);
                        int index_knb4 = d_index2d(i,j-1,ext_gmax);

                        LHS = ( (ext_array[index_nb1] - ext_array[index])*ext_D[index_knb1]
                            +  (ext_array[index_nb2] - ext_array[index])*ext_D[index_knb2] )/h2[0]
                            +  ( (ext_array[index_nb3] - ext_array[index])*ext_D[index_knb3]
                            +  (ext_array[index_nb4] - ext_array[index])*ext_D[index_knb4] )/h2[1];
                        

                        /*
                        LHS = (ext_array[index_nb1] + ext_array[index_nb2]
                            + ext_array[index_nb3] + ext_array[index_nb4] 
                            - 4.0*ext_array[index])/4.0/top_h[1]/top_h[1];
                        */

                        RHS = 1.0;
                        //RHS = ext_source[index];
                    }
                    else
                    {
                        int index_nb1 = d_index2d(i+1,j,ext_gmax);
                        int index_nb2 = d_index2d(i-1,j,ext_gmax);
                        int index_nb3 = d_index2d(i,j+1,ext_gmax);
                        int index_nb4 = d_index2d(i,j-1,ext_gmax);

                        //TODO: Evaluate rho at 1/2 index in the single scheme?
                        //      If so, need rho itself not D = 1/rho
                        /*
                        double k0 = ext_D[index];
                        double knb1 = 0.5*(k0 + ext_D[index_nb1]);
                        double knb2 = 0.5*(k0 + ext_D[index_nb2]);
                        double knb3 = 0.5*(k0 + ext_D[index_nb3]);
                        double knb4 = 0.5*(k0 + ext_D[index_nb4]);

                        LHS = (ext_array[index_nb1] - ext_array[index])/knb1
                            +  (ext_array[index_nb2] - ext_array[index])/knb2
                            +  (ext_array[index_nb3] - ext_array[index])/knb3
                            +  (ext_array[index_nb4] - ext_array[index])/knb4;
                        */
                        
                        //double k0 = 1.0;
                        double k0 = ext_D[index];
                        
                        LHS = ( (ext_array[index_nb1] - ext_array[index])
                            +  (ext_array[index_nb2] - ext_array[index]) )*k0/sqr(top_h[0])
                            +  ( (ext_array[index_nb3] - ext_array[index])
                            +  (ext_array[index_nb4] - ext_array[index]) )*k0/sqr(top_h[1]);


                        /*
                        LHS = (ext_array[index_nb1] + ext_array[index_nb2]
                            + ext_array[index_nb3] + ext_array[index_nb4] 
                            - 4.0*ext_array[index])/top_h[1]/top_h[1];
                        */

                        RHS = 0.0;
                        //RHS = ext_source[index];
                    }

                    double residual = LHS - RHS;

                    printf("i = %d  LHS = %g  RHS = %g  res = %8.6g\n",
                                    i,LHS,RHS,residual);
                }
            //End prototype for dchecksolve()
            }
        }
    }


    //////////////////
    if (extended_domain)
    {
        printf("line 452 in double_ellip.cpp...exiting\n");
        clean_up(0);
    }
    //////////////////
    
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
            (void) printf("In dual elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d)\n",icrds_max[0],icrds_max[1]);
            error = dcheckSolver(icrds_max,YES);
        }
	if (debugging("check_div"))
            printf("Leaving dsolve2d()\n");

	FT_FreeThese(1,x);
}	/* end dsolve2d */

//TODO: Implement dsolve3d when dsolve2d complete
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
            (void) printf("In dual elliptic solver:\n");
            (void) printf("Max relative elliptic error: %20.14f\n",max_error);
            (void) printf("Occuring at (%d %d %d)\n",icrds_max[0],
                                icrds_max[1],icrds_max[2]);
            error = dcheckSolver(icrds_max,YES);
        }

	FT_FreeThese(1,x);
}	/* end dsolve3d */

boolean DOUBLE_ELLIPTIC_SOLVER::icoordsInterior(int *icoords)
{
    for (int i = 0; i < dim; ++i)
    {
        if (ext_l[i] != 0 && icoords[i] < ext_l[i])
            return NO;
        if (ext_u[i] != 0 && icoords[i] > ext_imax[i]-ext_u[i])
            return NO;
    }
    return YES;
}

/*
double DOUBLE_ELLIPTIC_SOLVER::dcheckSolver(
        int *icoords,
        boolean print_details)
{
    if( icoordsInterior(icoords) )
        return dcheckSolverInterior(icoords,print_details);
    else
        return dcheckSolverExtended(icoords,print_details);
}
*/

//TODO: Implement dcheckSolver
double DOUBLE_ELLIPTIC_SOLVER::dcheckSolver(
	int *icoords,
	boolean print_details)
{
	int i,j,l,m;
	int comp;
	double w[2];
	int id0,index_nb;
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
	    (void) printf("\nEntering dcheckSolver()\n");
	    (void) printf("icoords = ");
	    for (i = 0; i < dim; ++i)
	    	(void) printf("%d ",icoords[i]);
	    (void) printf("\n");
	}

	id0 = d_index(icoords,ext_gmax,dim);
	comp = dtop_comp[id0];
	lhs = 0.0;

	for (l = 0; l < dim; ++l)
	{
	    if (print_details)
	    	printf("Direction %d:\n",l);

	    for (i = 0; i < dim; ++i)
        {
            icnb[i] = icoords[i];
            iknb[i] = icoords[i];
        }
	    
        for (m = 0; m < 2; ++m)
	    {
            iknb[l] = (m == 0) ? icoords[l]-1 : icoords[l]+1;
            icnb[l] = (m == 0) ? icoords[l]-2 : icoords[l]+2;
		
            status = (*findStateAtCrossing)(front,icoords,dir[l][m],comp,
                                &intfc_state,&hs,crx_coords);

            //TODO: keep working from here
            if (status == CONST_V_PDE_BOUNDARY)
            {
                if (wave_type(hs) == NEUMANN_BOUNDARY)
                {
                    if (print_details)
                        (void) printf("Side %d-0 CONST_V_PDE_BOUNDARY\n",m);

                    icnb[l] = (m == 0) ? icoords[l] + 1 : icoords[l] - 1;
                    index_nb = d_index(icnb,ext_gmax,dim);

                    coefs[m] = ext_D[id0];
                    w[0] = ext_array[id0];
                    w[1] = ext_array[index_nb];
                }

            }

////////////////////////////
		if (status == NO_PDE_BOUNDARY)
		{
		    icnb[l] = (m == 0) ? icoords[l] - 1 : icoords[l] + 1;
		    index_nb = d_index(icnb,top_gmax,dim);
		    status = (*findStateAtCrossing)(front,icnb,dir[l][m],comp,
                                &intfc_state,&hs,crx_coords);
		
            if (status == NO_PDE_BOUNDARY)
		    {
		    	if (print_details)
		    	    (void) printf("Side %d NO_PDE_BOUNDARY\n",m);
		    	
                coefs[m] = D[index_nb];
		    	icnb[l] = (m == 0) ? icoords[l] - 2 : icoords[l] + 2;
		    	index_nb = d_index(icnb,top_gmax,dim);
			
                w[0] = soln[id0];
			    w[1] = soln[index_nb];
		    }
		    else if (status == CONST_V_PDE_BOUNDARY)
		    {
		    	if (print_details)
		    	    (void) printf("Side %d-1 CONST_V_PDE_BOUNDARY\n",m);
		    	coefs[m] = 0.5*(D[id0] + D[index_nb]);
			w[0] = soln[id0];
			w[1] = soln[index_nb];
		    }
		    else if (status == CONST_P_PDE_BOUNDARY)
		    {
		    	if (print_details)
		    	    (void) printf("Side %d-1 CONST_P_PDE_BOUNDARY\n",m);
		    	coefs[m] = 0.5*(D[id0] + D[index_nb]);
			w[0] = soln[id0];
		    	w[1] = getStateVar(intfc_state);
		    }
		}
		else if (status == CONST_V_PDE_BOUNDARY)
		{
		    if (print_details)
		    	(void) printf("Side %d-0 CONST_V_PDE_BOUNDARY\n",m);
		    coefs[m] = D[id0];
		    w[0] = soln[id0];
		    w[1] = soln[id0];
		}
		else if (status == CONST_P_PDE_BOUNDARY)
		{
		    if (print_details)
		    	(void) printf("Side %d-0 CONST_P_PDE_BOUNDARY\n",m);
		    icnb[l] = (m == 0) ? icoords[l] + 1 : icoords[l] - 1;
		    index_nb = d_index(icnb,top_gmax,dim);
		    coefs[m] = D[id0];
		    w[0] = soln[index_nb];
		    w[1] = getStateVar(intfc_state);
		}
		dw[m] = (w[1] - w[0])/2.0/top_h[l];
		if (denom < fabs(coefs[m]*dw[m]/2.0/top_h[l]))
		    denom = fabs(coefs[m]*dw[m]/2.0/top_h[l]);
	    
        }
	    if (print_details)
            {
	    	(void) printf("Coefs: %f %f\n",coefs[0],coefs[1]);
	    	(void) printf("C*dw: %f %f\n",coefs[0]*dw[0],coefs[1]*dw[1]);
	    }
	
        lhs += (coefs[1]*dw[1] + coefs[0]*dw[0])/2.0/top_h[l];
	
    }

	rhs = source[id0];
	if (print_details)
        {
	    (void) printf("Solution = %20.14f\n",soln[id0]);
	    (void) printf("LHS = %20.14f  RHS = %20.14f\n",lhs,rhs);
	    (void) printf("LHS - RHS = %20.14f\n",lhs-rhs);
	    (void) printf("Relative error = %20.14g\n",fabs(lhs-rhs)/denom);
	    (void) printf("Leaving dcheckSolver()\n\n");
	}
	return fabs(lhs-rhs)/denom;
}	/* end dcheckSolver */

