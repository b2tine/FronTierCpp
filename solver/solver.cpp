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

/*		PETSc.c
 *  Only for one node.
 *      This class PETSc is created to be used a handy interface
 *  to the function calls to PETSc. For each algebric equation 
 *  Ax=b, one PETSc instance is needed, since this instance has 
 *  the storage of these three variables. 
*/ 
#include "solver.h"

PETSc::PETSc()
{
	x = nullptr;			/* approx solution, RHS*/
	b = nullptr;
  	A = nullptr;            		/* linear system matrix */
  	
  	ksp = nullptr;        		/* Krylov subspace method context */
	nullsp = nullptr;
	pc = nullptr;
	
	KSPCreate(PETSC_COMM_WORLD,&ksp);
}

PETSc::PETSc(int ilower, int iupper, int d_nz, int o_nz)
{	
	x = nullptr;      			/* approx solution, RHS*/
	b = nullptr;
  	A = nullptr;            		/* linear system matrix */
  	
  	ksp = nullptr;          		/* Krylov subspace method context */
	nullsp = nullptr;
	pc = nullptr;
	Create(ilower, iupper, d_nz, o_nz);	
	KSPCreate(PETSC_COMM_WORLD,&ksp);
}

void PETSc::Create(int ilower, int iupper, int d_nz, int o_nz)
{	
	Create(PETSC_COMM_WORLD, ilower, iupper, d_nz, o_nz);	
}

void PETSc::Create(
	MPI_Comm Comm, 
	int ilower, 
	int iupper, 
	int d_nz, 
	int o_nz)
{	

    //TODO: Figure out if there is a good reason iupper is passed
    //      in as 'iupper -1' (see solve() functions in ellip.cpp)
    //      just to have 1 added back here, everywhere else a range
    //      of indices is calculated (e.g. see get_X()).
    //      If not, simplify. We should not have to do index bookkeeping
    //      both inside and outside the class. Preferably only inside
    //      class member functions.
	int n	= iupper - ilower + 1;
	
	comm 	= Comm;
	iLower	= ilower;	
	iUpper 	= iupper;	
	
    /*MatCreateAIJ(PETSC_COMM_WORLD,n,n,PETSC_DETERMINE,PETSC_DETERMINE,
				d_nz,PETSC_NULL,o_nz,PETSC_NULL,&A);*/
    MatCreateAIJ(PETSC_COMM_WORLD,n,n,PETSC_DECIDE,PETSC_DECIDE,
       d_nz,PETSC_NULL,o_nz,PETSC_NULL,&A);

    //TODO: See petsc manual for optimal sparse matrix
    //      creation and preallocation paradigm, similar
    //      to the following:
    //
    //       MatCreate(...,&A);
    //       MatSetType();
    //       MatSetFromOptions();
    //       MatXXXXSetPreallocation(); one of the below functions
    //
    //       MatSeqAIJSetPreallocation(); if serial
    //       MatMPIAIJSetPreallocation(); if parallel
	
    PetscObjectSetName((PetscObject) A, "A");
    MatSetFromOptions(A);
        //MatSetUp(A);//TODO: Need this???
	
    //TODO: use VecCreateMpi()?
	VecCreate(PETSC_COMM_WORLD, &b);
	PetscObjectSetName((PetscObject) b, "b");
	VecSetSizes(b, n, PETSC_DECIDE);
    VecSetFromOptions(b);
	
	VecCreate(PETSC_COMM_WORLD,&x);
	PetscObjectSetName((PetscObject) x, "X"); 
	VecSetSizes(x, n, PETSC_DECIDE);
    VecSetFromOptions(x);
}

PETSc::~PETSc()
{
	if(x!=nullptr)
	{
		VecDestroy(&x);
		x = nullptr;
	}
	if(b!=nullptr)
	{
		ierr =VecDestroy(&b);
		b = nullptr;
	}
	if(A!=nullptr)
	{
		MatDestroy(&A);
		A = nullptr;
	}
	if(ksp!=nullptr)
	{
		KSPDestroy(&ksp);
		ksp = nullptr;
	}
	if(nullsp!=nullptr)
	{
		MatNullSpaceDestroy(&nullsp);
		nullsp = nullptr;
	}
}

void PETSc::Reset_A()	// Reset all entries to zero ;
{
	MatZeroEntries(A);
}
void PETSc::Reset_b()  //  Reset all entries to zero ;
{
    VecZeroEntries(b);
}
void PETSc::Reset_x()
{
    VecZeroEntries(x);
}

void PETSc::FlushMatAssembly_A()
{
    ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY);
    ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY);
}

/*
void PETSc::Set_A(PetscInt m, PetscInt* Iids, PetscInt n, PetscInt* Jids, double* vals)
{
	ierr = MatSetValues(A,m,Iids,n,Jindices,vals,INSERT_VALUES);
}
*/

void PETSc::Set_A(PetscInt i, PetscInt j, double val)	// A[i][j]=val;
{
	ierr = MatSetValues(A,1,&i,1,&j,&val,INSERT_VALUES);
}

//TODO: add calls to FlushMatAssembly_A()
void PETSc::Add_A(PetscInt i, PetscInt j, double val)	// A[i][j]+=val;
{	
    //FlushMatAssembly_A();
	ierr = MatSetValues(A,1,&i,1,&j,&val,ADD_VALUES);
    //FlushMatAssembly_A();
}

void PETSc::Get_row_of_A(PetscInt i, PetscInt *ncol, PetscInt **cols, double **row)
{	
	ierr = MatGetRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
    //TODO: can't call MatRestoreRow() until after using the values in row.
	ierr = MatRestoreRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
    //At this point all the values in row have been zeroed out,
    //and the calling function will always get an array of zeros.
}

// x
void PETSc::Set_x(PetscInt i, double val)	// x[i]=val;
{
	ierr = VecSetValues(x,1,&i,&val,INSERT_VALUES);	
}

void PETSc::Add_x(PetscInt i, double val)	// x[i]+=val;
{
	ierr = VecSetValues(x,1,&i,&val,ADD_VALUES);
}

void PETSc::Set_b(PetscInt i, double val)	// x[i]=val;
{
	ierr = VecSetValues(b,1,&i,&val,INSERT_VALUES);
}

void PETSc::Add_b(PetscInt i, double val)	// x[i]+=val;
{
	ierr = VecSetValues(b,1,&i,&val,ADD_VALUES);
}

void PETSc::Get_x(double *p)
{
	PetscScalar      *values;
	VecGetArray(x,&values);
	for(int i = 0; i < iUpper-iLower+1; i++)
		p[i] = values[i];	
    VecRestoreArray(x,&values); 
}

void PETSc::Get_b(double *p)
{
	PetscScalar      *values;
	VecGetArray(b,&values);
	for(int i = 0; i < iUpper-iLower+1; i++)
		p[i] = values[i];	
    VecRestoreArray(b,&values); 
}

void PETSc::Get_x(double *p, 
	int n, 
	int *global_index)
{
    printf("ERROR: not implemented\n");
    LOC(); clean_up(EXIT_FAILURE);
}

void PETSc::SetPrevSolnInitialGuess()
{
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
}

void PETSc::SetMaxIter(int val)
{
	PetscInt maxits;
	double rtol, atol, dtol;
	
	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);
	ierr = KSPSetTolerances(ksp, rtol, atol, dtol, val);
}	/* end SetMaxIter */

void PETSc::SetTolerances(double rel_tol, double abs_tol, double div_tol)
{
	PetscInt maxits;
	double rtol, atol, dtol;
	KSPGetTolerances(ksp,&rtol,&atol,&dtol,&maxits);
    ierr = KSPSetTolerances(ksp,rel_tol,abs_tol,div_tol,maxits);
}
	
void PETSc::SetTol(double rel_tol)
{
	PetscInt maxits;
	double rtol, atol, dtol;
	
	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);

    //TODO: this only sets rtol (rel tol). Incorrect use in code???
    //      elliptic solver diverges during projection method
    //      if we use actual abs tolerance like below.
	
    ierr = KSPSetTolerances(ksp, rel_tol, atol, dtol, maxits);
	
    //TODO: Crashes projection method elliptic solver for poisson eqn.
    //      if the absolute tolerance, atol, is used.
    //
    //ierr = KSPSetTolerances(ksp, rtol, val, dtol, maxits);
}

void PETSc::SetKDim(int val)
{
    printf("ERROR: not implemented\n");
    LOC(); clean_up(EXIT_FAILURE);
}

void PETSc::GetNumIterations(PetscInt *num_iterations)
{
	KSPGetIterationNumber(ksp,num_iterations);        
}	/* end GetNumIterations */

void PETSc::GetResidualNorm(double *resid_norm)
{
	KSPGetResidualNorm(ksp,resid_norm);
}	/* end GetResidualNorm */

//TODO: This is also not the relative residual norm,
//      which is defined as ||Ax-b||/||b||
void PETSc::GetFinalRelativeResidualNorm(double *rel_resid_norm)
{
	KSPGetResidualNorm(ksp,rel_resid_norm);
}	/* end GetFinalRelativeResidualNorm */


void PETSc::Solve_PetscDecide()
{
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    KSPSetOperators(ksp,A,A);
    KSPSolve(ksp,b,x);
}

//TODO: Try tuning parameters, restart etc.
void PETSc::Solve_GMRES(void)
{
    start_clock("Assemble matrix and vector");
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	stop_clock("Assembly matrix and vector");

    KSPSetOperators(ksp,A,A);
	KSPSetType(ksp,KSPGMRES);

    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

	start_clock("KSPSolve");
    KSPSolve(ksp,b,x);
	stop_clock("KSPSolve");

}	/* end Solve_GMRES */

void PETSc::Solve(void)
{
#if defined HAVE_HYPRE
	Solve_HYPRE();
#else // defined HAVE_HYPRE*/
	Solve_BCGSL();
#endif // defined HAVE_HYPRE
}	/* end Solve */

/*
From the Hypre user manual:

For three-dimensional diffusion problems, it is recommended
to choose a lower complexity coarsening like HMIS or PMIS (coarsening 10 or 8) and combine it
with a distance-two interpolation (interpolation 6 or 7), that is also truncated to 4 or 5 elements
per row. Additional reduction in complexity and increased scalability can often be achieved using
one or two levels of aggressive coarsening
*/
#if defined HAVE_HYPRE
void PETSc::Solve_HYPRE(void)
{
    PC pc;
    start_clock("Assemble matrix and vector");
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    ierr = VecAssemblyBegin(x);
    ierr = VecAssemblyEnd(x);

    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);
    stop_clock("Assembly matrix and vector");

    //TODO: Can we still use KSPBCGSL with the
    //      boomeramg preconditioner? 
	KSPSetType(ksp,KSPBCGS);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
	PCSetType(pc,PCHYPRE);
    PCHYPRESetType(pc,"boomeramg");
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

    start_clock("KSPSolve");
    KSPSolve(ksp,b,x);
    stop_clock("KSPSolve");
}
#endif // defined HAVE_HYPRE

void PETSc::Solve_BCGSL(void)
{
    start_clock("Assemble matrix and vector");
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	stop_clock("Assembly matrix and vector");

    KSPSetOperators(ksp,A,A);
    KSPSetType(ksp,KSPBCGSL);
    
    //sets the number of search directions for BiCGStab(L)
	KSPBCGSLSetEll(ksp,2);

    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

	start_clock("KSPSolve");
    KSPSolve(ksp,b,x);
	stop_clock("KSPSolve");
}

void PETSc::Solve_withPureNeumann_GMRES(void)
{
	if (debugging("solver"))
	    printf("Entering Solve_withPureNeumann_GMRES()\n");
    	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
  	
    //NOTE: This sets the nullspace to contain the constant vector
    //      (second argument is PETSC_TRUE)
	MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
    MatSetNullSpace(A,nullsp);
	MatNullSpaceRemove(nullsp,b);
	
    KSPSetOperators(ksp,A,A);
	KSPSetType(ksp,KSPGMRES);

    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

	start_clock("Petsc Solve in pure neumann solver");
    KSPSolve(ksp,b,x);
	stop_clock("Petsc Solve in pure neumann solver");

    if (debugging("solver"))
	    printf("Leaving Solve_withPureNeumann_GMRES()\n");
}	/* end Solve_withPureNeumann_GMRES */

void PETSc::Solve_withPureNeumann(void)
{
#ifdef HAVE_HYPRE
	Solve_withPureNeumann_HYPRE();
#else
	Solve_withPureNeumann_BCGSL();
#endif
}

#if defined HAVE_HYPRE
void PETSc::Solve_withPureNeumann_HYPRE(void)
{
    if (debugging("solver"))
        printf("Entering Solve_withPureNeumann_HYPRE()\n");

    PC pc;
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    ierr = VecAssemblyBegin(x);
    ierr = VecAssemblyEnd(x);

    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);

    //NOTE: This sets the nullspace to contain the constant vector
    //      (second argument is PETSC_TRUE)
    MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
    MatSetNullSpace(A,nullsp);
    MatNullSpaceRemove(nullsp,b);

    KSPSetType(ksp,KSPBCGS);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCHYPRE);
    PCHYPRESetType(pc,"boomeramg");
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    
    start_clock("Petsc Solve in pure neumann solver");
    KSPSolve(ksp,b,x);
    stop_clock("Petsc Solve in pure neumann solver");

    if (debugging("solver"))
        printf("Leaving Solve_withPureNeumann_HYPRE()\n");
}
#endif // defined HAVE_HYPRE

void PETSc::Solve_withPureNeumann_BCGSL(void)
{
	if (debugging("solver"))
	    printf("Entering Solve_withPureNeumann_BCGSL()\n");
	
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  	
  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);
  	
  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	
    //NOTE: This sets the nullspace to contain the constant vector
    //      (second argument is PETSC_TRUE)
	MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
    MatSetNullSpace(A,nullsp);
    MatNullSpaceRemove(nullsp,b);
	
    KSPSetOperators(ksp,A,A);
	KSPSetType(ksp,KSPBCGSL);

    //sets the number of search directions for BiCGStab(L)
	KSPBCGSLSetEll(ksp,2);

    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

	start_clock("Petsc Solve in pure neumann solver");
    KSPSolve(ksp,b,x);
	stop_clock("Petsc Solve in pure neumann solver");
	
    if (debugging("solver"))
	    printf("Leaving Solve_withPureNeumann_BCGSL()\n");
}	/* end Solve_withPureNeumann_BCGSL */

void PETSc::Print_A(const char *filename)
{
	PetscViewer viewer;
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);
}	/* end Print_A */

void PETSc::Print_b(const char *filename)
{
    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,
            PETSC_VIEWER_ASCII_MATLAB);
    VecView(b, PETSC_VIEWER_STDOUT_WORLD);
}	/* end Print_b */

extern void viewTopVariable(
	Front *front,
	double *var,
	boolean set_bounds,
	double var_min,
	double var_max,
	char *dirname,
	char *var_name)
{
	HDF_MOVIE_VAR hdf_movie_var;
	HDF_MOVIE_VAR *hdf_movie_var_save = front->hdf_movie_var;
	front->hdf_movie_var = &hdf_movie_var;
	hdf_movie_var.num_var = 1;
	FT_MatrixMemoryAlloc((POINTER*)&hdf_movie_var.var_name,1,100,
				sizeof(char));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.top_var,1,
				sizeof(double*));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.preset_bound,1,
				sizeof(boolean));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.var_min,1,
				sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&hdf_movie_var.var_max,1,
				sizeof(double));
	sprintf(hdf_movie_var.var_name[0],"%s",var_name);
	hdf_movie_var.preset_bound[0] = set_bounds;
	hdf_movie_var.var_min[0] = var_min;
	hdf_movie_var.var_max[0] = var_max;
	hdf_movie_var.top_var[0] = var;
	gview_var2d_on_top_grid(front,dirname);

	FT_FreeThese(5,hdf_movie_var.var_name,hdf_movie_var.top_var,
				hdf_movie_var.preset_bound,
				hdf_movie_var.var_min,hdf_movie_var.var_max);
	front->hdf_movie_var = hdf_movie_var_save;
}	/* end viewTopVariable */

void PETSc::Solve_LU(void)
{
	PC pc;
    start_clock("Assemble matrix and vector");
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    ierr = VecAssemblyBegin(x);
    ierr = VecAssemblyEnd(x);

    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);
    stop_clock("Assembly matrix and vector");

    KSPSetType(ksp,KSPPREONLY);
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCLU);
    KSPSetOperators(ksp,A,A);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

    start_clock("KSPSolve");
    KSPSolve(ksp,b,x);
    stop_clock("KSPSolve");
} /*direct solver, usually give exact solution for comparison*/
