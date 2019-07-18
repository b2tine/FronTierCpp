/**********************************************************************
 * 		solver.h
 **********************************************************************/

#ifndef FT_IFLUID_SOLVER_H
#define FT_IFLUID_SOLVER_H

#include <FronTier.h>
#include <vector>
#include <assert.h>

#include <petscksp.h>
#include <petscmat.h>
#include <petscpc.h>

#include <algorithm>

enum
{
    NO_PDE_BOUNDARY = 0,
    CONST_V_PDE_BOUNDARY = 1,
    CONST_P_PDE_BOUNDARY,
    NEUMANN_PDE_BOUNDARY,
    DIRICHLET_PDE_BOUNDARY,
    MOVING_BOUNDARY,
    MIXED_PDE_BOUNDARY
};

class SOLVER
{
public:
	SOLVER(){};
	SOLVER(int ilower, int iupper, int d_nz, int o_nz){};
	virtual ~SOLVER(){};
	virtual void Create(int ilower, int iupper, int d_nz, int o_nz){};
	
	virtual void Set_A(PetscInt i, PetscInt j, double val){};
						// A[i][j]=val;
	virtual void Add_A(PetscInt i, PetscInt j, double val){};
						// A[i][j]=A[i][j]+val;
	virtual void Set_x(PetscInt i, double val){};	// x[i]=val;
	virtual void Set_x(double *p){};		// x[i]=p[i];
	virtual void Add_x(PetscInt i, double val){};	// x[i]=x[i]+val;
	virtual void Get_x(double *p){};	// get the x from ij_x to p.		
	virtual void Get_x(double *p, int n, int *global_index){};
	virtual void Set_b(PetscInt i, double val){};	// b[i]=val;
	virtual void Set_b(double *b){};	
	virtual void Add_b(PetscInt i, double val){};	// b[i]=b[i]+val;

	virtual void SetMaxIter(int val){};	
	virtual void GetFinalRelativeResidualNorm(double *rel_resid_norm){};
	virtual void GetNumIterations(int *num_iterations){};

	virtual void Solve(void){};	
	virtual void Solve_withPureNeumann(void){};	
	virtual void Read_A(char *filename){};
	virtual void Print_A(char *filename){};
	virtual void Read_b(char *filename){};
	virtual void Print_b(char *filename){};
	virtual void Read_x(char *filename){};
	virtual void Print_x(char *filename){};
	virtual void test(void){};
};

//class PETSc: public SOLVER
class PETSc
{
public:	
	MPI_Comm  comm;			// set to be MPI_COMM_WORLD.
	int iLower;
	int iUpper;			// global row range
	
	Vec x;      			/* approx solution, RHS*/
	Vec b;
  	Mat A;          		/* linear system matrix */
  	
  	KSP   ksp;          		/* Krylov subspace method context */
	PC    pc;
	MatNullSpace	nullsp;
		
	PetscErrorCode ierr;
	int its;			// numer of iterations;

public:
	PETSc();
	PETSc(int ilower, int iupper, int d_nz, int o_nz);		
		// global row range of A, x, b on this processor
	~PETSc();
	void Create(int ilower, int iupper, int d_nz, int o_nz);	
		// same as Hypre(int, int)
	void Create(MPI_Comm Comm, int ilower, int iupper, int d_nz, int o_nz);
		// same as Hypre(int, int)
	
	void Reset_A();				// Set A[i][j]=0.0;
	void Reset_b();
	void Reset_x();
	void Set_A(PetscInt i, PetscInt j, double val);	// A[i][j]=val;
    //void Set_A(PetscInt m, PetscInt* Iids, PetscInt n, PetscInt* Jids, double* vals);
	void Add_A(PetscInt i, PetscInt j, double val);	// A[i][j]=A[i][j]+val;
	void Get_row_of_A(PetscInt i, PetscInt *ncol, PetscInt **cols, double **row);
	void Set_x(PetscInt i, double val);		// x[i]=val;
	void Add_x(PetscInt i, double val);		// x[i]=x[i]+val;
	void Set_b(PetscInt i, double val);		// b[i]=val;
	void Add_b(PetscInt i, double val);		// b[i]=b[i]+val;
	void Get_x(double *p);		// get the x from ij_x to p.	
	void Get_b(double *p);		// get the b from ij_x to p.
	void Get_x(double *p, int n, int *global_index);
	
	void SetMaxIter(int val); 	// Set maximum number of iterations 
	void SetTol(double val);	// Set the convergence tolerance 
	void SetKDim(int k_dim);	
			// Set the maximum size of the Krylov space 
	void GetNumIterations(PetscInt *num_iterations);	
			// Return the number of iterations taken 
	void GetFinalRelativeResidualNorm(double *rel_resid_norm);
	
    void Solve_PetscDecide();
    void Solve(void);
	void Solve_GMRES(void);
	void Solve_BCGSL(void);
	void Solve_LU(void);
#if defined(HAVE_HYPRE)
	void Solve_HYPRE(void);
#endif
	void Solve_withPureNeumann(void);
	void Solve_withPureNeumann_GMRES(void);
	void Solve_withPureNeumann_BCGSL(void);
	void Solve_withPureNeumann_HYPRE(void);
	virtual void Print_A(const char *filename);
        virtual void Print_b(const char *filename);
};

class PARABOLIC_SOLVER{
        Front *front;
public:
        PARABOLIC_SOLVER(Front &front);

	// Time step for one call
	double dt;
	double sub_dt;

        // On topological grid
        boolean first;
	int *i_to_I;
	int **ij_to_I;
	int ***ijk_to_I;
	int ilower;
	int iupper;
	int order;		/* order of Runge-Kutta */
	CELL_PART *cell_part;	/* Cell partition structure */

	COMPONENT soln_comp;
	COMPONENT obst_comp;
	double *var;		/* field variable of old step */
	double *soln;		/* field variable of new step */
	double *source;		/* source field */
	double **a;		/* advection field */
	double D;
	double *nu;             /* Variable diffusion coefficient */
	double var_obst;	/* default solution in obst_comp */
	double (*var_obst_func)(POINTER,double*); 
				/* function for obst_comp solution */
	POINTER obst_func_params; 
				/* parameters for var_obst_func */
	double bdry_influx;
	double (*getStateVarFunc)(POINTER);
	double (*getStateVel[3])(POINTER);
	void set_solver_domain(void);

	void solveEX(void);
	void solveCN(void);
	void solve(double *var_in,double *var_out);
	void solve1d(double *var_in,double *var_out);
	void solve2d(double *var_in,double *var_out);
	void solve3d(double *var_in,double *var_out);
	void solveIM(void);

	void solveCEX(void);
	void solveCCN(void);
	void solveCIM(void);
	void Csolve(double *var_in,double *var_out);
	void Csolve1d(double *var_in,double *var_out);
	void Csolve2d(double *var_in,double *var_out);
	void Csolve3d(double *var_in,double *var_out);
	void CCNsolve();
	void CCNsolve1d();
	void CCNsolve2d();
	void CCNsolve3d();

	double compBdryFlux(double*,double);
	void addMeshState(double *ans,double C1,double *var1,double C2,
				double *var2);
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	int (*findCrossingInfo)(Front*,int*,GRID_DIRECTION,int,
                                double*,	// crx_old
				double*,	// crx_new
				double*);	// flux
private:
        // Dimension
        int dim;
        COMPONENT *top_comp;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;
	int array_size;
        double *array;          // for scatter states;
	double *top_h;
	double *top_L;
	double var_min;
        double var_max;
        double checkSolver(int *icoords,boolean print_details,double *var_in);
};

class DOUBLE_ELLIPTIC_SOLVER{
        Front *front;
public:
        DOUBLE_ELLIPTIC_SOLVER(Front &front);
        ~DOUBLE_ELLIPTIC_SOLVER();

        // On topological grid
	int **dij_to_I;
	int ***dijk_to_I;
	int eilower;
	int eiupper;
        int *ext_gmax;
        int *ext_l,*ext_u;
        int *ext_imin,*ext_imax;
        COMPONENT *ext_comp;

        double dt;              //time step
	double porosity;
	double *soln;	        // field variable of new step
	double *source;	        // source field
        double *D;              // div(D*grad)phi = source,  where D = 1.0/rho
        double *ext_array;

	void set_solver_domain(void);
        void set_extension(void);
	void dsolve(double *soln);
	
        double (*getStateVar)(POINTER);
        double (*getStateVel[3])(POINTER);

	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	double checkSolver(int *icoords,boolean print_details);
	int skip_neumann_solver;
private:
        // Dimension
        int dim;
        COMPONENT *top_comp;
	double *top_h;
        int *top_gmax;
	int imin[MAXD],imax[MAXD];
        double *ext_source;             // for extended source;
        double *ext_D;                  // for extended D;
	int array_size;
	double max_soln;
	double min_soln;
        double max_rhs;
	void dsolve2d(double *soln);
	void dsolve3d(double *soln);
        bool icoordsInterior(int*);
        bool icoordsBoundary(int*);
        double dcheckSolver(int*,boolean);
        double dcheckSolverInterior(int*,boolean);
        double dcheckSolverExtended(int*,boolean);
};

class ELLIPTIC_SOLVER{
        Front *front;
public:
        ELLIPTIC_SOLVER(Front &front);

        // On topological grid
	int *i_to_I;
	int **ij_to_I;
	int ***ijk_to_I;
	int ilower;
	int iupper;

        double dt;          //time step
	double porosity;
	double *soln;		/* field variable of new step */
	double *source;		/* source field */
        double **vel;       /* velocity field */
        double *D;          /* div(D*grad)phi = source,  where D = 1.0/rho */

	void set_solver_domain(void);
	void solve(double *soln);
	void dsolve(double *soln);
	
        double (*getStateVar)(POINTER);
        double (*getStateVel[3])(POINTER);

	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	double checkSolver(int *icoords,boolean print_details);
        void printIsolatedCells();
	int skip_neumann_solver;
private:
        // Dimension
        int dim;
        COMPONENT *top_comp;
	double *top_h;
	double *top_L;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;
        double *array;          // for scatter states;
	int array_size;
	double max_soln;
	double min_soln;
	void solve1d(double *soln);
	void solve2d(double *soln);
	void solve3d(double *soln);
};

struct _SWEEP {
	double **vel;
	double *rho;
};
typedef struct _SWEEP SWEEP;

struct _FSWEEP {
	double **vel_flux;
};
typedef struct _FSWEEP FSWEEP;

class HYPERB_SOLVER{
        Front *front;
public:
        HYPERB_SOLVER(Front &front);

	// Time step for one call
	double dt;

	int order;		/* order of Runge-Kutta */
	int size;

	COMPONENT soln_comp1;
	COMPONENT soln_comp2;
	COMPONENT obst_comp;
	double porosity;
	double **var;		/* field variable of old step */
	double **soln;		/* field variable of new step */
	double **source;	/* source field */
	double *rho;
	double rho1;
	double rho2;
	double var_obst;	/* default solution in obst_comp */
	double max_speed;
	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
	void (*numericalFlux)(SWEEP*,FSWEEP*,double,int,int,int);
	void solveRungeKutta();
	double (*getStateVel[3])(POINTER);
private:
        // Dimension
        int dim;
        COMPONENT *top_comp;
	double *top_L,*top_h;
        int *top_gmax;
	int imin,jmin,kmin;
	int imax,jmax,kmax;
	int nrad;
        double *array;          // for scatter states;
	SWEEP *st_field,st_tmp;
        FSWEEP *st_flux;
        double **a,*b;
	void setSolverDomain(void);
	void allocMeshVst(SWEEP*);
	void allocMeshFlux(FSWEEP*);
	void allocDirectionalVstFlux(SWEEP*,FSWEEP*);
	void resetFlux(FSWEEP*);
	void copyToMeshVst(SWEEP*);
	void copyMeshVst(SWEEP,SWEEP*);
	void computeMeshFlux(SWEEP,FSWEEP*);
	void addMeshFluxToVst(SWEEP*,FSWEEP,double);
	void copyFromMeshVst(SWEEP);
	void addFluxInDirection(int,SWEEP*,FSWEEP*);
	void addFluxInDirection1d(int,SWEEP*,FSWEEP*);
	void addFluxInDirection2d(int,SWEEP*,FSWEEP*);
	void addFluxInDirection3d(int,SWEEP*,FSWEEP*);
	void appendGhostBuffer(SWEEP*,SWEEP*,int,int*,int,int);
	void setNeumannStates(SWEEP*,SWEEP*,HYPER_SURF*,POINTER,int*,int,
				int,int,int,int);
	void setDirichletStates(SWEEP*,SWEEP*,HYPER_SURF*,POINTER,int*,int,
				int,int,int,int);
	void setElasticStates(SWEEP*,SWEEP*,HYPER_SURF*,POINTER,int*,int,
				int,int,int,int);
	void addSourceTerm(SWEEP*,FSWEEP*);
};


extern	void upwind_flux(SWEEP*,FSWEEP*,double,int,int,int);
extern	void weno5_flux(SWEEP*,FSWEEP*,double,int,int,int);

extern void viewTopVariable(Front*,double*,boolean,double,double,char*,char*);
extern double   compBdryFlux(Front*,double*,double,int,double,
		double(*)(POINTER));
#endif
