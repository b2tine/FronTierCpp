/**********************************************************************
 * 		iFluid.h
 **********************************************************************/

#ifndef FT_IFLUID_H
#define FT_IFLUID_H

#include <vector>
#include <assert.h>

#include <solver.h>
#include "ifluid_state.h"
#include "rigidbody.h"

#define SOLID_COMP		0
#define LIQUID_COMP1	2
#define LIQUID_COMP2	3
#define LIQUID_COMP		3
#define	FILL_COMP		10

#define	ifluid_comp(comp) (((comp) == LIQUID_COMP1 || 	\
		comp == LIQUID_COMP2) ? YES : NO)

enum _IF_PROB_TYPE {
        ERROR_TYPE = -1,
	BEE_3D = 1,
        BUBBLE_SURFACE,
	CHANNEL_FLOW,
        FLUID_CRYSTAL,
        FLUID_RIGID_BODY,
        FLUID_SOLID_CIRCLE,
	FLUID_SOLID_CONE,
        FLUID_SOLID_CYLINDER,
        FLUID_SOLID_RECT,
        FLUID_SOLID_TRIANGLE,
	HELICOPTER_3D,
	OPEN_ROTOR,
	PRESSURE_PUMP,
	RANDOM_FLOW,
        ROTOR_ONE_FLUID,
        ROTOR_TWO_FLUID,
	TAYLOR_GREEN_VORTEX,
        TWO_FLUID_BUBBLE,
        TWO_FLUID_KH,
        TWO_FLUID_RT,
	WINDMILL_2D,
        WINDMILL_3D,
	HUMAN_BODY_3D
};
typedef enum _IF_PROB_TYPE IF_PROB_TYPE;

enum EBM_COORD
{
    COORD_X = 0,  COORD_Y = 1,  COORD_Z = 2
};

enum _DOMAIN_STATUS {
        NOT_SOLVED      =       0,
        TO_SOLVE,
        SOLVED
};
typedef enum _DOMAIN_STATUS DOMAIN_STATUS;

struct IF_FIELD {
	double **vel;			/* Velocities */
	double **vorticity;		/* 3d Vorticity vector */
	double *temperature;            /* Temperature */
	double *phi;
	double **grad_phi;
	double *q;
	double *pres;			/* Pressure */
	double *vort;			/* Magnitude of Vorticity in 2D */
	double *mu;
	double *rho;
	double **grad_q;
	double **f_surf;		// Surface force (such as tension)
	double **old_var;		// For debugging purpose

        //double *d_phi;          //Dual grid phi
	double *div_U;
	double *nu_t;			/* Turbulent viscosity */
	double **ext_accel;		/*external forcing from other field*/
};

enum _PROJC_METHOD {
	ERROR_PROJC_SCHEME		= -1,
        SIMPLE			=  1,
        BELL_COLELLA,
        KIM_MOIN,
        PEROT_BOTELLA,
        PMI,
        PMII,
        PMIII
};
typedef enum _PROJC_METHOD PROJC_METHOD;

enum _ADVEC_METHOD {
	ERROR_ADVEC_SCHEME		= -1,
        UPWIND			=  1,
        WENO
};
typedef enum _ADVEC_METHOD ADVEC_METHOD;

enum _ELLIP_METHOD {
	ERROR_ELLIP_SCHEME		= -1,
	SIMPLE_ELLIP		= 1,
	DOUBLE_ELLIP,
    DUAL_ELLIP
};

typedef enum _ELLIP_METHOD ELLIP_METHOD;

enum _DOMAIN_METHOD {
	BY_COMPONENT		= 1,
	BY_CROSSING,
	BY_WALL
};
typedef enum _DOMAIN_METHOD DOMAIN_METHOD;

enum _EDDY_VISC {
	BALDWIN_LOMAX		= 1,
	VREMAN,
	SMAGORINSKY,
    KEPSILON
};
typedef enum _EDDY_VISC EDDY_VISC;

struct _NS_SCHEME {
	PROJC_METHOD projc_method;
	ADVEC_METHOD advec_method;
	ELLIP_METHOD ellip_method;
};
typedef struct _NS_SCHEME NS_SCHEME;

struct FINITE_STRING {         // For fluid drag on string chord
    double radius;
    double dens;
    double c_drag;
    double ampFluidFactor;
};

//vortex params
struct VPARAMS {
    double center[MAXD];            // center of vortex
    double D;                       // size of vortex
    double A;                       // intensity of vortex
};

struct IF_PARAMS
{
    int dim;
    POINTER level_func_params;
	NS_SCHEME num_scheme;
    double rho1;
    double rho2;
	double mu1;
	double mu2;
	double U1[MAXD];
	double U2[MAXD];
	double gravity[MAXD];
	double U_ambient[MAXD];

	double surf_tension;
	double smoothing_radius {1.0};
	
    double ub_speed;
	double min_speed;	/* Limit time step in zero ambient velocity */
	COMPONENT m_comp1;
	COMPONENT m_comp2;
	
    IF_FIELD *field;

	int adv_order;
	boolean total_div_cancellation;
	boolean buoyancy_flow {NO};
	boolean if_buoyancy {NO};
	double  ref_temp;
	boolean if_ref_pres {NO};
	double  ref_pres;
	double  Amplitute; 	/*Amplitute of velocity*/
	boolean  with_porosity;    /*porosity: 1/0 with/without porosity*/
        
    double  porous_coeff[2];   /*dp = a*v + b*v^2*/
	double	porosity;
	char base_dir_name[200];
    int base_step;
	boolean scalar_field; /*include scalar field or not*/
	boolean skip_neumann_solver;
    int fsi_startstep;

    //TODO: factor out turbulence params into separate data structure, eddy_params
	//POINTER eddy_params;
	
    EDDY_VISC eddy_visc_model;
	boolean use_eddy_visc;	/* Yes if to use eddy viscosity */
	double	ymax {0};	   	/* Maximum distance in Baldwin-Lomax model */
    double C_s;     //Smagorinsky model constant
    double C_v;     //Vreman model constant
};

struct _FLOW_THROUGH_PARAMS {
        POINT *oldp;
        COMPONENT comp;
};
typedef struct _FLOW_THROUGH_PARAMS FLOW_THROUGH_PARAMS;

enum _TIME_FUNC_TYPE {
	CONSTANT		=  1,
	PULSE_FUNC,
	SINE_FUNC	
};
typedef enum _TIME_FUNC_TYPE TIME_FUNC_TYPE;

struct _TIME_DEPENDENT_PARAMS {
	TIME_FUNC_TYPE td_type;
	double v_base[MAXD],p_base;
	double v_peak[MAXD],p_peak;
	double v_tail[MAXD],p_tail;
	double v_amp[MAXD],p_amp;
	double omega,phase;
	double T[10];
};
typedef struct _TIME_DEPENDENT_PARAMS TIME_DEPENDENT_PARAMS;

struct _SPLIT_STATE_PARAMS {
	int dir;
	double split_coord;
	STATE left_state;
	STATE right_state;
};
typedef struct _SPLIT_STATE_PARAMS SPLIT_STATE_PARAMS;

struct _PARABOLIC_STATE_PARAMS {
	int dir;
	double v_peak[MAXD];
	STATE state;
};
typedef struct _PARABOLIC_STATE_PARAMS PARABOLIC_STATE_PARAMS;

struct CUBIC_STATE_PARAMS {
    STATE state;
    double h;
    enum FLOW_TYPE {LAMINAR, TURBULENCE};
    FLOW_TYPE flow_type;
};

struct _OPEN_PIPE_PARAMS
{
        int dir;
        int side;
        double center[MAXD];
        double radius;
        int in_pipe_bdry;
        int out_pipe_bdry;
        boolean in_flow_through;
        boolean out_flow_through;
        STATE state[2];
};
typedef struct _OPEN_PIPE_PARAMS OPEN_PIPE_PARAMS;

/******************************************************************************
 * 		lcartsn.h
 * A simple incompressible flow solver using the ghost fluid method and the
 * projection method.
 *
 * the main function is 
 * 	L_CARTESIAN::solve().
 *
 * References:
 ******************************************************************************/

class SOLVER;
class Incompress_Solver_Basis;

class KE_PARAMS;

class IF_RECTANGLE
{
    public:
	
        int comp;			 
        int index;
        double area;
        double coords[MAXD]; // x,y,z data of mesh block's center
        int icoords[MAXD];     // i,j,k indices of mesh block

        IF_RECTANGLE();

        void setCoords(double*,int);
        std::vector<double> getCoords();
};

class Incompress_Solver_Basis{
public:
       	Incompress_Solver_Basis() {}; // constructor
	virtual ~Incompress_Solver_Basis() {};

};

class Incompress_Solver_Smooth_Basis : public Incompress_Solver_Basis
{
public:
        //constructor
	Incompress_Solver_Smooth_Basis(Front &front);
	virtual ~Incompress_Solver_Smooth_Basis() {};

	double m_dt;
	double accum_dt;
	double max_speed;
	double min_pressure;
        double max_pressure;
        double min_value; //for debugging
	double max_value; //for debugging
	double max_dt;
	double min_dt;
	double *top_h;
	double vmin[MAXD],vmax[MAXD];
	int dim;
	int icrds_max[MAXD];

	void initMesh(void);
	void computeMaxSpeed(void);
	void setAdvectionDt(void);
	
    void readFrontInteriorStates(char *state_name);
	void printFrontInteriorStates(char *state_name);
	void initMovieVariables(void);
	void getVelocity(double *p, double *U);
	void initSampleVelocity(char *in_name);
    
    void printEnstrophy();
    void printEnstrophy2d();
    void printEnstrophy3d();

	//Initialization of States
	void (*getInitialState) (COMPONENT,double*,IF_FIELD*,int,int,
				IF_PARAMS*);
	/*set initial velocity with one function, no loop needed*/
	void (*setInitialVelocity)(COMPONENT,int*,double*,double*,double*,
                                RECT_GRID*,IF_PARAMS*);

 	int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
				POINTER*,HYPER_SURF**,double*);
	int (*findStateAtCGCrossing)(Front*,int*,GRID_DIRECTION,int,
				POINTER*,HYPER_SURF**,double*);
	
    void applicationSetComponent();
	void applicationSetStates();
	
    double computeFieldPointPressureJump(int*,double,double);     
    void computeFieldPointGradJump(int*,double*,double*);

    
    void setSlipBoundary(int* icoords, int idir, int nb, int comp,
            HYPER_SURF* hs, POINTER state, double** vel, double* v_slip);
        
    void setSlipBoundaryNIP(int* icoords, int idir, int nb, int comp,
            HYPER_SURF* hs, POINTER state, double** vel, double* v_slip);

        void setSlipBoundaryGNOR(int* icoords, int idir, int nb, int comp,
                HYPER_SURF* hs, POINTER state, double** vel, double* v_slip);
    
    //For debugging test
	void compareWithBaseSoln(void);
        void readBaseFront(IF_PARAMS *,int i);
        void readBaseStates(char *restart_name);
	void solveTest(const char *msg);
        void addVortexDisturbance(const VPARAMS&);

	//User interface
	int skip_neumann_solver;
	virtual void setInitialCondition(void) = 0;
	virtual void setParallelVelocity(void) = 0;
	virtual void solve(double dt) = 0; // main step function
        virtual void vtk_plot_scalar(char*, const char*) = 0;
    virtual void writeMeshFileVTK();

    //std::priority_queue<IF_Injection*> InjectionEvents;
    //void scheduleInjectionEvent(IF_Injection*);
    //void consumeInjectionEvent(IF_Injection*);


protected:
	Front *front;
	Front *base_front;	// for convergence tests
	IF_PARAMS *iFparams;
	IF_FIELD  *field;
	
	// On dual topological grid
	RECT_GRID *top_grid;
	double *array;
	double *source;
	double *diff_coeff;
	COMPONENT *top_comp;
	int *top_gmax;
	int *lbuf, *ubuf;
	double *top_L, *top_U;
	int **ij_to_I, ***ijk_to_I;
	int *domain_status;
	int smin[MAXD],smax[MAXD];
	
    // Sweeping limits
	int imin, jmin, kmin;
	int imax, jmax, kmax;
	// for parallel partition
	int NLblocks, ilower, iupper;
	int *n_dist;
	
    // On Double solver
	COMPONENT *ext_comp;
	int ext_gmax[MAXD];
        int ext_l[MAXD],ext_u[MAXD];
        int D_extension;
	int **dij_to_I,***dijk_to_I;
	// Sweeping limites
	int ext_imin[MAXD];
	int ext_imax[MAXD];
	// for parallel partition
	int dNLblocks, eilower, eiupper;
	int *dn_dist;

	// Index shift between dual and comp grids 
	int ishift[MAXD];

    //TODO: should rename this to avoid confusion/collision with the macro in geom.h
	//member data: mesh storage
	std::vector<IF_RECTANGLE>   cell_center;

	//member data:
	int    m_comp[2];
	double m_mu[2];
	double m_rho[2];// two components at most
	double m_sigma; //surface tension
	double m_smoothing_radius;// used by smoothing function

	double hmin; //smallest spacing
	double mu_min; //smallest viscocity
	double rho_min;// smallest density
	double m_t;
	double m_t_old, m_t_int, m_t_new;

	//following is for new painting algorithm
	std::vector<int> color_map; 		//big array
	std::vector<int> color_map_copy; 	//big array
	int drawColorMap();
	void paintInterior(int&);
	void makeGlobalColorMap(int&);
	void paintConnectedRegion(int,int);
	boolean paintToSolveGridPoint2(int);
protected:
	void setComponent(void); //init components;
	void setDomain();
	void setDoubleDomain();

	// parallelization related functions
	void scatMeshArray(void);
	void setGlobalIndex(void);
	void setDoubleGlobalIndex(void);
	void setIndexMap(void);
	void setDoubleIndexMap(void);
	void paintAllGridPoint(int status);
	void paintSolvedGridPoint();
	void setReferencePressure();
        void setIsolatedSoln(int,double*);
	boolean paintToSolveGridPoint();
	boolean nextConnectedPoint(int*,GRID_DIRECTION,int*,int,int*,int*);

/****************  Functions related to solve() *********/

	virtual void copyMeshStates(void) = 0;
	virtual void computeGradientQ(void) = 0;
	virtual void computeVelDivergence(void) = 0;
	virtual void surfaceTension(double*, HYPER_SURF_ELEMENT*,
		HYPER_SURF*, double*, double) = 0;

	
/***********************  Utility functions  *******************/

	void   getRectangleIndex(int indexRectangle, int &i, int &j);
	void   getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int    getRectangleComponent(int index);	// the center component
	void   getRectangleCenter(int index, double *coords);
	double getDistance(double *coords0, double *coords1);
	int    getComponent(int *icoords);	
	int    getComponent(double *coords);	
	void   save(char *filename);
	double computeFieldPointDiv(int*, double**);
	double computeFieldPointDivSimple(int*, double**);
	double computeFieldPointDivDouble(int*, double**);
	double computeMuOfBaldwinLomax(int*, double, boolean);
	double computeMuOfVremanModel(int*);
	double computeMuofSmagorinskyModel(int*);
	KE_PARAMS* computeMuOfKepsModel();
	
    std::vector<std::vector<double>> computeVelocityGradient(int* icoords);

    void computeFieldPointGrad(int* icoords, double* field, double* grad_field);
    void computeFieldPointGradQ(int* icoords, double* field, double* grad_field);

	void checkVelocityDiv(const char*);
/************* TMP Functions which are not implemented or used ***********/

	void computeSubgridModel(void);    // subgrid model by Hyunkyung Lim
	void getNearestInterfacePoint(COMPONENT,double*,double*,double*,
					double*); 
};

///////////////Interface for Embedded Boundary Method////////////////////

class Incompress_Solver_EBM:public Incompress_Solver_Basis{
public:
        Incompress_Solver_EBM(Front &front) {};//constructor
	~Incompress_Solver_EBM() {};
};
/////////////////////////////////////////////////////////////////////////////////
//
class Incompress_Solver_Smooth_2D_Basis:
public 	Incompress_Solver_Smooth_Basis{
public:
        Incompress_Solver_Smooth_2D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_2D_Basis() {};

protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();
	void setSmoothedProperties(void);
};


class Incompress_Solver_Smooth_3D_Basis:
public 	Incompress_Solver_Smooth_Basis{
public:
        Incompress_Solver_Smooth_3D_Basis(Front &front):
	Incompress_Solver_Smooth_Basis(front) {};
	virtual ~Incompress_Solver_Smooth_3D_Basis() {};
protected:
	double getSmoothingFunction(double r);
	double getSmoothingFunctionD(double*, double*);
	double smoothedDeltaFunction(double*, double*);
	double smoothedStepFunction(double*, double*, int);
	void sampleVelocity();
	void setSmoothedProperties(void);
    void addImmersedForce();
};

class Incompress_Solver_Smooth_2D_Cartesian
    : public 	Incompress_Solver_Smooth_2D_Basis
{
public:
        Incompress_Solver_Smooth_2D_Cartesian(Front &front)
            : Incompress_Solver_Smooth_2D_Basis(front)
        {}
	
        ~Incompress_Solver_Smooth_2D_Cartesian()
        {}

	void setInitialCondition(void);
	void setParallelVelocity(void);
	void solve(double dt);
        void vtk_plot_scalar(char*, const char*);
protected:
    
    void computeAdvection(void);
	
    void computeDiffusion(void);
	void computeDiffusionCN(void);
	void computeDiffusionExplicit(void);
	void computeDiffusionImplicit(void);
	void computeDiffusionParab(void);
	
    void computeProjection(void);
	void computeProjectionCim(void);
	void computeProjectionSimple(void);
	void computeProjectionDouble(void);
	    //void computeProjectionDual(void);
	
    void computePressure(void);
	void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computePressureSimple(void);
	
    void computeGradientQ(void);
	void computeNewVelocity(void);
	void extractFlowThroughVelocity(void);
	void computeSourceTerm(double *coords, double *source);
	void surfaceTension(double*, HYPER_SURF_ELEMENT*, HYPER_SURF*, 
				double*, double);
	void computeVarIncrement(double*,double*,boolean);
	void computeVelDivergence();

    void copyMeshStates(void);
	
    /***************   Low level computation functions  *************/
	double getVorticity(int i, int j);
};

class Incompress_Solver_Smooth_3D_Cartesian:
public 	Incompress_Solver_Smooth_3D_Basis{
public:
        Incompress_Solver_Smooth_3D_Cartesian(Front &front):
	Incompress_Solver_Smooth_3D_Basis(front) {};
	~Incompress_Solver_Smooth_3D_Cartesian() {};

	void setInitialCondition(void);
	void setParallelVelocity(void);
	void solve(double dt);
	void solveTest(const char *msg);
        void vtk_plot_scalar(char*, const char*);

protected:
	
    void copyMeshStates(void);
	
    void computeAdvection(void);
	
    void computeDiffusion(void);
	void computeDiffusionCN(void);
	void computeDiffusionExplicit(void);
	void computeDiffusionImplicit(void);
	void computeDiffusionParab(void);
	
    void computeProjection(void);
	void computeProjectionCim(void);
	void computeProjectionSimple(void);
	void computeProjectionDouble(void);
	
    void computePressure(void);
    void computePressurePmI(void);
	void computePressurePmII(void);
	void computePressurePmIII(void);
	void computePressureSimple(void);

	void computeGradientQ(void);
	    //void computeGradientPhi();
	void computeNewVelocity(void);
	void updateComponent(void);
	void extractFlowThroughVelocity(void);
	void computeSourceTerm(double *coords, double *source);
	void surfaceTension(double*, HYPER_SURF_ELEMENT*, HYPER_SURF*, 
				double*, double);
	void computeVarIncrement(double*,double*,boolean);
	void computeVelDivergence();
	void appendOpenEndStates();

    void computeVorticity();
    std::vector<double> computePointVorticity(int* icoords, double** vel);
};


extern int next_index_in_dir(int* icoords, GRID_DIRECTION dir, int dim, int* top_gmax);

extern double getStateVort(POINTER);
extern double getStateXvel(POINTER);
extern double getStateYvel(POINTER);
extern double getStateZvel(POINTER);
extern double getStateXimp(POINTER);
extern double getStateYimp(POINTER);
extern double getStateZimp(POINTER);
extern double getStateComp(POINTER);
extern double getStateMu(POINTER);
extern double getStateDens(POINTER);
extern double getStateTemp(POINTER);

extern double getStatePres(POINTER);
extern double getStatePhi(POINTER);
extern double getStateQ(POINTER);
extern double getStateGradPhiX(POINTER);
extern double getStateGradPhiY(POINTER);
extern double getStateGradPhiZ(POINTER);

extern double getPressure(Front*,double*,double*);
extern double getPhiFromPres(Front* front, double pres);
extern double getQFromPres(Front* front, double pres);

extern double burger_flux(double,double,double);
extern double linear_flux(double,double,double,double);

extern void fluid_print_front_states(FILE*,Front*);
extern void fluid_read_front_states(FILE*,Front*);

extern void read_iF_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern boolean isDirichletPresetBdry(Front*,int*,GRID_DIRECTION,COMPONENT);

extern int ifluid_find_state_at_crossing(Front*,int*,GRID_DIRECTION,
			int,POINTER*,HYPER_SURF**,double*);
extern int ifluid_find_state_at_cg_crossing(Front*,int*,GRID_DIRECTION,
			int,POINTER*,HYPER_SURF**,double*);
extern int ifluid_find_state_at_dual_crossing(Front*,int*,GRID_DIRECTION,
			int,POINTER*,HYPER_SURF**,double*);

extern double p_jump(POINTER,int,double*);
extern double grad_p_jump_n(POINTER,int,double*,double*);
extern double grad_p_jump_t(POINTER,int,int,double*,double*);

extern boolean neumann_type_bdry(int);

/*extern void setSlipBoundary(int* icoords, int idir, int nb, int comp,
        HYPER_SURF* hs, POINTER state, double** vel, double* vtmp);*/

extern void ifluid_compute_force_and_torque(Front*,HYPER_SURF*,double,double*,
                        double*);

extern void initInnerBoundary(Front*,LEVEL_FUNC_PACK*);
extern void restart_set_dirichlet_bdry_function(Front*);

extern void iF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void iF_timeDependBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);

extern void ifluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void ifluid_compute_force_and_torque(Front*,CURVE*,double,double*,
                        double*);

extern void setInitialIntfc(Front*,LEVEL_FUNC_PACK*,char*,IF_PROB_TYPE);
extern void init_fluid_state_func(Incompress_Solver_Smooth_Basis*,IF_PROB_TYPE);
extern void read_iFparams(char*,IF_PARAMS*);
extern void read_iF_prob_type(char*,IF_PROB_TYPE*);
extern void recordBdryEnergyFlux(Front*,char*);

extern void read_open_end_bdry_data(char*,Front*);
extern void setContactNodeType(Front*);
extern int contact_node_propagate(Front*,POINTER,NODE*,NODE*,RPROBLEM**,
                               double,double*,NODE_FLAG);
#endif
