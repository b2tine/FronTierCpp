#ifndef CFLUID_H
#define CFLUID_H

#include <FronTier.h>
#include "state.h"

#include <vector>
#include <list>
#include <assert.h>

#define EXT_COMP		0
#define SOLID_COMP		0
#define GAS_COMP1		2
#define GAS_COMP2		3
#define MAX_COMP       10

#define	gas_comp(comp) (((comp) == GAS_COMP1 || \
            comp == GAS_COMP2) ? YES : NO)

enum PROB_TYPE
{
        ERROR_TYPE = -1,
        TWO_FLUID_BUBBLE = 1,
        TWO_FLUID_RT,
        FLUID_SOLID_CIRCLE,
        BUBBLE_SURFACE,
        FLUID_RIGID_BODY,
	ROTOR_ONE_FLUID,
	ROTOR_TWO_FLUID,
        TWO_FLUID_RM,
        TWO_FLUID_RM_RAND,
        IMPLOSION,
        MT_FUSION,
        PROJECTILE,
        RIEMANN_PROB,
	FLUID_CRYSTAL,
	ONED_BLAST,
	ONED_SSINE,
	ONED_ASINE,
        FLUID_SOLID_RECT,
        FLUID_SOLID_TRIANGLE,
        FLUID_SOLID_CYLINDER,
        OBLIQUE_SHOCK_REFLECT,
    FABRIC
};

struct EOS_PARAMS
{
    double  gamma;
    double  pinf;
    double  einf;
};

//Now in state.h
//
//struct STATE {
//	double dens;			/* density */
//	double engy;			/* energy density */
//	double momn[MAXD];		/* momentum deisnty */
//	double pres;			/* Pressure */
//	double vel[MAXD];		/* Velocities */
//	double vort;			/* Vorticity */
//	EOS_PARAMS *eos;
//	int dim;
//};

enum NUM_SCHEME
{
	TVD_FIRST_ORDER		=	1,
        TVD_SECOND_ORDER,
        TVD_FOURTH_ORDER,
        WENO_FIRST_ORDER,
        WENO_SECOND_ORDER,
        WENO_FOURTH_ORDER
};

enum POINT_PROP_SCHEME
{
	FIRST_ORDER		=	1,
    SECOND_ORDER,
    FOURTH_ORDER
};

enum SHOCK_PARAMETER
{
	BEHIND_PRESSURE 	=	1,
	BEHIND_VELOCITY,
	SHOCK_SPEED,
	SHOCK_MACH_NUMBER
};

struct EQN_PARAMS
{
    int dim;
    PROB_TYPE prob_type;
    POINTER level_func_params;
	NUM_SCHEME num_scheme;
    POINT_PROP_SCHEME point_prop_scheme;
	EOS_PARAMS eos[MAX_COMP];
	boolean tracked;
	boolean articomp;
	boolean contact_stationary;
	int idir;
	int shock_side;
	
    double p0;
        double p1;
        double p2;
        double p3;
	double rho0;
        double rho1;
        double rho2;
        double rho3;
	double v0[MAXD];
        double v1[MAXD];
        double v2[MAXD];
        double v3[MAXD];

    double min_dens;
    double min_pres;
	double mu1;
	double mu2;
	double gamma;
	double gravity[MAXD];
	double Mach_number;
	double shock_position;
	double contact_vel;
	double **vel;
	double **mom;
	double *dens;
	double *engy;
	double *pres;
	double *vort;       //2d
	double **vorticity; //3d

    //GFM
	double **gnor;
	double **Gdens;
	double ***Gvel;
	double **Gpres;

	// Base front for comparison
	boolean use_base_soln;
    char base_dir_name[200];
    int num_step;
    int *steps;
    F_BASIC_DATA *f_basic;

    boolean with_porosity;
    double porosity;
};

struct SCHEME_PARAMS
{
	boolean artificial_compression;
    double lambda;
    double beta;
	double gamma, einf;
};

struct FLOW_THROUGH_PARAMS
{
    POINT *oldp;
    COMPONENT comp;
    EQN_PARAMS *eqn_params;
};

//Now in rigid_body.h
//
//struct RG_PARAMS
//{
//    int dim;
//	boolean no_fluid;               /* For benchmark tests */
//        double  total_mass;             /* Total mass */
//        double  moment_of_inertial;     /* Moment of inertial about the axis */
//        double  center_of_mass[MAXD];   /* Center of mass */
//        double  rotation_dir[MAXD];     /* Direction of rotation */
//        double  translation_dir[MAXD];  /* Restricted direction of motion */
//        double  rotation_cen[MAXD];     /* Center of rotation */
//        double  cen_of_mass_velo[MAXD]; /* Center of mass velocity */
//        double  angular_velo;           /* Angular velocity of rotation */
//        double  p_moment_of_inertial[MAXD];
//        double  p_angular_velo[MAXD];
//        double  euler_params[4];
//        double  old_euler_params[4];
//        MOTION_TYPE motion_type;
//};


struct VAR_BDRY_PARAMS
{
	int dim;
    double center[MAXD];        /* Center of disk/sphere */
    double *angles_pistons;     /* Angles to the pistons' centers */
    double half_angular_width; /* Half-angle of the piston's surface */
    double bdry_vel;            /* Boundary velocity */
    double bdry_dens;           /* Boundary density */
    double bdry_pres;           /* Boundary pressure */
    int number_pistons;         /* Number of pistons */
	double jet_duration_time;   /* Time duration for the jet */
};

/*********************************************************************
*	Definition and structures for Riemann solution.              *
**********************************************************************/

#define         GAMMA_PLUS              2
#define         GAMMA_MINUS             3
#define         LF_SHOCK                4
#define         RF_SHOCK                5
#define         CONTACT                 6
#define         VACUUM                  7

struct RIEM_STATE
{
    double d,p,u;
    double gamma;
};

struct CENTERED_WAVE
{
    int wave_type;                   /* simple wave or shock */
    double speed_leading_edge;       /* for simple wave */
    double speed_trailing_edge;      
    double  speed_shock;             /* for shock */
    double  speed_contact;           /* for contact line */
};

struct RIEMANN_SOLN
{
    RIEM_STATE left_state;
    RIEM_STATE right_state;
    RIEM_STATE left_center_state;
    RIEM_STATE right_center_state;
    CENTERED_WAVE left_wave;
    CENTERED_WAVE right_wave;
    CENTERED_WAVE contact;
};

struct RIEMANN_INPUT
{
    RIEM_STATE left_state;
    RIEM_STATE right_state;
};

struct OPEN_PIPE_PARAMS
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

struct L_RECTANGLE
{
    int comp {-1};			 
	int m_index {-1};
	double m_coords[MAXD];	
	int icoords[MAXD];

	void setCoords(double* coords, int dim)
    {
        for (int i = 0; i < dim; ++i)
            m_coords[i] = coords[i];
    }
};

struct FIELD
{
	double **vel;
	double **momn;
	double *dens;
	double *engy;
	double *pres;
	double *vort;       //2d
	double **vorticity; //3d
};

struct SWEEP
{
    double *dens;       /* density vector */
    double **momn;      /* momentum vector */
    double *engy;       /* internal energy vector */
    double *pres;       /* used for EOS */
};

struct FSWEEP
{
    double *dens_flux;  /* density flux */
    double **momn_flux; /* momentum flux */
    double *engy_flux;  /* internal energy flux */
};


class G_CARTESIAN
{
public:

    G_CARTESIAN(Front* ft)
        : front{ft}
    {}

    //~G_CARTESIAN() = default;

    int dim;
    double m_dt;			// time increment
    double max_dt;			// max_dt from cartesian
    double hmin;			// smallest spacing

    void setInitialIntfc(LEVEL_FUNC_PACK*,char*);// setup initial geometry
    void setInitialStates(); 	// setup initial state
    //void setProbParams(); 	// setup initial state
    void initMesh();		// setup the cartesian grid
    void readGasStates(char* restart_state_name);
    void readFrontStates(Front* frt, char* restart_state_name);
    void readInteriorStates(char* restart_state_name);
    void printGasStates(char* out_name);
    void initMovieVariables();
    void getVelocity(double*,double*);
    void initSampleVelocity(char *in_name);
    void compareWithBaseData(char *out_name);
    void freeBaseFront();
    void errFunction();

    //cFbase.cpp
    void applicationSetComponent();
    void applicationSetStates();

    // main step function
    void solve(double dt);		

private:

	Front* front;

	// On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	COMPONENT *top_comp;
	EQN_PARAMS *eqn_params;
	FIELD field;
	FIELD *base_field;
	Front *base_front;

	int top_gmax[MAXD];
	int lbuf[MAXD],ubuf[MAXD];
	double top_L[MAXD],top_U[MAXD],top_h[MAXD];
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D
	int nrad;			// Buffer size for a given solver

	// Sweeping limits
	int imin[MAXD];
	int imax[MAXD];

	// member data: mesh storage
	std::vector<L_RECTANGLE> cell_center;

	// member data: 
	int m_comp[2];
	double m_mu[2];
	double m_dens[2];		// two component at most
	double m_smoothing_radius;	// used by getSmoothingFunction()

	double m_t;                     // time
	double max_speed;		// for stability of convection
	double min_dens,min_pres;	// minimum physical variables

	// for parallel partition
	int NLblocks,ilower,iupper;
    int* n_dist;

	// mesh: full cells mesh
	void setComponent(void);	// init components	
	void setDomain();
	void augmentMovieVariables(void);
	void copyMeshStates();
	void sampleVelocity();
	void sampleVelocity2d();
	void sampleVelocity3d();

	/*TMP*/
	void checkVst(SWEEP*);
	void checkFlux(FSWEEP*);

	// parallelization related functions
	void scatMeshStates();
	void scatMeshVst(SWEEP*);
	void scatMeshFlux(FSWEEP*);
	void appendOpenEndStates();

	// compressible solver functions
	void setAdvectionDt(void);
	void computeAdvection(void);

	/* Mesh memory management */
	bool withinStencilLen(int*,int);
	void allocMeshVst(SWEEP*);
	void allocMeshFlux(FSWEEP*);
	void allocDirVstFlux(SWEEP*,FSWEEP*);
	void freeDirVstFlux(SWEEP*,FSWEEP*);
	void freeVst(SWEEP*);
	void freeFlux(FSWEEP*);

	/* Mesh operations */
	void solveRungeKutta(int);
	void addMeshFluxToVst(SWEEP*,const FSWEEP&,double);
	void computeMeshFlux(SWEEP,FSWEEP*,double);
	void copyMeshVst(const SWEEP&,SWEEP*);
	void copyFromMeshVst(const SWEEP&);
	void copyToMeshVst(SWEEP*);
	void addSourceTerm(const SWEEP&,FSWEEP*,double);

	/* Directional flux solver */
	void resetFlux(FSWEEP*);
	void addFluxInDirection(int,SWEEP*,FSWEEP*,double);
	void addFluxAlongGridLine(int,int*,double,SWEEP*,FSWEEP*);
	void augmentOneDimBuffer(int,int);
	void numericalFlux(POINTER,SWEEP*,FSWEEP*,int);
	void appendStencilBuffer2d(SWEEP*,SWEEP*,int,int);
	void appendStencilBuffer3d(SWEEP*,SWEEP*,int,int,int);
	void appendGhostBuffer(SWEEP*,SWEEP*,int,int*,int,int);
	
    void initFabricStates();

    //TODO: move out of class
    // initialization functions
	
    void initSinePertIntfc(LEVEL_FUNC_PACK*,char*);
	void initRandPertIntfc(LEVEL_FUNC_PACK*,char*);
	void initCirclePlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initImplosionIntfc(LEVEL_FUNC_PACK*,char*);
	void initProjectileIntfc(LEVEL_FUNC_PACK*,char*);
	void initProjectileIntfc2d(LEVEL_FUNC_PACK*,char*);
	void initMTFusionIntfc(LEVEL_FUNC_PACK*,char*);
	void initRiemannProb(LEVEL_FUNC_PACK*,char*);
	void initRayleiTaylorStates();
	void initRichtmyerMeshkovStates();
	void initBubbleStates();
	void initImplosionStates();
	void initMTFusionStates();
	void initProjectileStates();
	void initRiemProbStates();
	void initBlastWaveStates();
	void initShockSineWaveStates();
	void initAccuracySineWaveStates();
	void initRectPlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initTrianglePlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initCylinderPlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initObliqueIntfc(LEVEL_FUNC_PACK*,char*);
	void initObliqueStates();
	
	void readBaseFront(int i);
	void readBaseStates(char *restart_name);

	void compSGS(void);             // subgrid model by Hyunkyung Lim

	void getPressureJumpParameter(
            double *coords0, double *coords1, 
			double &theta, double &jumpPressure, 
			double &jumpDerivative);

	// velocity field query
	void getVelocityGradient(double *p,double *gradU,double *gradV);
    void computeVorticity();

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	int getRectangleComponent(int index);	// the center component
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	
	int getInteger(double i);
	boolean isInteger(double i);

	double getVorticity(int i, int j);
    double getVorticityX(int i, int j, int k);
    double getVorticityY(int i, int j, int k);
    double getVorticityZ(int i, int j, int k);
	double getDistance(double *coords0, double *coords1);

			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	void save(char *filename);

	void setDirichletStates(STATE*,SWEEP*,SWEEP*,HYPER_SURF*,
                            int*,int, int,int,int);

	void setNeumannStates(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

	void setElasticStates(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

	void setElasticStatesRiem(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

	//GFM
	void solve_exp_value();
	boolean get_ave_normal(int*,int***);
	boolean get_ave_state(const SWEEP&, int*,int***,int,int);
	boolean needBufferFromIntfc(COMPONENT,COMPONENT);
	void get_normal_from_front();
	void get_ghost_state(const SWEEP&, int,int);
	void tecplot_interior_states(char*);
	void scatMeshGhost();
	void GFMGhostState(int*,int,STATE*);
	void checkCorrectForTolerance(STATE*);
	void adjustGFMStates();
};

// cFphys.cpp
extern void read_cFluid_params(char*,EQN_PARAMS*);

// cFsub.cpp
extern double burger_flux(double,double,double);
extern double linear_flux(double,double,double,double);
extern void read_dirichlet_bdry_data(char*,Front*);
extern void read_open_end_bdry_data(char*,Front*);
extern void restart_set_dirichlet_bdry_function(Front*);
extern void reflectVectorThroughPlane(double*,double*,double*,int);
extern boolean reflectNeumannState(Front*,HYPER_SURF*,double*,COMPONENT,SWEEP*,
			STATE*);
extern void cfluid_compute_force_and_torque(Front*,HYPER_SURF*,double,double*,
			double*);
extern void cF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void cF_variableBoundaryState(double*,HYPER_SURF*,Front*,POINTER,
                        POINTER);
extern void cFluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);


/*	rgbody.c functions */
extern void init_moving_bodies(EQN_PARAMS,char*,Front*);
extern void record_moving_body_data(char*,Front*);

// cFeos.cpp
extern double EosPressure(STATE*);
extern double EosInternalEnergy(STATE*);
extern double EosEnergy(STATE*);
extern double EosSoundSpeed(STATE*);
extern double EosSoundSpeedSqr(STATE*);
extern double EosMaxBehindShockPres(double,STATE*);
extern void   EosSetTVDParams(SCHEME_PARAMS*,EOS_PARAMS*);
extern void   CovertVstToState(STATE*,SWEEP*,EOS_PARAMS*,int,int);
extern void   findGhostState(STATE,STATE,STATE*);

	/* Riemann solution functions */
/* cFriem.cpp */
extern boolean RiemannSolution(RIEMANN_INPUT,RIEMANN_SOLN*);
extern boolean RiemannSolnAtXi(RIEMANN_SOLN*,RIEM_STATE*,double);

	/* Structures and functions for TVD scheme */

extern	void TVD_flux(POINTER,SWEEP*,FSWEEP*,int);
extern	void WENO_flux(POINTER,SWEEP*,FSWEEP*,int);

class EOS{
	EOS_PARAMS *params;
public:
	EOS(EOS_PARAMS &params);
	double Pressure(STATE);
	double Energy(STATE);
};

class COMPRESSIBLE_GAS_SOLVER{
	Front *front;
public: 
	COMPRESSIBLE_GAS_SOLVER(Front &front);

	void set_solver_domain(void);
        void runge_kutta(void);
        void solve(void);
        void solve1d(void);
        void solve2d(void);
        void solve3d(void);
private:
	int dim;
        COMPONENT *top_comp;
        double *top_h;
        int *top_gmax;
        int imin,jmin,kmin;
        int imax,jmax,kmax;
};

//extern void prompt_for_rigid_body_params(int,char*,RG_PARAMS*);
//extern void set_rgbody_params(RG_PARAMS,HYPER_SURF*);

// cFinit.cpp
extern void insert_objects(Front*);

void set_cFluid_params(char*, EQN_PARAMS*);
void setFabricParams(char*,EQN_PARAMS*);

void setRayleiTaylorParams(char*,EQN_PARAMS*);
void setRichtmyerMeshkovParams(char*,EQN_PARAMS*);
void setBubbleParams(char*,EQN_PARAMS*);
void setImplosionParams(char*,EQN_PARAMS*);
void setMTFusionParams(char*,EQN_PARAMS*);
void setProjectileParams(char*,EQN_PARAMS*);
void setRiemProbParams(char*,EQN_PARAMS*);
void setRiemProbParams1d(char*,EQN_PARAMS*);
void setRiemProbParams2d(char*,EQN_PARAMS*);
void setOnedParams(char*,EQN_PARAMS*);
//void setObliqueParams(char*,EQN_PARAMS*);

#endif
