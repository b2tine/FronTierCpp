#ifndef CFLUID_H
#define CFLUID_H

/**********************************************************************
 * 		cFluid.h
 **********************************************************************/

#include <FronTier.h>
#include "state.h"
#include "rigidbody.h"
#include "cFvisc.h"

#include <vector>
#include <string>
#include <list>
#include <assert.h>

#define         EXT_COMP		0
#define         SOLID_COMP		0
#define         GAS_COMP1		2
#define         GAS_COMP2		3
#define         GAS_COMP3		4
#define         MAX_COMP        10

#define	gas_comp(comp) (((comp) == GAS_COMP1 || \
            (comp) == GAS_COMP2 || (comp) == GAS_COMP3) ? YES : NO)

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
    NACA0012_AIRFOIL,
    CHANNEL_FLOW,
    NO_FLUID
};

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

//TODO: Move to cAirfoil
enum class PORO_SCHEME
{
    NORMAL_REFLECTION,
    REFLECTION,
    RIEMANN,
    DARCY
};

//TODO:Add function addImmersedForce() for drag on parachute strings
//
// For fluid drag on string chord
struct FINITE_STRING 
{
    double radius;
    double dens;
    double c_drag;
    double ampFluidFactor;
};

enum class EDDY_VISC_MODEL
{
    NONE,
    VREMAN,
    WALE
};

//For specifying freestream condition to initialize flow field,
//and for external flow problems to apply far-field (outflow) boundary conditions.
struct FREESTREAM_PARAMS
{
    double M;   //mach number
    double P;   //pressure
    double T;   //temperature
    double Re;  //reynolds number
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
    double shock_speed;
	
    double T0;
    double T1;
    double T2;
    double T3;

    double p0;
    double p1;
    double p2;
    double p3;

    double p_operating; //environment reference pressure
	
    double rho0;
    double rho1;
    double rho2;
    double rho3;
	
    double v0[MAXD];
    double v1[MAXD];
    double v2[MAXD];
    double v3[MAXD];
	
    double mu1;
	double mu2;
	
    double gamma;
        //double R_specific;

	double Mach_number;
	
    double gravity[MAXD];
    
    double min_dens;
    double min_pres;

	double shock_position;
	double contact_vel;
    
	double **vel;
	double *vort;    /* Vorticity-2d */
    double **vorticity; /* Vorticity-3d */

    double** f_surf;

	double **mom;
	double *dens;
	double *engy;
	double *pres;
    double *mu;
    double *mu_turb;
    double *k_turb;
    double *temp;   //temperature

    //GFM
	double **gnor;
	double **Gdens;
	double ***Gvel;
	double **Gpres;

    bool is_inviscid {false}; //euler eqns solution

	//LES Turbulence
    bool use_eddy_viscosity {false};
    EDDY_VISC_MODEL eddy_viscosity_model;
    double C_v; //Vreman model coefficient

    bool perturb_const_inlet_bdry {false};

    ////////////////////////////////////////////////////////////////////////////
    //TODO: Package all these into the FREESTREAM_PARAMS struct
    //
    //FREESTREAM VARIABLES
    int dir_freestream;     //freestream orientation ({x,y,z} = {0,1,2})
    double alpha {0.0};     //angle of attack wrt to the freestream orientation
    double beta {0.0};      //sideslip angle wrt to the freestream orientation
    
    double Mach_freestream; //freestream mach number
    double c_freestream;    //freestream sound speed
    double dens_freestream; //freestream density
    double pres_freestream; //freestream pressure
    double vel_freestream;  //freestream velocity
    ////////////////////////////////////////////////////////////////////////////

    bool no_slip_wall {false};
    bool use_fixed_wall_temp {false};
    double fixed_wall_temp {294.0};

    bool use_preset_dist_reflect {false};
    double dist_reflect;


    /////////////////////////////    
    bool no_fluid {false};
    int fsi_startstep;
    
    boolean with_porosity;
    double porosity;
    PORO_SCHEME poro_scheme;
    double porous_coeff[2];
    double permeability {7.6033e-12};

    bool with_string_fsi {false};
    //FINITE_STRING* string_fsi_params;
    double string_radius {0.0};
    double string_density {0.0};
    double string_drag_coef {0.0};
    double string_ampFactor {0.0};
    //////////////////////////////
    
    //Base front for comparison
	boolean use_base_soln;
    char base_dir_name[200];
    int num_step;
    int *steps;
    F_BASIC_DATA *f_basic;
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

struct WHITE_NOISE_PARAMS
{
    int dim;
    double amplitude;
    STATE base_state;
};

//TODO: This should be named PISTON_BDRY_PARAMS or similar
//      -- not generic variable bdry
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

/*********************************************************************
*	Definition and structures for Riemann solution.              *
**********************************************************************/

#define         NO_WAVE                 1
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
    int wave_type;                      /* simple wave or shock */
    double  speed_leading_edge;         /* for simple wave */
    double speed_trailing_edge;
    double  speed_shock;                /* for shock */
    double  speed_contact;              /* for contact line */
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
/*********************************************************************/


/*
enum VISITED_TYPE 
{
    UNVISITED,
    VISITED,
    PARTIAL_VISITED,
    FULL_VISITED
};
*/

struct FIELD
{
    double **f_surf;
	double **vel;
	double **momn;
	double *dens;
	double *engy;
	double *pres;
    double *vort;
    double **vorticity;
    double *mu;
    double *mu_turb;
    double *k_turb;
    double *temp;
};

struct SWEEP
{
    double *dens;           /* density vector */
    double **momn;          /* momentum vector */
    double **vel;
    double *engy;           /* internal energy vector */
    double *pres;           /* used for EOS */
    double *mu;
    double *mu_turb;
    double *k_turb;
    double *temp;
};

struct FSWEEP
{
    double *dens_flux;      /* density flux */
    double **momn_flux;     /* momentum flux */
    double *engy_flux;      /* internal energy flux */
};

class L_RECTANGLE 
{
public:
	int comp {-1};
	//int old_comp {-1};
	int m_index {-1};
	double m_coords[MAXD];
	int icoords[MAXD];

public:

	void setCoords(double* coords, int dim)
    {
        for (int i = 0; i < dim; ++i)
            m_coords[i] = coords[i];
    }

    void getCoords(double* coords)
    {
        for (int i = 0; i < MAXD; ++i)
            coords[i] = m_coords[i];
    }

    std::vector<double> getCoords()
    {
        return std::vector<double>(m_coords,m_coords+3);
    }
};

class G_CARTESIAN
{
public:
	
    EQN_PARAMS* eqn_params;

	explicit G_CARTESIAN(Front* ft)
        : front{ft}
    {
        //eqn_params = (EQN_PARAMS*)front->extra1;
    }
	
    virtual ~G_CARTESIAN() = default;
	
    int dim;
	double m_dt;			// time increment
	double max_dt;			// max_dt from cartesian (advection)
    double visc_max_dt;     
    double min_dt;
	double hmin;			// smallest spacing

    void setInitialIntfc(LEVEL_FUNC_PACK*,char*);// setup initial geometry

	void initMesh(void);		// setup the cartesian grid
	void setInitialStates(); 	// setup initial state

	//void readFrontInteriorStates(char *restart_state_name); //TODO: write this

    int (*findStateAtCrossing)(Front*,int*,GRID_DIRECTION,int,
            POINTER*,HYPER_SURF**,double*);

    void readInteriorStates(char*);
	void printFrontInteriorStates(char*);
	void initMovieVariables();
	void getVelocity(double*,double*);
	void initSampleVelocity(char *in_name);
	void compareWithBaseData(char *out_name);
	void freeBaseFront();
	void errFunction();


    void writeMeshFileVTK();
    void writeCompGridMeshFileVTK();
    void writeMeshComponentsVTK();  //parallel version not implemented yet; calls serial version
    
    void serialWriteMeshFileVTK();          //should be safe to remove
    void serialWriteCompGridMeshFileVTK();  //should be safe to remove
    void serialWriteMeshComponentsVTK();
    

	void solve(double dt);

//private:
protected:
	
    Front *front;
	
    // On topological grid
	RECT_GRID *top_grid;
	double *array;		// for scatter states;
	COMPONENT *top_comp;
	    //EQN_PARAMS *eqn_params;
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
	double m_dens[2];		    // two component at most
	double m_smoothing_radius;	// used by getSmoothingFunction()

	double m_t;         // time
	double max_speed;	// for stability of convection

    //for viscous flux time step restriction
    double mu_max {-1};
	double rho_min {-1};

    //User defined minimum physical variables
    double min_dens;
    double min_pres;

	// for parallel partition
	int NLblocks, ilower, iupper;
    int *n_dist;

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
	//
	void scatMeshStates();
	void scatMeshVst(SWEEP*);
	void scatMeshFlux(FSWEEP*);
	void appendOpenEndStates();

	// -------------------------------------------------------
	// 		compressible solver functions
	// -------------------------------------------------------
	void setMaxTimestep(void);
	void advanceSolution(void);

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
    
    void computeMeshFlux(SWEEP&,FSWEEP*,double); //void computeMeshFlux(SWEEP,FSWEEP*,double);

	void copyMeshVst(const SWEEP&,SWEEP*);
	void copyFromMeshVst(const SWEEP&);
	void copyToMeshVst(SWEEP*);

	void addSourceTerm(SWEEP*,FSWEEP*,double); //void addSourceTerm(const SWEEP&,FSWEEP*,double);
	
    void addImmersedForce(); //for string fluid interaction, and potentially other forces

	/* Directional flux solver */
	void resetFlux(FSWEEP*);
	void addFluxInDirection(int,SWEEP*,FSWEEP*,double);
	
    virtual void addFluxAlongGridLine(int,int*,double,SWEEP*,FSWEEP*);
	
    /* Advective flux */
	void numericalFlux(POINTER,SWEEP*,FSWEEP*,int);
    
    virtual void appendGhostBuffer(SWEEP*,SWEEP*,int,int*,int,int);

	void appendStencilBuffer2d(SWEEP*,SWEEP*,int,int);//UNUSED
	void appendStencilBuffer3d(SWEEP*,SWEEP*,int,int,int);//UNUSED
	
	void setDirichletStates(STATE*,SWEEP*,SWEEP*,HYPER_SURF*,int*,int,int,int,int,int);

	void setNeumannStates(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,int*,int,int,int,int,int);
	void setNeumannStatesNEW(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,int*,int,int,int,int,int);
	void setNeumannStatesOLD(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,int*,int,int,int,int,int);
	
	/* Viscous flux */
    void addViscousFlux(SWEEP* m_vst, FSWEEP* m_flux, double delta_t);
    
    void fillViscousFluxStencil2d(int* icoords, SWEEP* m_vst, VStencil2d* vsten);
    void fillViscousFluxStencil2d_5pt(int* icoords, SWEEP* m_vst, VStencil2d_5pt* vsten);

    void fillViscousFluxStencil3d(int* icoords, SWEEP* m_vst, VStencil3d* vsten);
    void fillViscousFluxStencil3dNEW(int* icoords, SWEEP* m_vst, VStencil3d* vsten);
        //void fillViscousFluxStencil3d_5pt(int* icoords, SWEEP* m_vst, VStencil3d_5pt* vsten);
    
    virtual void setViscousGhostState(int* icoords, COMPONENT comp, VSWEEP* vs, SWEEP* m_vst);
    
    void setDirichletViscousGhostState(VSWEEP* vs, COMPONENT comp, double* intrp_coeffs,
            HYPER_SURF_ELEMENT* hse, HYPER_SURF* hs);

    void setNeumannViscousGhostState(int* iccords, SWEEP* m_vst, VSWEEP* vs, double* ghost_coords,
            double* crx_coords, COMPONENT comp, double* intrp_coeffs,
            HYPER_SURF_ELEMENT* hse, HYPER_SURF* hs);
    
    void computeViscousFlux2d(int* icoords, SWEEP* m_vst, VFLUX* v_flux,
            double delta_t, VStencil2d* vsten);
    void computeViscousFlux2d_5pt(int* icoords, SWEEP* m_vst, VFLUX* v_flux,
            double delta_t, VStencil2d_5pt* vsten);

    void computeViscousFlux3d(int* icoords, SWEEP* m_vst, VFLUX* v_flux,
            double delta_t, VStencil3d* vsten);
    

    //Sutherland's Law
    void computeDynamicViscosity();
    void computeDynamicViscosity2d();
    void computeDynamicViscosity3d();

    //For LES turbulence
    void computeSGSTerms();

    void computeEddyViscosity();
    void computeEddyViscosity2d();
    void computeEddyViscosity3d();
    
    double computeEddyViscosityVremanModel(int *icoords);
    double computeEddyViscosityVremanModel_BdryAware(int *icoords);

    double computeEddyViscosityWALE(int *icoords);
    
    std::vector<std::vector<double>> computeVelocityGradient(int *icoords);
    
    void setSlipBoundary(int* icoords, int idir, int nb, int comp,
            HYPER_SURF* hs, POINTER state, double** vel, double* v_slip);

    void setSlipBoundaryNIP(int* icoords, int idir, int nb, int comp,
            HYPER_SURF* hs, POINTER state, double** vel, double* v_slip);

    void setPoroSlipBoundaryNIP(int* icoords, int idir, int nb, int comp,
            HYPER_SURF* hs, POINTER state, double** vel, double* v_slip);


    // -------------------------------------------------------
	// 		initialization functions
	// -------------------------------------------------------
	void initChannelFlow(LEVEL_FUNC_PACK*,char*);
	void initSinePertIntfc(LEVEL_FUNC_PACK*,char*);
	void initRandPertIntfc(LEVEL_FUNC_PACK*,char*);
	void initCirclePlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initImplosionIntfc(LEVEL_FUNC_PACK*,char*);
	void initProjectileIntfc(LEVEL_FUNC_PACK*,char*);
	void initProjectileIntfc2d(LEVEL_FUNC_PACK*,char*);
	void initMTFusionIntfc(LEVEL_FUNC_PACK*,char*);
	void initRiemannProb(LEVEL_FUNC_PACK*,char*);

	void initChannelFlowStates();
	void initRayleiTaylorStates();
	void initRichtmyerMeshkovStates();
	void initBubbleStates();
	void initImplosionStates();
	void initMTFusionStates();
	void initRiemProbStates();
	void initBlastWaveStates();
	void initShockSineWaveStates();
	void initAccuracySineWaveStates();
	
    void initNACA0012PlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initRectPlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initTrianglePlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initCylinderPlaneIntfc(LEVEL_FUNC_PACK*,char*);
	void initObliqueIntfc(LEVEL_FUNC_PACK*,char*);
	void initObliqueStates();
	
    void readBaseFront(int i);
	void readBaseStates(char *restart_state_name);

	void compSGS(void); //UNUSED FUNCTION (subgrid model by Hyunkyung Lim)

	void getPressureJumpParameter(double *coords0, double *coords1, 
			double &theta, double &jumpPressure, 
			double &jumpDerivative);//NO DEFINITION

	// velocity field query
	void getVelocityGradient(double *p,double *gradU,double *gradV);//NO DEFINITION

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

    void computeVorticity();
	double getVorticity(int i, int j);
    double getVorticityX(int i, int j, int k);
    double getVorticityY(int i, int j, int k);
    double getVorticityZ(int i, int j, int k);
	
    double getDistance(double *coords0, double *coords1);

	void getNearestInterfacePoint(double *q,double *p); // incompletely implemented
		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	void save(char *filename);


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


// cFinit.cpp
extern void read_cFluid_params(char*,EQN_PARAMS*);
extern void getChannelInletState(STATE* state, EQN_PARAMS *eqn_params, COMPONENT comp);
extern void insert_objects(Front*);
extern void insert_boundary_objects(Front* front);
extern void initSplitState(Front* front);
extern void initBackwardFacingStep(Front* front);
extern void initBump(Front* front);
extern void initRamp(Front* front);


// cFcartsn.cpp
extern void readFrontStates(Front*,char*);

// cFsub.cpp
extern double burger_flux(double,double,double);
extern double linear_flux(double,double,double,double);
extern void read_dirichlet_bdry_data(char*,Front*);
extern void read_open_end_bdry_data(char*,Front*);
extern void restart_set_dirichlet_bdry_function(Front*);

extern void cF_constantWithWhiteNoise(double*,HYPER_SURF*,Front*,POINTER,POINTER);
extern void cF_variableBoundaryState(double*,HYPER_SURF*,Front*,POINTER,POINTER);

extern void cF_farfieldBoundaryState(double*,HYPER_SURF*,Front*,POINTER,POINTER);
extern void cF_supersonicOutflowState(double*,HYPER_SURF*,Front*,POINTER,POINTER);
extern void cF_symmetryBoundaryState(double*,HYPER_SURF*,Front*,POINTER,POINTER);

extern void cF_flowThroughBoundaryState(double*,HYPER_SURF*,Front*,POINTER,POINTER);
extern void cF_flowThroughBoundaryState2d(double*,HYPER_SURF*,Front*,POINTER,POINTER);
extern void cF_flowThroughBoundaryState3d(double*,HYPER_SURF*,Front*,POINTER,POINTER);

extern void cFluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);

extern boolean reflectNeumannState(Front*,HYPER_SURF*,double*,COMPONENT,SWEEP*,STATE*);
extern void reflectVectorThroughPlane(double*,double*,double*,int);

extern void findGhostState(STATE,STATE,STATE*);

// cFeos.cpp
extern double EosViscosity(STATE*);
extern double EosTemperature(STATE*);
extern double EosDensity(STATE*);
extern double EosPressure(STATE*);
extern double EosInternalEnergy(STATE*);
extern double EosEnergy(STATE*);
extern double EosSoundSpeed(STATE*);
extern double EosSoundSpeedSqr(STATE*);
extern double EosMaxBehindShockPres(double,STATE*);
extern void EosSetTVDParams(SCHEME_PARAMS*,EOS_PARAMS*);
extern void ConvertVstToState(STATE*,SWEEP*,EOS_PARAMS*,int,int);

/* Riemann solution functions */
// cFriem.cpp
extern boolean RiemannSolution(RIEMANN_INPUT,RIEMANN_SOLN*);
extern boolean RiemannSolnAtXi(RIEMANN_SOLN*,RIEM_STATE*,double);

/* Structures and functions for TVD scheme */
extern void TVD_flux(POINTER,SWEEP*,FSWEEP*,int);
extern void WENO_flux(POINTER,SWEEP*,FSWEEP*,int);

/*	rgbody.c functions */
//TODO: Move the above function declarations and definitions into rigidbody.h/.cpp
extern void cfluid_compute_force_and_torque(Front*,HYPER_SURF*,double,double*,double*);
extern void record_moving_body_data(char*,Front*);//NO DEFINITION


//Unused?
class EOS
{
public:
	
    EOS(EOS_PARAMS &params);
	double Pressure(STATE);
	double Energy(STATE);

private:

	EOS_PARAMS *params;
};



#endif
