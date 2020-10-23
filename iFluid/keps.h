#ifndef KEPS_H
#define KEPS_H

#include "iFluid.h"
/*Following is for turbulence model RNG k-eps model*/
    
struct KE_FIELD {
    double *gamma;//TODO: can probably remove since not the same for k and eps anymore
	double *k;
	double *eps;
	double *Pk;
	double *mu_t;
	double *Cmu;
	double *temp; /*temporary field for debugging, not temperature*/
	double **vel;
};

struct KE_PARAMS {
    int dim;
	double delta_k;
	double delta_eps;
	double C2;
	double C1;
	double B; /*constant in log law, 5.2 for smooth wall*/
	double E; /*alt constant in log law, 9.0 for smooth wall*/
	double Cbc;
	double Cmu;
	double mu0;
	double mu;
	double rho;
	double l0;      //turbulence length scale
	double I;       //turbulence intensity
    double ViscRatioFar;      //eddy viscosity ratio far-field
    double ViscRatioNear;      //eddy viscosity ratio near-field
    double U0;      //freestream average velocity
	double lmax;
	double delta;   //specify boundary layer width
	double y_p;     //non-dimensional wall distance
	double t0;
    double k0;
    double eps0;
    double k_inlet;
    double eps_inlet;
	KE_FIELD* field;
};

enum KEPS_MODEL {
    STANDARD = 0,
    RNG,
    REALIZABLE
};

class KE_RECTANGLE {
public:
        int index;                      // rectangle index
        int comp;
        double area;
        double coords[MAXD];
        int icoords[MAXD];

        KE_RECTANGLE();

        void setCoords(double*,int);
};

class KE_CARTESIAN{

private:

    Front *front;

    static bool activated;
    void applyInitialConditions();

public:
		
    KE_CARTESIAN(Front &front);
	~KE_CARTESIAN();

	KEPS_MODEL keps_model;
	
	int dim;

    double *array;		// for scatter states;
	double *Earray;		// for scatter states;
	double *Karray;		// for scatter states;

	// On topological grid
	RECT_GRID *top_grid;
	double *source;		// for source of parabolic solver;
	double *top_L,*top_U,*top_h,hmin;
	int *top_gmax;
    double lmin, lmax;
	COMPONENT *top_comp;
	KE_PARAMS *eqn_params;
	KE_FIELD *field;
	int comp_size;

	int *lbuf,*ubuf,*gmax;
	int *i_to_I,*I_to_i;		// Index mapping for 1D
	int **ij_to_I,**I_to_ij;	// Index mapping for 2D
	int ***ijk_to_I,**I_to_ijk;	// Index mapping for 3D

	// Sweeping limits
	int imin,jmin,kmin;
	int imax,jmax,kmax;

	enum BC_TYPE { 									// used by ADVECTION
		BC_PERIODIC = 1,
		BC_Extrapolation0 = 2,
		BC_Extrapolation1 = 3,
		BC_InflowOutflow = 4
    };	
	BC_TYPE m_bc[4];								// down, right, up, left 		

	// member data: mesh storage
	std::vector<KE_RECTANGLE> cell_center;

	double m_t;                     // time
	double m_dt;			// time increment
	double min_dt;			// Minimum dt to use non-explicit

	// for parallel partition
	int             NLblocks,ilower,iupper;
    int             *n_dist;

    static void activateKE();
    static void deactivateKE();
    
	// mesh: full cells mesh
	void initMesh(void);		// setup the cartesian grid
	void setDomain(void);
	void setComponent(void);	// init components	
	void setAdvectionDt(void);	// compute time step 	
	//void setBoundary(void);		// set up boundary conditions 	
	void readFrontInteriorState(char*);
	void printFrontInteriorState(char*);

    //TODO: remove these since gamma no longer the same for k and eps
    /*
	    void updateGamma();
	    void updateGamma2d();
	    void updateGamma3d();

        void implicitComputeK(COMPONENT sub_comp);
        void implicitComputeK2d(COMPONENT sub_comp);
        void implicitComputeE(COMPONENT sub_comp);
        void implicitComputeE2d(COMPONENT sub_comp);
    
        void explicitComputeKE(COMPONENT sub_comp);
    
        void explicitComputeKE2d(COMPONENT sub_comp);
        void explicitComputeK2d(COMPONENT sub_comp);
        void explicitComputeE2d(COMPONENT sub_comp);
    
        void explicitComputeKE3d(COMPONENT sub_comp);
        void explicitComputeK3d(COMPONENT sub_comp);
        void explicitComputeE3d(COMPONENT sub_comp);
    */

	//void computeKE();
	void computeAdvection();
	void computeAdvectionK(COMPONENT);
	void computeAdvectionE_STD(COMPONENT);
	void computeAdvectionE_RNG(COMPONENT);
	void computeAdvectionE_REAL(COMPONENT);
	void computeMuTurb();
	double computePointFieldCmu(int*);
	double computePointFieldStrain(int*);
	double computePointFieldC2_RNG(int*);
	double computePointFieldC1_REAL(int*,double);
	void findBdryPoint();

	void computeSource();
	double computeWallPk(int*,int,int,int,
			     HYPER_SURF*,POINTER,double*);
	    /*double computeWallPk(int*,int,int,int,
		    	     HYPER_SURF*,POINTER,double**);*/
	void setSlipBoundary(int*,int,int,int,
			     HYPER_SURF*,POINTER,double**,double*);
	void setTKEatWall(int*,int,int,int,
			     HYPER_SURF*,POINTER,double*,double*);

	// interface functions
	void makeGridIntfc();//doesn't exist
	void deleteGridIntfc();

	// Extra plot functions
	void oneDimPlot(char*);
	void xgraphOneDimPlot(char*);
	void initMovieVariables();
	void augmentMovieVariables(const char*);
	void vtk_plot_temperature2d(char*);
    void vtk_plot3d(const char*);

	// Extra movie functions
	void temperatureMovie(char*);

	void checkStates();

	// parallelization related functions
	//
	void scatMeshArray();
	void setGlobalIndex(COMPONENT);

	// physics calculation
	void setInitialCondition(void);

	void setIndexMap(COMPONENT);
		// for compProjWithSmoothProperty(), 
		// should be changed to use setIndexMap() only.


	// main step function
	void solve(double dt);		

	void getVelocity(double *p, double *U);

	// ----------------------------------------------------------
	// 		utility functions
	// ----------------------------------------------------------

	void getRectangleIndex(int indexRectangle, int &i, int &j);
	void getRectangleIndex(int indexRectangle, int &i, int &j, int &k);
	void getRectangleCenter(int index, double *coords);
	void getRectangleCenter(int index0, int index1, double *coords);
	int getRectangleComponent(int index);	// the center component
	
	double getDistance(double *coords0, double *coords1);
	
			// incompletely implemented
	void getNearestInterfacePoint(double *q,double *p); 
		
	int  getComponent(int *icoords);	
	int  getComponent(double *coords);	
	void save(char *filename);
	void read_params(char*,KE_PARAMS*);
	void computeLiftDrag(Front*);
};


extern double getStateTemp(POINTER);

#endif
