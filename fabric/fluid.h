#ifndef FLUID_H
#define FLUID_H

#include <FronTier.h>
#include "state.h"

#include <algorithm>


#define SOLID_COMP 0
#define LIQUID_COMP1 2
#define LIQUID_COMP2 3
#define LIQUID_COMP	3
#define	FILL_COMP 10

#define fluid_comp(comp) (((comp) == LIQUID_COMP1 || \
            comp == LIQUID_COMP2) ? YES : NO)

enum PROJC_METHOD
{
    ERROR_PROJC_SCHEME = -1,
    SIMPLE =  1,
    BELL_COLELLA,
    KIM_MOIN,
    PEROT_BOTELLA
};

enum ADVEC_METHOD
{
	ERROR_ADVEC_SCHEME = -1,
    UPWIND =  1,
    WENO
};

enum ELLIP_METHOD
{
	ERROR_ELLIP_SCHEME = -1,
	SIMPLE_ELLIP = 1,
	DOUBLE_ELLIP,
};

struct NS_SCHEME
{
    PROJC_METHOD projc_method;
    ADVEC_METHOD advec_method;
    ELLIP_METHOD ellip_method;
};

enum EDDY_VISC
{
    BALDWIN_LOMAX = 1,
    MOIN,
    SMAGORINSKY
};

struct F_FIELD
{
	double **vel;			/* Velocities */
	double *temperature;    /* Temperature */
	double *phi;
	double *q;
	double *pres;			/* Pressure */
	double *vort;			/* Vorticity in 2D */
	double *mu;
	double *rho;
	double **grad_q;
	double **f_surf;		// Surface force (such as tension)
	double **old_var;		// For debugging purpose

	double *div_U;
	double *nu_t;			/* Turbulent viscosity */
	double **ext_accel;		/*external forcing from other field*/
};

struct F_PARAMS
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
	double smoothing_radius;
	double ub_speed;
	double min_speed;	/* Limit time step in zero ambient velocity */
	COMPONENT m_comp1;
	COMPONENT m_comp2;
	F_FIELD *field;
	int adv_order;
    boolean total_div_cancellation;
	boolean buoyancy_flow;
	boolean if_buoyancy;
	double  ref_temp;
	boolean if_ref_pres;
	boolean use_eddy_visc;	/* Yes if to use eddy viscosity */
	double  ref_pres;
	EDDY_VISC eddy_visc_model;
	POINTER eddy_params;
	double  Amplitute; 	/*Amplitute of velocity*/
	double	ymax;	   	/* Maximum distance in Baldwin-Lomax model */
	boolean  with_porosity;    /*porosity: 1/0 with/without porosity*/
    double  porous_coeff[2];   /*dp = a*v + b*v^2*/
	double	porosity;
	char base_dir_name[200];
    int base_step;
	boolean scalar_field; /*include scalar field or not*/
	boolean skip_neumann_solver;
};


// fluidinit.cpp
extern void read_Fparams(char*,F_PARAMS*);
extern void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);

// fluidprop.cpp
extern double getStatePres(POINTER);
extern double getStateVort(POINTER);
extern double getStateXvel(POINTER);
extern double getStateYvel(POINTER);
extern double getStateZvel(POINTER);
extern double getStateXimp(POINTER);
extern double getStateYimp(POINTER);
extern double getStateZimp(POINTER);
extern double getStateComp(POINTER);
extern double getStateMu(POINTER);
extern double getStateTemp(POINTER);
extern double getPressure(Front*,double*,double*);
extern void fluid_compute_force_and_torque(Front*,HYPER_SURF*,
        double,double*,double*);
extern void fluid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);



#endif
