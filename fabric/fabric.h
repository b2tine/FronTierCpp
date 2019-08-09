#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <FronTier.h>

#include "fluid.h"
#include "state.h"
#include "rigidbody.h"

#include "airfoil_sv.h"

#include <cassert>
#include <vector>
#include <map>


struct CONSTR_PARAMS
{
    double N[MAXD];         /* normal of the plane */
    double P[MAXD];         /* a point on the plane */
    double cen[MAXD];       /* center of the vent */
    double radius;          /* radius of the vent */
};

enum PERTURBATION_TYPE
{
	NO_PERT	= 1,
	PARALLEL_RAND_PERT,
	ORTHOGONAL_RAND_PERT,
	LINEAR_PERT,
	RADIAL_PERT,
	SINE_PERT
};

enum STRING_NODE_TYPE
{
	FIXED_END =	1,
	FREE_END,
	LOADED_END
};

enum AF_NODE_TYPE
{
	UNKNOWN_AF_NODE = -1,
	LOAD_NODE = 1,
	RG_STRING_NODE,
	GORE_NODE,
	STRING_NODE,
	PRESET_NODE,
	SEC_LOAD_NODE,
    THR_LOAD_NODE
};

enum SPRING_MODEL
{
	UNKNOWN_MODEL 	= -1,
	MODEL1	= 1,
	MODEL2,
	MODEL3
};

struct PERT_PARAMS
{
	PERTURBATION_TYPE pert_type;
	int dir;
	double x0,xl,xu;
	double pert_amp;
	double cen[MAXD];
	double pert_radius;
};

struct AF_PARAMS
{
        int dim;
        SPRING_MODEL spring_model;
	boolean no_fluid;
	boolean is_parachute_system;
	boolean attach_gores;
	boolean attach_fixer;
	boolean cut_vent;
        boolean use_total_mass;
	boolean use_gpu;
	PERT_PARAMS pert_params;
	STRING_NODE_TYPE start_type;
	STRING_NODE_TYPE end_type;
	
        double gore_len_fac;
        double gravity[MAXD];		/* gravitational force */
	double payload;
	double ks;			/* spring constant of surface */
	double kl;			/* spring constant of string curves */
	double kg;                      /* spring constant of gore curves */
	double lambda_s;		/* damping/friction factor of surface */
	double lambda_l;		/* damping/friction factor of string curves */
	double lambda_g;                /* damping factor of gore curves */
	double m_s;			/* point mass of surface */
	double m_l;			/* point mass of string curves */
	double m_g;                     /* point mass of gore curves */
	double total_string_mass;	/* Total mass of string chord */
	double total_canopy_mass;	/* Total mass of string chord */
        double total_gore_mass;         /* Total mass of gore */
	boolean with_porosity;          /* with or without porosity*/
	double porous_coeff[2];         /* viscous and inertial coefficients*/
	double gamma;			/* canopy porosity */
	double area_dens;		/* canopy area density */
	int    n_sub;			/* number of sub-steps for tan prop */
	int    num_opt_round;		/* number of mesh optimizations rounds*/
	int    num_smooth_layers;	/* number of layer to smooth high frequency velocity */
	int    num_np;			/* number of master node to run spring model */
	int    node_id[10];		/* master node id */
	
        double break_strings_time;	/* time to break some strings */
	int    break_strings_num;	/* number of strings to break */
	int    *break_strings_gindex;	/* gindex of strings to break */
	double unequal_coeff;		/* the unequal coefficient */
	int    unequal_strings_num;	/* number of unequal strings */
	int    *unequal_strings_gindex; /* gindex of unequal strings */

	std::vector<CURVE*> string_curves;	/* string curves in order */
	std::map<int,int> string_hash;	/* map from string gindex to string 
					   id, for users' convenience */
}; 


/*	airfoil.cpp functions */

struct DOUBLE_VORTEX_PARAMS
{
    double i1,i2;
    double cen1[2],cen2[2];
};

struct TOROIDAL_PARAMS
{
    double tcen[MAXD]; /* toroidal center */
	double R0;	/* distance between poloidal and toroidal centers */
	double v0;	/* amplitude of velocity */
	double stop_time;
};

struct PARABOLIC_PARAMS
{
    double cen[MAXD]; /* parabolic center */
	double v0;	  /* velocity at center */
	double a;	  /* concavity (downward) */
	double stop_time;
};

struct SINGULAR_PARAMS
{
    double cen[MAXD]; /* parabolic center */
	double R;	  /* radius of for v0 */
	double v0;	  /* velocity at center */
	double stop_time;
};

struct STRING_PARAMS
{
	int num_strings;
	double start_angle;
	double coords_load[MAXD];
	double cen[MAXD];
	double shift[MAXD];
	double theta;
	double phi;
	double L[MAXD],U[MAXD];
	double P[MAXD];
};

struct PARALLEL_GORE_PARAMS
{
    int gores_n;
    double gores_start_x;
    double gores_dis;
    double coords_load[MAXD];
};

enum LOAD_TYPE
{
	NO_LOAD = 0,
	FREE_LOAD,
	RIGID_LOAD
};

struct BDRY_PARAMS
{
    boolean lower_bdry[MAXD];
    boolean upper_bdry[MAXD];
    double L[MAXD];         /* Lower bounds of box */
    double U[MAXD];         /* Upper bounds of box */
	LOAD_TYPE lower_side[MAXD];
	LOAD_TYPE upper_side[MAXD];
	double lower_mass[MAXD];
	double upper_mass[MAXD];
	double lower_force[MAXD][MAXD];
	double upper_force[MAXD][MAXD];
}; 

struct C_PARAMS
{
	LOAD_TYPE load_type;
	int dir;
	double load_mass;
	double point_mass;
	double force[MAXD];
	double ave_accel;
};

struct AF_NODE_EXTRA
{
	AF_NODE_TYPE af_node_type;
	double accum_impulse;
};

struct REGISTERED_PTS
{
	int num_pts;
	int *global_ids;
};

struct ELASTIC_SET
{
	Front *front;
    NODE *load_node;
	SURFACE *surfs[100];
	CURVE *curves[1000];
	NODE *nodes[1000];
	int num_surfs;
	int num_curves;
	int num_nodes;
	double ks;
	double kl;
	double kg;
	double lambda_s;
	double lambda_l;
	double lambda_g;
	double m_s;
	double m_l;
	double m_g;
	int num_verts;		/* Total number of spring-mass points */
	double dt_tol;
	double dt;
};

//NOTE: these never had a definition
/*
void read_movie_options(char*,F_PARAMS*);
void liquid_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
*/

// afinit.cpp
extern void read_Fparams(char*,F_PARAMS*);
extern void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
extern void restart_set_dirichlet_bdry_function(Front*);
extern void setInitialIntfcAF(Front*,LEVEL_FUNC_PACK*,char*);
extern void initPerturbation3d(Front*);

/* afinit3d.cpp */
extern void initEllipticSurf(FILE*,Front*,LEVEL_FUNC_PACK*);
extern void initParabolicSurf(FILE*,Front*,LEVEL_FUNC_PACK*);
extern void initPlaneSurf(FILE*,Front*,LEVEL_FUNC_PACK*);
extern void initAirbag(FILE*,Front*,LEVEL_FUNC_PACK*);
extern void initIsolated3dCurves(Front*);

// afprop.cpp
extern void airfoil_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void elastic_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
extern void fourth_order_elastic_set_propagate(Front*,double);
extern void airfoil_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
extern int airfoil_node_propagate(Front*,POINTER,NODE*,NODE*,RPROBLEM**,
                               double,double*,NODE_FLAG,POINTER);
extern void coating_mono_hyper_surf(Front*);
extern int airfoil_node_propagate(Front*,POINTER,NODE*,NODE*,RPROBLEM**,
			double,double*,NODE_FLAG,POINTER);

// aftest.cpp
extern void second_order_elastic_curve_propagate(Front*,Front*,INTERFACE*,
                                CURVE*,CURVE*,double);
extern void second_order_elastic_surf_propagate(Front*,double);
extern void set_equilibrium_mesh(Front*);
extern void print_airfoil_stat(Front*,char*);
extern void fixed_length_tan_curve_propagate(Front*,Front*,INTERFACE*,
                                CURVE*,CURVE*,double);
extern void fourth_order_elastic_curve_propagate(Front*,Front*,INTERFACE*,
                                CURVE*,CURVE*,double);
extern void fourth_order_elastic_surf_propagate(Front*,double);
extern void resolve_wall_collision(Front*,SPRING_VERTEX*,int);

// afcnpy.cpp
extern void compute_total_canopy_force(Front*,double*,double*);
extern int airfoil_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
                                double*);
extern int af_find_state_at_crossing(Front*,int*,GRID_DIRECTION,
                        int,POINTER*,HYPER_SURF**,double*);
extern void assign_node_field(NODE*,double**,double**,int*);
extern void assign_curve_field(CURVE*,double**,double**,int*);
extern void assign_surf_field(SURFACE*,double**,double**,int*);
extern void compute_surf_accel1(ELASTIC_SET*,SURFACE*,double**,double**,
				double**,int*);
extern void compute_surf_accel2(ELASTIC_SET*,SURFACE*,double**,double**,
				double**,int*);
extern void compute_curve_accel1(ELASTIC_SET*,CURVE*,double**,double**,
				double**,int*);
extern void compute_node_accel1(ELASTIC_SET*,NODE*,double**,double**,double**,
				int*);
extern void compute_curve_accel2(ELASTIC_SET*,CURVE*,double**,double**,
				double**,int*);
extern void compute_node_accel2(ELASTIC_SET*,NODE*,double**,double**,double**,
				int*);
extern void compute_curve_accel3(ELASTIC_SET*,CURVE*,double**,double**,
				double**,int*);
extern void compute_node_accel3(ELASTIC_SET*,NODE*,double**,double**,double**,
				int*);
extern void compute_center_of_mass_velo(ELASTIC_SET*);
extern void propagate_surface(ELASTIC_SET*,SURFACE*,double**,int*);
extern void propagate_curve(ELASTIC_SET*,CURVE*,double**,int*);
extern void propagate_node(ELASTIC_SET*,NODE*,double**,int*);
extern boolean is_registered_point(SURFACE*,POINT*);
extern void scatterAirfoilExtra(Front*);
extern void setSpecialNodeForce(Front*,double);
extern void break_strings(Front*);
extern void record_break_strings_gindex(Front*);
extern void set_unequal_strings(Front*);
extern void collectNodeExtra(Front*,INTERFACE*,int);

// afsetd.cpp
extern void assembleParachuteSet(INTERFACE*,ELASTIC_SET*);
extern void count_vertex_neighbors(ELASTIC_SET*,SPRING_VERTEX*);
extern void set_spring_vertex_memory(SPRING_VERTEX*,int);
extern void compute_spring_accel1(SPRING_VERTEX*,double*,int);
extern void generic_spring_solver(SPRING_VERTEX*,int,int,int,double);
extern void set_vertex_impulse(ELASTIC_SET*,SPRING_VERTEX*);
extern void set_geomset_velocity(ELASTIC_SET*,SPRING_VERTEX*);
extern void link_point_set(ELASTIC_SET*,GLOBAL_POINT**,GLOBAL_POINT*);
extern void set_vertex_neighbors(ELASTIC_SET*,SPRING_VERTEX*,GLOBAL_POINT**);
extern void set_node_spring_vertex(ELASTIC_SET*,NODE*,SPRING_VERTEX*,
				int*,GLOBAL_POINT**);
extern void set_curve_spring_vertex(ELASTIC_SET*,CURVE*,SPRING_VERTEX*,
				int*,GLOBAL_POINT**);
extern void set_surf_spring_vertex(ELASTIC_SET*,SURFACE*,SPRING_VERTEX*,
				int*,GLOBAL_POINT**);
extern void get_point_set_from(ELASTIC_SET*,GLOBAL_POINT**);
extern void put_point_set_to(ELASTIC_SET*,GLOBAL_POINT**);
extern void set_elastic_params(ELASTIC_SET*,double);
extern void merge_global_point_set(GLOBAL_POINT**,GLOBAL_POINT*,int);
extern void copy_from_client_point_set(GLOBAL_POINT**,GLOBAL_POINT*,int,
				double*,double*);
extern void copy_to_client_point_set(GLOBAL_POINT**,GLOBAL_POINT*,int);
extern void set_vertex_impulse(ELASTIC_SET*,GLOBAL_POINT**);
extern void set_geomset_velocity(ELASTIC_SET*,GLOBAL_POINT**);

// afvelo.cpp
extern void setMotionParams(Front*);
extern void resetFrontVelocity(Front*);
extern void convert_to_point_mass(Front*, AF_PARAMS*);
extern void initVelocityFunc(FILE*,Front*);
extern int countSurfPoints(INTERFACE*);
extern int countStringPoints(INTERFACE*,boolean);

// afmodule.cpp -- in parachute/modules.cpp
extern void initParachuteDefault(Front*);
extern void initParachuteModules(Front*);

// afdata.cpp
extern void printAfExtraData(Front*,char*);
extern void readAfExtraData(Front*,char*);
extern void printHyperSurfQuality(Front*);
extern void print_elastic_params(ELASTIC_SET);
extern void optimizeElasticMesh(Front*);
extern void modifyInitialization(Front*);
extern void setStressColor(Front*);
extern void vtkPlotSurfaceStress(Front*);
extern void poisson_ratio(Front*);
extern void initMovieStress(char*,Front*);
extern int numOfGoreHsbdry(INTERFACE*);
extern int numOfMonoHsbdry(INTERFACE*);
extern int numOfGoreNodes(INTERFACE*);
extern int arrayOfMonoHsbdry(INTERFACE*,CURVE**);
extern int arrayOfGoreHsbdry(INTERFACE*,CURVE**);
extern int getGoreNodes(INTERFACE*,NODE**);
extern boolean is_load_node(NODE*);
extern boolean is_rg_string_node(NODE*);
extern boolean is_gore_node(NODE*);
extern boolean is_bdry_node(NODE*);
extern double springCharTimeStep(Front*);	// spring characteristic time

// cgal.cpp
extern void CgalCanopySurface(FILE*,Front*,SURFACE**);
extern void InstallNewLoadNode(Front*,int);

// fabric_api_internal.cpp
extern void SMM_InitFronTier(Front*,F_BASIC_DATA*);
extern void SMM_InitModules(Front*);
extern void SMM_InitPerturbation(Front*);
extern void SMM_InitPropagator(Front*);
extern void SMM_InitSpringMassParams(Front*);
extern void SMM_InitTestVelFunc(Front*);
extern void SMM_InitTestTimeContrl(Front*);

#endif
