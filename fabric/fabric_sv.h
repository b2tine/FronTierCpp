#ifndef FABRIC_SV
#define FABRIC_SV

struct SPRING_VERTEX
{
    int ix;                 // local point index
    double *x;              // position
    double *v;              // velocity
	double *f;              // force
    
    int num_nb;             // number of neighboring vertices
	int *ix_nb;             // neighbor local indices
    double **x_nb;          // neighbor position
    double **v_nb;          // neigbor velocity
    double **x_ajl;         // left adjacent side nb vertex
    double **x_ajr;         // right adjacent side nb vertex
    int *ix_ajl;            // index of x_ajl
    int *ix_ajr;            // index of x_ajr

    double m;               // mass
    double lambda;          // damping coefficient
    double *k;              // spring constant
    double *len0;           // equilibrium lengths to neighbors
    
    double *gam_adj00;      // gamma between x and x_ajl
	double *gam_adj01;		// gamma between x and x_ajr
	double *gam_adj10;		// gamma between x_nb and x_ajl
	double *gam_adj11;		// gamma between x_nb and x_ajr
	double *len0_adj00;		// equilibrium length x - a_ajl
	double *len0_adj01;		// equilibrium length x - a_ajr
	double *len0_adj10;		// equilibrium length xnb - a_ajl
	double *len0_adj11;		// equilibrium length xnb - a_ajr
    
    double ext_accel[3];    // gravity
    double *ext_impul;      // velocity comp due to impulse of external force
	double *fluid_accel;    // acceleration due to fluid pressue difference on fabric
	double *other_accel;    // acceleration for special nodes
};

struct GLOBAL_POINT
{
    double x[3];            // position
    double v[3];            // velocity
	double f[3];            // force
    double impuls[3];
	double fluid_accel[3];  // acceleration due to fluid pressue difference on fabric
	double other_accel[3];  // acceleration for special nodes
    long gindex;            // global index
};

#endif
