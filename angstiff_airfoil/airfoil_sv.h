#ifndef AIRFOIL_SV
#define AIRFOIL_SV

//TODO: Comments with description of each variable.

struct SPRING_VERTEX
{
    double *x;
    double *v;
	double *f;
	
    double m;           //mass
    double lambda;      //damping coefficient
    
    double *k;
    double *len0;

    int ix;
    int num_nb;

    int *ix_nb;
    double **x_nb;

    int *ix_ajl;            /* index of x_ajl */
    int *ix_ajr;            /* index of x_ajr */
    double **x_ajl;         /* left adjacent side nb vertex */
    double **x_ajr;         /* right adjacent side nb vertex */
	
    double *gam_adj00;      /* gamma between x and x_ajl */
    double *gam_adj01;      /* gamma between x and x_ajr */
    double *gam_adj10;      /* gamma between x_nb and x_ajl */
    double *gam_adj11;      /* gamma between x_nb and x_ajr */

    double *len0_adj00;     /* equilibrium length x-a_ajl */
    double *len0_adj01;     /* equilibrium length x-a_ajr */
    double *len0_adj10;     /* equilibrium length xnb-a_ajl */
    double *len0_adj11;     /* equilibrium length xnb-a_ajr */
    
    double ext_accel[3];
    double *ext_impul; //links to GLOBAL_POINT::impuls
	double *fluid_accel;
    double *bendforce;
	double *other_accel;
};

struct GLOBAL_POINT
{
    long gindex;

    double x[3];
    double v[3];
	double f[3];
	
    double impuls[3]; //links to STATE::impulse
	double fluid_accel[3];
    double bendforce[3];
	double other_accel[3];
};

#endif
