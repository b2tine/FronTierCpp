#ifndef AIRFOIL_SV_H
#define AIRFOIL_SV_H

struct _SPRING_VERTEX 
{
    double *x;
    double *v;
	double *f;
    
    int num_nb;
    double m;
    double lambda;
    
    double **x_nb;
    double **x_ajl;			/* left adjacent side nb vertex */
    double **x_ajr;			/* right adjacent side nb vertex */
	
    int *ix_nb;
	int *ix_ajl;			/* index of x_ajl */
	int *ix_ajr;			/* index of x_ajr */
    
    double *k;
    double *len0;
	
    double *gam_adj00;		/* gamma between x and x_ajl */
	double *gam_adj01;		/* gamma between x and x_ajr */
	double *gam_adj10;		/* gamma between x_nb and x_ajl */
	double *gam_adj11;		/* gamma between x_nb and x_ajr */

    double *len0_adj00;		/* equilibrium length x-a_ajl */
	double *len0_adj01;		/* equilibrium length x-a_ajr */
	double *len0_adj10;		/* equilibrium length xnb-a_ajl */
	double *len0_adj11;		/* equilibrium length xnb-a_ajr */

    double ext_accel[3];
};
typedef struct _SPRING_VERTEX SPRING_VERTEX;

#endif
