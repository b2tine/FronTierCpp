#ifndef AIRFOIL_SV
#define AIRFOIL_SV

//TODO: Comments with description of each variable.

struct SPRING_VERTEX
{
    double *x;
    double *v;
	double *f;
    double *ext_impul; //links to GLOBAL_POINT::impuls
	int ix;
    int num_nb;
    double m;           //mass
    double lambda;      //friction/damping coefficient
    double **x_nb;
    double **v_nb;
	int *ix_nb;
    double *k;
    double *len0;
    double ext_accel[3];
	double *fluid_accel;
	double *other_accel;
};

struct GLOBAL_POINT
{
    double x[3];
    double v[3];
	double f[3];
    double impuls[3]; //links to STATE::impulse
    long gindex;
	double fluid_accel[3];
	double other_accel[3];
};

#endif
