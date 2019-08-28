#ifndef AIRFOIL_SV
#define AIRFOIL_SV

//TODO: Comments with description of each variable.
//      Or even better, rename variables.

        
struct SPRING_VERTEX
{
        double *x;
        double *v;
	double *f;
        double *ext_impul;
	int ix;
        int num_nb;
        double m;
        double lambda;
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
        double impuls[3];
        long gindex;
	double fluid_accel[3];
	double other_accel[3];
};

#endif
