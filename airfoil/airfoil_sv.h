#ifndef AIRFOIL_SV
#define AIRFOIL_SV

//TODO: Comments with description of each variable.
//      Or even better, rename variables.

//TODO: To implement artificial viscosity do we need pointers to
//      the neighboring vertices so the velocity components due to
//      the impulse of external force can be subtracted off.
//      Then we would be using the relative velocity of only the
//      internal internal component resulting from the impulse of
//      the spring force between vertices. See afsetd.cpp line 289,
//      and compare with line 307.

struct SPRING_VERTEX
{
    double *x;
    double *v;
	double *f;
    double *ext_impul;
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
    double impuls[3];
    long gindex;
	double fluid_accel[3];
	double other_accel[3];
};

#endif
