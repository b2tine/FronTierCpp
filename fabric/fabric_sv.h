#ifndef AIRFOIL_SV
#define AIRFOIL_SV

//TODO: Comments with description of each variable.
//      Or even better, rename variables.

        
struct SPRING_VERTEX
{
    double *x;              //position
    double *v;              //velocity
	double *f;              //force
    double *ext_impul;      //velocity comp due to impulse of external force
	int ix;
    int num_nb;             //number of neighboring vertices
    double m;               //mass
    double lambda;          //damping coefficient
    double **x_nb;
    double **v_nb;
	int *ix_nb;
    double *k;              //spring constant
    double *len0;           //equilibrium lengths to neighbors
    double ext_accel[3];    //gravity
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
