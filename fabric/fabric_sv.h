#ifndef FABRIC_SV
#define FABRIC_SV

struct SPRING_VERTEX
{
    double *x;              //position
    double *v;              //velocity
	double *f;              //force
	double *bendforce;      //bending force
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
	double *fluid_accel;    //acceleration due to fluid pressue difference on fabric
	double *other_accel;    //acceleration for special nodes
};

struct GLOBAL_POINT
{
    double x[3];            //position
    double v[3];            //velocity
	double f[3];            //force
	double bendforce[3];    //bending force
    double impuls[3];
	double fluid_accel[3];  //acceleration due to fluid pressue difference on fabric
	double other_accel[3];  //acceleration for special nodes
    long gindex;            //global index
};

#endif
