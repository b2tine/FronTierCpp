#ifndef STATE_H
#define STATE_H

struct STATE
{
    double dens;                    /* Density */
    double pres;                    /* Pressure */
    double vel[MAXD];               /* Velocities */
    double vort;                    /* Vorticity in 2D */
	double temperature;             /* For melting with flow problem */
	double mu;			            /* For eddy viscosity */
	double fluid_accel[MAXD];       /* acceleration from fluid force */
    double other_accel[MAXD];       /* acceleration for special nodes */
    double impulse[MAXD];           /* Accum impact from external force */
    
    double bendforce[MAXD];         /* bending force */


    //TODO: Move collision state data into a proxy class for POINTs.
	
    /* for collision */
	struct UF   
    {
        int num_pts;
        POINT* root;
        POINT* tail;
        POINT* next_pt;
    };
    
    UF impZone;
    double collsn_dt;
    double collsnImpulse[3];
    double collsnImpulse_RG[3];
    double strainImpulse[3];
    double friction[3];
    double avgVel[3];
    double avgVel_old[3];
    double avgVel_postprox[3];
    double x_old[3];
    int strain_num;
    int collsn_num;
    int collsn_num_RG;
    bool has_proximity;
    bool has_strainlim_prox;
    bool has_collsn;
    bool has_strainlim_collsn;
    bool is_fixed;
    bool is_movableRG;
    bool is_stringpt;
};

#endif
