#ifndef IFLUID_STATE_H
#define IFLUID_STATE_H


struct STATE
{
    double dens;                /* Density */
    double pres;                /* Pressure */
    double phi;                 /* Potential */
    double vel[MAXD];           /* Velocity */
    double vort;                /* Magnitude of Vorticity in 2D */
	double solute;			    /* For subsurface problem */

	double temperature;         /* For melting with flow problem */
	double mu;			        /* For eddy viscosity TODO: is this actually used?*/
    double tan_stress[MAXD];    /* tangential stress on fluid computed from
                                        turb model + wall functions */
	double fluid_accel[MAXD];   /* acceleration from fluid force */
    double other_accel[MAXD];   /* acceleration for special nodes */

    double impulse[MAXD];       /* Accum impact from external force */


	/* for collision */
	struct UF
    {
        int num_pts;
        POINT* root;
        POINT* tail;
        POINT* next_pt;
    };
    
    UF     impZone;
    double collsn_dt;
    double collsnImpulse[3];
	double collsnImpulse_RG[3];
    double strainImpulse[3];
    double friction[3];
    double avgVel[3];
    double avgVel_old[3];
    double x_old[3];
    int strain_num;
    int collsn_num;
	int collsn_num_RG;
    bool has_collsn;
    bool has_strainlim;
    bool is_fixed;
	bool is_movableRG;
    bool is_stringpt;
};



#endif
