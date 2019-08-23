#ifndef IFLUID_STATE_H
#define IFLUID_STATE_H

//TODO: Split into individual classes for iFluid and Collision.
//      Then can make an FSI class for airfoil via composition
//      of the iFluid and Collision STATE structures.


struct STATE {
        double dens;                    /* Density */
        double pres;                    /* Pressure */
        double phi;                     /* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
	double solute;			/* For subsurface problem */

	double temperature;             /* For melting with flow problem */
	double vapor;                   /* For climate problem */
        double supersat;		/* For climate problem */
	double mu;			/* For eddy viscosity */
	double fluid_accel[MAXD];       /* acceleration from fluid force */
        double other_accel[MAXD];       /* acceleration for special nodes */

        double impulse[MAXD];            /* Accum impact from external force */


	/* for collision */
	struct UF   {
        POINT* next_pt;
        POINT* root;
        POINT* tail;
        int num_pts;
        };
    
        UF     impZone;
        double collsnImpulse[3];
	double collsnImpulse_RG[3];
        double friction[3];
        double avgVel[3];
        double avgVel_old[3];
        double x_old[3];
        int    collsn_num;
	int    collsn_num_RG;
        bool   has_collsn;
        bool   is_fixed;
	bool   is_movableRG;
};



#endif
