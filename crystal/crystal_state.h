/**********************************************************************
 * 		crystal_state.h					      *
 **********************************************************************/


#ifndef CRYSTAL_STATE_H
#define CRYSTAL_STATE_H


struct Crystal_STATE {
        double dens;                    /* Density */
        double pres;                    /* Pressure */
        double phi;                     /* Potential */
        double vel[MAXD];               /* Velocities */
        double vort;                    /* Vorticity in 2D */
        double solute;                  /* For subsurface problem */
        double mu;                      /* For eddy viscosity */
};

#endif

