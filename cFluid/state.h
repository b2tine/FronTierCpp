#ifndef STATE_H
#define STATE_H

#include <FronTier.h>


struct FABRIC_STATE
{   
    /* Fluid-Structure Interaction */
    double fluid_accel[MAXD];       /* acceleration from fluid force */
    double other_accel[MAXD];       /* acceleration for special nodes */
    double impulse[MAXD];            /* Accum impact from external force */
    
    /* Collision Handling */
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


//TODO: make EOS_PARAMS non pointer in structs?
struct EOS_PARAMS
{
    double  gamma;
    double  pinf;
    double  einf;
};


struct STATE : public FABRIC_STATE 
{
    double dens;                    /* Density */
    double pres;                    /* Pressure */
    double engy;                    /* energy density */
    double momn[MAXD];              /* momentum density */
    double vel[MAXD];               /* Velocities */
    double vort;                    /* Vorticity in 2D */
    //double vorticity[MAXDF]; //TODO: 3d vorticity vector
    
    EOS_PARAMS* eos;
    int dim;//TODO: can get rid of this?
    
    double temperature;
    double mu;
};  


// state.cpp
extern double getStateDens(POINTER);
extern double getStatePres(POINTER);
extern double getStateEngy(POINTER);
extern double getStateMu(POINTER);
extern double getStateTemp(POINTER);

extern double getStateXmom(POINTER);
extern double getStateYmom(POINTER);
extern double getStateZmom(POINTER);

extern double getStateXvel(POINTER);
extern double getStateYvel(POINTER);
extern double getStateZvel(POINTER);

extern double getStateVort(POINTER);

extern double getStateXimp(POINTER);
extern double getStateYimp(POINTER);
extern double getStateZimp(POINTER);


#endif
