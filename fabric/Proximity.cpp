#include "Proximity.h"


void updateAverageVelocityProximity(POINT** pts)
{
	POINT *p;
	STATE *sl;
	
    for (int i = 0; i < 4; ++i)
    {
        p = pts[i];
    
        if (isStaticRigidBody(p))
            continue;

        sl = (STATE*)left_state(p);

        sl->has_collsn = true;
        sl->avgVel_old = sl->avgVel;
    
        //TODO: does friction impulse also need to be divided
        //      by the number of collisions?
        for (int k = 0; k < 3; ++k)
        {
            sl->avgVel[k] += sl->collsnImpulse[k]/sl->collsn_num;
            sl->avgVel[k] += sl->friction[k]/sl->collsn_num;
            //sl->avgVel[k] += sl->friction[k];
        
            if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
            {
                printf("inf/nan vel[%d]: impulse = %f, friction = %f, collsn_num = %d\n",
                k,sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
                clean_up(ERROR);
            }

            sl->collsnImpulse[k] = 0.0;
            sl->friction[k] = 0.0;
        }

        sl->collsn_num = 0;
    
        // test for RG
        if (sl->collsn_num_RG > 0)
        {
            for (int k = 0; k < 3; ++k)
                sl->avgVel[k] += sl->collsnImpulse_RG[k]/sl->collsn_num_RG;
        
            sl->collsn_num_RG = 0;
        }

    }
	
    /*
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects
    */
}



