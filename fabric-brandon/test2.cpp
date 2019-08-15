/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/


/*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include "fabric.h"

static  void fabric_driver(Front*);

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;

	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

        SMM_InitFronTier(&front,&f_basic);
        SMM_InitModules(&front);
        SMM_InitSpringMassParams(&front);
        SMM_InitPropagator(&front);
        SMM_InitTestVelFunc(&front);
        SMM_InitTestTimeContrl(&front);
	    
        FT_Draw(&front);

        fabric_driver(&front);
	clean_up(0);
}

static  void fabric_driver(Front *front)
{
        double CFL;
        int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

        CFL = Time_step_factor(front);

	FT_ResetTime(front);

	FT_Save(front);
        FT_Draw(front);

	FT_Propagate(front);
	FT_InteriorPropagate(front);

        FT_SetOutputCounter(front);
	FT_SetTimeStep(front);

        FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);
	
        // For restart debugging 
	if (FT_TimeLimitReached(front) && debugging("restart")) 
	{
	    FT_Save(front);
	    return;
	}
    
        for (;;)
        {
            /* Propagating interface for time step dt */
            if(debugging("CLOCK"))
                reset_clock();

	    start_clock("timeStep");

            FT_Propagate(front);
            FT_InteriorPropagate(front);

            if (debugging("trace"))
            {
                (void) printf("After solve()\n");
                (void) print_storage("at end of time step","trace");
            }

	    FT_AddTimeStepToCounter(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
            if (debugging("step_size"))
            {
                (void) printf("Time step from FrontHypTimeStep(): %f\n",
                                    front->dt);
            }
            
	    /* Output section */

            if (FT_IsSaveTime(front))
	    {
                FT_Save(front);
	    }
        
            if (FT_IsDrawTime(front))
	    {
                FT_Draw(front);
	    }

            if (FT_TimeLimitReached(front))
	    {
                FT_PrintTimeStamp(front);
	    	stop_clock("timeStep");
                break;
	    }

	    /* Time and step control section */

	    FT_TimeControlFilter(front);
	    print_storage("after time loop","trace");

	    FT_PrintTimeStamp(front);
            fflush(stdout);
	    stop_clock("timeStep");
        }
    
        FT_FreeMainIntfc(front);
}       /* end fabric_driver */
