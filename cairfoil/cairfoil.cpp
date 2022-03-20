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

#include "airfoil.h"

static void airfoil_driver(Front*,CFABRIC_CARTESIAN*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReSetTime;
int RestartStep;
int constrained_propagate;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static EQN_PARAMS eqn_params;
	static AF_PARAMS af_params;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);
	
    if (f_basic.dim != 3)
    {
        printf("\nERROR: dim must be equal to 3\n");
        LOC(); clean_up(EXIT_FAILURE);
    }

	CFABRIC_CARTESIAN g_cartesian(&front);


	/* Initialize basic computational data */

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;
        ReSetTime             	= f_basic.ReSetTime;

        sprintf(restart_state_name,"%s/state.ts%s",restart_name,
			right_flush(RestartStep,7));
        sprintf(restart_name,"%s/intfc-ts%s",restart_name,	
			right_flush(RestartStep,7));
	
    if (pp_numnodes() > 1)
    {
        sprintf(restart_name,"%s-nd%s",restart_name,
                right_flush(pp_mynode(),4));
        sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                right_flush(pp_mynode(),4));
	}

	af_params.num_np = 1;
    FT_VectorMemoryAlloc((POINTER*)&af_params.node_id,1,sizeof(int));
    af_params.node_id[0] = 0;

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	if (debugging("trace"))
        (void) printf("Passed FT_StartUp()\n");

    eqn_params.dim = f_basic.dim;
    front.extra1 = (POINTER)&eqn_params;
    front.extra2 = (POINTER)&af_params;
    
    read_cFluid_params(in_name,&eqn_params);
    
    if (debugging("trace")) 
        (void) printf("Passed read_cFluid_params()\n");

    level_func_pack.pos_component = GAS_COMP2;
	if (!RestartRun)
	{
	    FT_InitIntfc(&front,&level_func_pack);
        
        initFabricModules(&front);
	    
        if (f_basic.dim == 3)
            initIsolated3dCurves(&front);

        if (consistent_interface(front.interf) == NO)
        {
            printf("consistent_interface(front.interf) == NO\n");
            LOC(); clean_up(ERROR);
        }
	    
        if (f_basic.dim == 3 && debugging("init_intfc"))
	    {
            char gvdir[100];
            sprintf(gvdir,"%s/gv-init",out_name);
            gview_plot_interface(gvdir,front.interf);
	    }

	    read_dirichlet_bdry_data(in_name,&front);

	    if (f_basic.dim < 3)
            FT_ClipIntfcToSubdomain(&front);
	}
	else
	{
	    read_dirichlet_bdry_data(in_name,&front);
	}

	/* Time control */
	FT_ReadTimeControl(in_name,&front);
        
	if (!RestartRun)
	{
	    optimizeElasticMesh(&front);
	    set_equilibrium_mesh(&front);
	    FT_SetGlobalIndex(&front);
        static_mesh(front.interf) = YES;
	}

	/* Initialize velocity field function */
	setMotionParams(&front);

    front._scatter_front_extra = scatterAirfoilExtra;
	front._compute_force_and_torque = cfluid_compute_force_and_torque;
	
    g_cartesian.findStateAtCrossing = af_find_state_at_crossing;
	g_cartesian.initMesh();
    
    g_cartesian.writeMeshFileVTK();
    g_cartesian.writeCompGridMeshFileVTK();
        
	    //g_cartesian.getInitialState = zero_state;

    if (debugging("sample_velocity"))
        g_cartesian.initSampleVelocity(in_name);

    if (RestartRun)
	{
	    if (ReSetTime) 
	    {
		    /* forbidden if restart with inherited states */
	    	readAfExtraData(&front,restart_state_name);
	    	modifyInitialization(&front);
	    	read_dirichlet_bdry_data(in_name,&front);
		    g_cartesian.initMesh();
            g_cartesian.setInitialStates();
	    }
	    else
	    {
            if (!af_params.no_fluid)
            {
                FT_MakeGridIntfc(&front);
                coating_mono_hyper_surf(&front);
                g_cartesian.applicationSetComponent();
                FT_FreeGridIntfc(front);
            }

            readFrontStates(&front,restart_state_name);
            g_cartesian.readInteriorStates(restart_state_name);
	    	readAfExtraData(&front,restart_state_name);
	    }
	}
    else
    {
        g_cartesian.setInitialStates();
    }

    g_cartesian.initMovieVariables();
    initMovieStress(in_name,&front);

	if (!RestartRun || ReSetTime)
	    resetFrontVelocity(&front);

    if (debugging("trace"))
        (void) printf("Passed state initialization()\n");

	/* Propagate the front */
	airfoil_driver(&front,&g_cartesian);

	clean_up(0);
}

void airfoil_driver(Front *front,
        CFABRIC_CARTESIAN* g_cartesian)
{
    double CFL;
    int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

    CFL = Time_step_factor(front);
	Tracking_algorithm(front) = STRUCTURE_TRACKING;
    
        //TwoStepIntfc(front) = YES;

	(void) printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun || ReSetTime)
	{
	    FT_ResetTime(front);
	    setStressColor(front);
	    FT_Save(front);

        g_cartesian->printFrontInteriorStates(out_name);
        printAfExtraData(front,out_name);

        if (!RestartRun && !ReSetTime)
        {
            FT_Draw(front);
        }

	    FrontPreAdvance(front);
	    FT_Propagate(front);
        FT_RelinkGlobalIndex(front);

	    if (!af_params->no_fluid)
	    {
            if (debugging("trace")) printf("Calling ifluid solve()\n");
            g_cartesian->solve(front->dt);
            if (debugging("trace")) printf("Passed ifluid solve()\n");
	    }
	    print_airfoil_stat(front,out_name);

        FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);

	    if (!af_params->no_fluid)
	    {
	    	front->dt = std::min(front->dt,CFL*g_cartesian->max_dt);
	    }

        front->dt = std::min(front->dt,springCharTimeStep(front));
	}
	else
	{
	    FT_SetOutputCounter(front);
	}
	
    FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);

    for (;;)
    {
        /* Propagating interface for time step dt */
	    if (debugging("CLOCK"))
            reset_clock();

	    start_clock("time_step");
	    if (debugging("trace"))
            (void) printf("Before FT_Propagate()\n");

	    if (!af_params->no_fluid)
	    {
	    	coating_mono_hyper_surf(front);
	    	g_cartesian->applicationSetComponent();
	    }

	    break_strings(front);
	    
        FrontPreAdvance(front);
        FT_Propagate(front);
        FT_RelinkGlobalIndex(front);
        FT_InteriorPropagate(front);
	    
        if (debugging("trace"))
            printf("Passed FT_Propagate()\n");

        if (!af_params->no_fluid)
	    {
	    	coating_mono_hyper_surf(front);
            g_cartesian->applicationSetStatesNEW();
	    }
        
	    if (!af_params->no_fluid)
	    {
            g_cartesian->solve(front->dt);
            FT_FreeGridIntfc(front);
	        FT_MakeGridIntfc(front);
	    }
	    else
        {
            g_cartesian->max_dt = HUGE;
        }

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
            printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
        
        front->dt = std::min(front->dt,CFL*g_cartesian->max_dt);

        if (debugging("step_size"))
            printf("Time step from g_cartesian->max_dt(): %f\n",front->dt);

	    /* Output section */

	    print_airfoil_stat(front,out_name);

        /*
        //TODO: write printEnstrophy()
        //
        if (!af_params->no_fluid)
             g_cartesian->printEnstrophy();
        */

        setStressColor(front);

        if (FT_IsSaveTime(front))
	    {
            FT_Save(front);
            g_cartesian->printFrontInteriorStates(out_name);
	    	printAfExtraData(front,out_name);
	    }
        
        if (debugging("trace"))
            (void) printf("After print output()\n");
        
        if (FT_IsDrawTime(front))
	    {
            FT_Draw(front);
            g_cartesian->writeMeshComponentsVTK();
	    }

        if (FT_TimeLimitReached(front))
	    {
            FT_PrintTimeStamp(front);
            stop_clock("time_step");
            break;
	    }

	    /* Time and step control section */

	    FT_TimeControlFilter(front);
	    print_storage("after time loop","trace");

	    FT_PrintTimeStamp(front);
        fflush(stdout);
	    stop_clock("time_step");
    }

    FT_FreeMainIntfc(front);
}       /* end airfoil_driver */

