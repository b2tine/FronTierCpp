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

extern Front *SMM_GetFront()
{
        static Front front;
        return &front;
}       /* end SMM_GetFront */

extern F_BASIC_DATA *SMM_GetBasicData()
{
        static F_BASIC_DATA f_basic;
        return &f_basic;
}       /* end SMM_GetBasicData */

extern void SMM_InitCpp(int argc, char **argv)
{
        Front *front = SMM_GetFront();
        F_BASIC_DATA *f_basic = SMM_GetBasicData();
        static LEVEL_FUNC_PACK level_func_pack;
        static AF_PARAMS af_params;

        f_basic->size_of_intfc_state = sizeof(STATE);
        
        FT_Init(argc,argv,f_basic);
        FT_ReadSpaceDomain(f_basic->in_name,f_basic);
        
        af_params.num_np = 1;
        af_params.node_id[0] = 0;
        front->extra2 = (POINTER)&af_params;

        if (!f_basic->RestartRun)
        {
            FT_StartUp(front,f_basic);
            FT_InitDebug(f_basic->in_name);
        
            if (FT_Dimension() == 2) // initialization using old method
                setInitialIntfcAF(front,&level_func_pack,InName(front));
            else
                level_func_pack.pos_component = LIQUID_COMP2;
        
            FT_InitIntfc(front,&level_func_pack);
            FT_ResetTime(front);
        }
        else
        {
            SMM_Restart(front,f_basic);
        }
}       /* end SMM_InitCpp */

extern void SMM_Restart(Front *front, F_BASIC_DATA *f_basic)
{
        char *restart_name            = f_basic->restart_name;
        char *restart_state_name      = f_basic->restart_state_name;

        sprintf(restart_state_name,"%s/state.ts%s",restart_name,
			right_flush(f_basic->RestartStep,7));
        sprintf(restart_name,"%s/intfc-ts%s",restart_name,	
			right_flush(f_basic->RestartStep,7));
	
        if (pp_numnodes() > 1)
        {
            sprintf(restart_name,"%s-nd%s",restart_name,
                    right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                    right_flush(pp_mynode(),4));
        }

        FT_StartUp(front,f_basic);
        FT_InitDebug(f_basic->in_name);
        readAfExtraData(front,restart_state_name);
        FT_SetOutputCounter(front);
}

extern void SMM_StartUpStep()
{
        Front *front = SMM_GetFront();
        F_BASIC_DATA *f_basic = SMM_GetBasicData();
        FILE *infile = fopen(InName(front),"r");
        
        int  dim = front->rect_grid->dim;
	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
        double CFL = Time_step_factor(front);

	    SMM_Save();
        SMM_Plot();

        if (!f_basic->RestartRun)
        {
            FT_Propagate(front);
            FT_RelinkGlobalIndex(front);
	        FT_InteriorPropagate(front);

            FT_SetOutputCounter(front);
            FT_SetTimeStep(front);

            FT_TimeControlFilter(front);
        }
        else
        {
            setSpecialNodeForce(front,af_params->kl);
            FT_SetOutputCounter(front);
	        FT_TimeControlFilter(front);
        }

        FT_PrintTimeStamp(front);
}

extern void SMM_TimeMarch()
{
        Front *front = SMM_GetFront();
        F_BASIC_DATA *f_basic = SMM_GetBasicData();
        FILE *infile = fopen(InName(front),"r");
        
        int  dim = front->rect_grid->dim;
	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
        double CFL = Time_step_factor(front);

        for (;;)
        {
            /* Propagating interface for time step dt */
            if(debugging("CLOCK"))
                reset_clock();

	        start_clock("timeStep");

            FrontPreAdvance(front);
            FT_Propagate(front);
            FT_RelinkGlobalIndex(front);
            FT_InteriorPropagate(front);

            if (debugging("trace"))
            {
                printf("After solve()\n");
                print_storage("at end of time step","trace");
            }

            FT_AddTimeStepToCounter(front);

            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

            FT_SetTimeStep(front);
            if (debugging("step_size"))
            {
                printf("Time step from FrontHypTimeStep(): %f\n",
                        front->dt);
            }
                
            /* Output section */

            if (FT_IsSaveTime(front))
                SMM_Save();
        
            if (FT_IsDrawTime(front))
                SMM_Plot();

            if (FT_TimeLimitReached(front))
	        {
                FT_PrintTimeStamp(front);
                if (!FT_IsSaveTime(front))
                    SMM_Save();
                if (!FT_IsDrawTime(front))
                    SMM_Plot();
	    	    stop_clock("timeStep");
                break;
	        }

            /* Time and step control section */

            FT_TimeControlFilter(front);
            FT_PrintTimeStamp(front);

            fflush(stdout);
            stop_clock("timeStep");
        }
    
        FT_FreeMainIntfc(front);
}

#ifdef __cplusplus
extern "C" {
#endif

extern void SMM_Init(char* args)
{
        std::string s(args);
        std::vector<std::string> argstrings;
        argstrings.insert(argstrings.begin(),"dummyarg");
        int argc = 1;

        std::string arg;
        std::stringstream ss(s);
        while (ss >> arg)
        {
            argstrings.push_back(arg);
            argc++;
        }

        char* argv[argc];
        for (int i = 0; i < argc; ++i)
        {
            argv[i] = const_cast<char*>(argstrings[i].c_str());
        }

        SMM_InitCpp(argc,argv);
}

extern void SMM_InitModules()
{
        Front *front = SMM_GetFront();
        F_BASIC_DATA *fbasic = SMM_GetBasicData();

        if (!fbasic->RestartRun)
        {
            FILE *infile = fopen(InName(front),"r");
            printf("\n");
            CursorAfterString(infile,"Start parameters for modules");
            printf("\n");
            fclose(infile);

            if (FT_Dimension() == 3) // 2D initialization used old method
            {
                initParachuteModules(front);
                initPerturbation3d(front);
                initIsolated3dCurves(front);
                optimizeElasticMesh(front);
            }
            set_equilibrium_mesh(front);
            FT_SetGlobalIndex(front);
            static_mesh(front->interf) = YES;
        }
}       /* end SMM_InitModules */

extern void SMM_InitPropagator()
{
        Front *front = SMM_GetFront();
        int dim = front->rect_grid->dim;

        Tracking_algorithm(front) = STRUCTURE_TRACKING;
        front->_point_propagate = airfoil_point_propagate;
        front->curve_propagate = airfoil_curve_propagate;
        front->node_propagate = airfoil_node_propagate;
        front->interior_propagate = fourth_order_elastic_set_propagate;
        front->_compute_force_and_torque = fluid_compute_force_and_torque;
}       /* end SMM_InitPropagator */

extern void SMM_InitSpringMassParams()
{
        Front *front = SMM_GetFront();
	FILE *infile = fopen(InName(front),"r");
	int i,dim = front->rect_grid->dim;
	char string[100];
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

    printf("\n");
    CursorAfterString(infile,"Start parameters for spring-mass system");
    printf("\n");

    af_params->n_sub = 1;
	if (CursorAfterStringOpt(infile,"Enter interior sub step number:"))
        {
	    fscanf(infile,"%d",&af_params->n_sub);
	    (void) printf("%d\n",af_params->n_sub);
        }

        af_params->use_total_mass = NO;
        if (CursorAfterStringOpt(infile,"Enter yes to use total mass:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                af_params->use_total_mass = YES;
        }

	if (FT_FrontContainWaveType(front,ELASTIC_BOUNDARY))
	{
	    CursorAfterString(infile,"Enter fabric spring constant:");
            fscanf(infile,"%lf",&af_params->ks);
            (void) printf("%f\n",af_params->ks);
            
            CursorAfterString(infile,"Enter fabric damping constant:");
            fscanf(infile,"%lf",&af_params->lambda_s);
            (void) printf("%f\n",af_params->lambda_s);

            CursorAfterString(infile,"Enter fabric friction constant:");
            fscanf(infile,"%lf",&af_params->mu_s);
            (void) printf("%f\n",af_params->mu_s);

            if (af_params->use_total_mass)
            {
                CursorAfterString(infile,"Enter fabric total mass:");
                fscanf(infile,"%lf",&af_params->total_canopy_mass);
                (void) printf("%f\n",af_params->total_canopy_mass);
            }
            else
            {
	            CursorAfterString(infile,"Enter fabric point mass:");
                fscanf(infile,"%lf",&af_params->m_s);
                (void) printf("%f\n",af_params->m_s);
	        }
	}

	if ((dim == 2 && FT_FrontContainWaveType(front,ELASTIC_STRING)) || 
	    (dim == 3 && FT_FrontContainHsbdryType(front,STRING_HSBDRY)) )
	{
        CursorAfterString(infile,"Enter string spring constant:");
        fscanf(infile,"%lf",&af_params->kl);
        (void) printf("%f\n",af_params->kl);
        CursorAfterString(infile,"Enter string damping constant:");
        fscanf(infile,"%lf",&af_params->lambda_l);
        (void) printf("%f\n",af_params->lambda_l);
        if (af_params->use_total_mass)
        {
            CursorAfterString(infile,"Enter string total mass:");
            fscanf(infile,"%lf",&af_params->total_string_mass);
            (void) printf("%f\n",af_params->total_string_mass);
        }
        else
        {
            CursorAfterString(infile,"Enter string point mass:");
            fscanf(infile,"%lf",&af_params->m_l);
            (void) printf("%f\n",af_params->m_l);
        }
	}

	if (dim == 3 && af_params->is_parachute_system == YES)
	{
	    af_params->m_g = af_params->m_s;
            if (af_params->attach_gores == YES)
	    {
		CursorAfterString(infile,"Enter gore spring constant:");
        	fscanf(infile,"%lf",&af_params->kg);
        	(void) printf("%f\n",af_params->kg);
        	CursorAfterString(infile,"Enter gore damping constant:");
        	fscanf(infile,"%lf",&af_params->lambda_g);
        	(void) printf("%f\n",af_params->lambda_g);
		if (af_params->use_total_mass)
		{
		    CursorAfterString(infile,"Enter gore total mass:");
		    fscanf(infile,"%lf",&af_params->total_gore_mass);
		    (void) printf("%f\n",af_params->total_gore_mass);
		}
		else
		{
		    CursorAfterString(infile,"Enter gore point mass:");
		    fscanf(infile,"%lf",&af_params->m_g);
		    (void) printf("%f\n",af_params->m_g);
		}
            }
	}


        if (af_params->use_total_mass)
            convert_to_point_mass(front,af_params);
	if (af_params->is_parachute_system == NO)
	{
	    if (af_params->m_s == 0)
		af_params->m_s = af_params->m_l;
	    if (af_params->m_l == 0)
		af_params->m_l = af_params->m_s;
	}

        if (dim == 3)
        {
	    printf("canopy points count (fabric+gore)  = %d, "
		    "string points count = %d\n",
	    countSurfPoints(front->interf), 
	    countStringPoints(front->interf,af_params->is_parachute_system));
	    printf("fabric point mass  = %f,\nstring point mass = %f,\n"
		    "gore point mass = %f\n",af_params->m_s,af_params->m_l, 
		    af_params->m_g);
        }

	af_params->num_smooth_layers = 1;
	if (CursorAfterStringOpt(infile,"Enter number of smooth layers:"))
	{
            fscanf(infile,"%d",&af_params->num_smooth_layers);
            (void) printf("%d\n",af_params->num_smooth_layers);
	}
        af_params->payload = af_params->m_s;       // default
	if (CursorAfterStringOpt(infile,"Enter payload:"))
	{
            fscanf(infile,"%lf",&af_params->payload);
            (void) printf("%f\n",af_params->payload);
	}

    //TODO: this can be start of an InitAABBParams function
    if (CursorAfterStringOpt(infile,"Enter fabric thickness:"))
    {
        fscanf(infile,"%lf",&af_params->fabric_thickness);
        (void) printf("%f\n",af_params->fabric_thickness);
    }
    
    if (CursorAfterStringOpt(infile,"Enter vol_diff coefficient:"))
    {
        fscanf(infile,"%lf",&af_params->vol_diff);
        (void) printf("%f\n",af_params->vol_diff);
    }

	fclose(infile);
}	/* end SMM_InitSpringMassParams */

extern void SMM_InitTestVelFunc()
{
        Front *front = SMM_GetFront();
        AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	FILE *infile = fopen(InName(front),"r");
        af_params->no_fluid = YES;
        printf("\n");
        CursorAfterString(infile,"Start test velocity function set");
        printf("\n");
        initVelocityFunc(infile,front);
        fclose(infile);
}       /* end SMM_InitTestVelFunc */

extern void SMM_InitTestTimeControl()
{
        Front *front = SMM_GetFront();
	FILE *infile = fopen(InName(front),"r");
        printf("\n");
        CursorAfterString(infile,"Start test time control parameters");
        printf("\n");
        fclose(infile);
        FT_ReadTimeControl(InName(front),front);
}       /* end SMM_InitTestTimeContrl */

extern void SMM_Driver()
{
    SMM_StartUpStep();
    SMM_TimeMarch();
}

extern void SMM_TestDriver()
{
        Front *front = SMM_GetFront();
        F_BASIC_DATA *f_basic = SMM_GetBasicData();
	FILE *infile = fopen(InName(front),"r");
        double CFL;
        int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

        if (CursorAfterStringOpt(infile,"Turn off test run"))
            return;

        CFL = Time_step_factor(front);

	    SMM_Save();
        SMM_Plot();

        if (!f_basic->RestartRun)
        {
	    FT_Propagate(front);
            FT_RelinkGlobalIndex(front);
	    FT_InteriorPropagate(front);

            FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);

            FT_TimeControlFilter(front);
        }
        else
        {
            setSpecialNodeForce(front,af_params->kl);
            FT_SetOutputCounter(front);
	    FT_TimeControlFilter(front);
        }
        FT_PrintTimeStamp(front);
	
        // For restart debugging 
	if (FT_TimeLimitReached(front) && debugging("restart")) 
	{
	    SMM_Save();
	    return;
	}
    
        for (;;)
        {
            /* Propagating interface for time step dt */
            if(debugging("CLOCK"))
                reset_clock();

	    start_clock("timeStep");

            FrontPreAdvance(front);
            FT_Propagate(front);
            FT_RelinkGlobalIndex(front);
            FT_InteriorPropagate(front);

            if (debugging("trace"))
            {
                (void) printf("After solve()\n");
                (void) print_storage("at end of time step","trace");
            }

            // TMP
            (void) print_storage("at end of time step","trace");

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
                SMM_Save();
        
            if (FT_IsDrawTime(front))
                SMM_Plot();

            if (FT_TimeLimitReached(front))
	    {
                FT_PrintTimeStamp(front);
                if (!FT_IsSaveTime(front))
                    SMM_Save();
                if (!FT_IsDrawTime(front))
                    SMM_Plot();
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

extern void SMM_Plot()
{
        Front *front = SMM_GetFront();
        FT_Draw(front);
}       /* end SMM_Plot */

extern void SMM_Save()
{
        Front *front = SMM_GetFront();
        FT_Save(front);
        printAfExtraData(front,OutName(front));
}       /* end SMM_Save */

extern void SMM_CleanUp(int exitcode)
{
    clean_up(exitcode);
}

#ifdef __cplusplus
}
#endif
