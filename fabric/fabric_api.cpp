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

extern void SMM_InitCpp(int argc, char **argv)
{
        Front *front = SMM_GetFront();
        static F_BASIC_DATA f_basic;
        static LEVEL_FUNC_PACK level_func_pack;
        static AF_PARAMS af_params;

        FT_Init(argc,argv,&f_basic);
        front->extra2 = (POINTER)&af_params;
        f_basic.size_of_intfc_state = sizeof(STATE);
        af_params.node_id[0] = 0;

        FT_ReadSpaceDomain(f_basic.in_name,&f_basic);
        FT_StartUp(front,&f_basic);
        if (FT_Dimension() == 2) // initialization using old method
            setInitialIntfcAF(front,&level_func_pack,InName(front));
        FT_InitDebug(InName(front));
        FT_InitIntfc(front,&level_func_pack);
}       /* end SMM_InitCpp */

#ifdef __cplusplus
extern "C" {
#endif

extern void SMM_Init(char inname[])
{
        std::ifstream infile(inname);
        std::vector<std::string> argstrings;

        int argc = 0;
        while (!infile.eof())
        {
            std::string curr;
            infile >> curr;
            argstrings.push_back(curr);
            argc++;
        }

        argstrings.insert(argstrings.begin(),"dummyarg");
    
        char* argv[argc];
        for (int i = 0; i < argc; ++i)
        {
            argv[i] = new char[argstrings[i].length()+1];
            std::strcpy(argv[i],argstrings[i].c_str());
        }

        SMM_InitCpp(argc,argv);
}       /* end SMM_Init */

extern void SMM_InitModules()
{
        Front *front = SMM_GetFront();
        FILE *infile = fopen(InName(front),"r");

        if (FT_Dimension() == 3) // 2D initialization used old method
        {
            printf("\n");
            CursorAfterString(infile,"Start parameters for modules");
            printf("\n");
            initParachuteModules(front);
            initPerturbation3d(front);
            optimizeElasticMesh(front);
        }

        set_equilibrium_mesh(front);
        FT_SetGlobalIndex(front);

        fclose(infile);
}       /* end SMM_InitModules */

extern void SMM_InitPropagator()
{
        Front *front = SMM_GetFront();
        int dim = front->rect_grid->dim;

        Tracking_algorithm(front) = STRUCTURE_TRACKING;
        front->_point_propagate = elastic_point_propagate;
        front->curve_propagate = airfoil_curve_propagate;
        front->node_propagate = airfoil_node_propagate;
        front->interior_propagate = fourth_order_elastic_set_propagate;
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
        af_params->payload = 0.0;       // default
	if (CursorAfterStringOpt(infile,"Enter payload:"))
	{
            fscanf(infile,"%lf",&af_params->payload);
            (void) printf("%f\n",af_params->payload);
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

extern void SMM_InitTestTimeContrl()
{
        Front *front = SMM_GetFront();
	FILE *infile = fopen(InName(front),"r");
        printf("\n");
        CursorAfterString(infile,"Start test time control parameters");
        printf("\n");
        fclose(infile);
        FT_ReadTimeControl(InName(front),front);
}       /* end SMM_InitTestTimeContrl */

extern void SMM_TestDriver()
{
        Front *front = SMM_GetFront();
	FILE *infile = fopen(InName(front),"r");
        double CFL;
        int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

        if (CursorAfterStringOpt(infile,"Turn off test run"))
            return;

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
extern void SMM_Plot()
{
        Front *front = SMM_GetFront();
        FT_Draw(front);
}       /* end SMM_Plot */

#ifdef __cplusplus
}
#endif
