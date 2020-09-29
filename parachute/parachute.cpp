/****************************************************************
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
*****************************************************************/

/*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <iFluid.h>
#include <airfoil.h>

static void airfoil_driver(Front*,Incompress_Solver_Smooth_Basis*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void xgraph_front(Front*,char*);

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
	static IF_PARAMS iFparams;
	static AF_PARAMS af_params;
        static RG_PARAMS rgb_params;  

	FT_Init(argc,argv,&f_basic);
	f_basic.dim = 3;
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc before FT_StartUp()
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	Incompress_Solver_Smooth_Basis *l_cartesian = NULL;
	l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

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

    iFparams.dim = f_basic.dim;
    front.extra1 = (POINTER)&iFparams;
    front.extra2 = (POINTER)&af_params;
    front.extra3 = (POINTER)&rgb_params;
    read_iFparams(in_name,&iFparams);
	initParachuteDefault(&front);

	level_func_pack.pos_component = LIQUID_COMP2;
	if (!RestartRun)
	{
	    FT_InitIntfc(&front,&level_func_pack);
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    initParachuteModules(&front);
        
        if (consistent_interface(front.interf) == NO)
        {
            printf("consistent_interface(front.interf) == NO\n");
            clean_up(ERROR);
        }
	}
	else
	{
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	}

	FT_ReadTimeControl(in_name,&front);

    //FT_SetGlobalIndex() Must be called before setMotionParams()
    //in order to allow advanced restart scenarios
	if (!RestartRun)
    {
        optimizeElasticMesh(&front);
        set_equilibrium_mesh(&front);
        FT_SetGlobalIndex(&front);
        static_mesh(front.interf) = YES;
    }

    setMotionParams(&front);

	record_break_strings_gindex(&front);
	set_unequal_strings(&front);

	front._compute_force_and_torque = ifluid_compute_force_and_torque;
	l_cartesian->findStateAtCrossing = af_find_state_at_crossing;
	l_cartesian->getInitialState = zero_state;
    l_cartesian->initMesh();
	l_cartesian->skip_neumann_solver = YES;

	if (debugging("sample_velocity"))
        l_cartesian->initSampleVelocity(in_name);


    if (RestartRun)
    {
        if (ReSetTime)
        {
            readAfExtraData(&front,restart_state_name);
                //clearRegisteredPoints(&front);
                //resetRigidBodyVelocity(&front);
            modifyInitialization(&front);
                //rgb_init(&front,&rgb_params);
            read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
            l_cartesian->initMesh(); //TODO: may be able to remove this one
            l_cartesian->setInitialCondition();
        
            if (debugging("trace"))
            {
                if (consistent_interface(front.interf) == NO)
                    clean_up(ERROR);
                gview_plot_interface("gmodified",front.interf);
            }
        }
        else
        {
            l_cartesian->readFrontInteriorStates(restart_state_name);
            readAfExtraData(&front,restart_state_name);
        }

    }
    else
    {
        l_cartesian->setInitialCondition();
    }

	    //setMotionParams(&front);

    l_cartesian->initMovieVariables();
	initMovieStress(in_name,&front);
	    
	if (!RestartRun || ReSetTime)
        resetFrontVelocity(&front);

	/* Propagate the front */

	airfoil_driver(&front,l_cartesian);
	clean_up(0);
}

void airfoil_driver(Front *front,
        Incompress_Solver_Smooth_Basis *l_cartesian)
{
    int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    VPARAMS vort_params;
    bool inject_vortex = false;
    double start_time = HUGE;

    double CFL = Time_step_factor(front);
	Tracking_algorithm(front) = STRUCTURE_TRACKING;

        FILE *infile = fopen(InName(front),"r");
        if (CursorAfterStringOpt(infile,"Type yes to add vortex disturbance:"))
        {
            char answer[100];
            fscanf(infile,"%s",answer);
            printf("%s\n",answer);
            if (answer[0] == 'y')
            {
                inject_vortex = true;
                CursorAfterString(infile,"Enter center of vortex: ");
                fscanf(infile,"%lf %lf %lf",vort_params.center,vort_params.center+1,
                                    vort_params.center+2);
                printf("%f %f %f\n",vort_params.center[0],vort_params.center[1],
                                    vort_params.center[2]);
                CursorAfterString(infile,"Enter size of vortex: ");
                fscanf(infile,"%lf",&vort_params.D);
                printf("%f\n",vort_params.D);
                CursorAfterString(infile,"Enter intensity of vortex: ");
                fscanf(infile,"%lf",&vort_params.A);
                printf("%f\n",vort_params.A);
                CursorAfterString(infile,"Enter vortex start time: ");
                fscanf(infile,"%lf",&start_time);
                printf("%f\n",start_time);
            }
        }
        fclose(infile);

	if (!RestartRun || ReSetTime)
	{
	    FT_ResetTime(front);
	    setStressColor(front);
	    FT_Save(front);

        l_cartesian->printFrontInteriorStates(out_name);
	    printAfExtraData(front,out_name);
        FT_Draw(front);

	    FrontPreAdvance(front);
	    FT_Propagate(front);
        FT_RelinkGlobalIndex(front);
	    
        if (!af_params->no_fluid)
	    {
            l_cartesian->solve(front->dt);
	    }
	    print_airfoil_stat(front,out_name);
	    
        if (ReSetTime)
            setSpecialNodeForce(front,af_params->kl);

        FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);
        if (!af_params->no_fluid)
        {
	        front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
	        front->dt = std::min(front->dt,springCharTimeStep(front));
        }
	}
	else
	{
	    setSpecialNodeForce(front,af_params->kl);
	    FT_SetOutputCounter(front);
	}

    FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);

        // For restart debugging 
        if (FT_TimeLimitReached(front) && debugging("restart")) 
        {
            FT_Save(front);
            l_cartesian->printFrontInteriorStates(out_name);
            printAfExtraData(front,out_name);
            return;
        }
    
    for (;;)
    {
        /* Propagating interface for time step dt */
        if(debugging("CLOCK"))
            reset_clock();

        start_clock("timeStep");
        if (!af_params->no_fluid)
        {
            coating_mono_hyper_surf(front);
            l_cartesian->applicationSetComponent();
        }

        break_strings(front); // test

        FrontPreAdvance(front);
        FT_Propagate(front);
        FT_RelinkGlobalIndex(front);
        FT_InteriorPropagate(front);

        if (!af_params->no_fluid)
        {
            coating_mono_hyper_surf(front);
            l_cartesian->applicationSetStates();
        }

        if (!af_params->no_fluid)
        {
            if (inject_vortex && front->time >= start_time) 
            {
                inject_vortex = false;
                printf("Time to add vortex\n");
                l_cartesian->addVortexDisturbance(vort_params);
            }

            l_cartesian->solve(front->dt);
        }
        else
        {
            l_cartesian->max_dt = HUGE;
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
        {
            printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
        }
        
        if (!af_params->no_fluid)
            front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
   
        if (debugging("step_size"))
        {
            printf("Time step from l_cartesian->max_dt(): %f\n",front->dt);
        }

        /* Output section */

        print_airfoil_stat(front,out_name);

        if (!af_params->no_fluid)
            l_cartesian->printEnstrophy();

        if (FT_IsSaveTime(front))
        {
            setStressColor(front);
            FT_Save(front);
            l_cartesian->printFrontInteriorStates(out_name);
            printAfExtraData(front,out_name);
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
}       /* end airfoil_driver */


void    zero_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index, int dim,
        IF_PARAMS *iFparams)
{
        for (int i = 0; i < dim; ++i)
            field->vel[i][index] = 0.0;
        field->pres[index] = 0.0;
}


void xgraph_front(Front *front,	char *outname)
{
	char fname[100];
	sprintf(fname,"%s/intfc-%s",outname,right_flush(front->step,4));
	xgraph_2d_intfc(fname,front->interf);
}

