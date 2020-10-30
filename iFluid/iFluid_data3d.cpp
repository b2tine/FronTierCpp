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

#include "iFluid.h"
#include <cstdio>

	/*  Function Declarations */
static void ifluid_driver(Front*,Incompress_Solver_Smooth_Basis*);
static int l_cartesian_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	static IF_PARAMS iFparams;
    static RG_PARAMS rgb_params;
	IF_PROB_TYPE prob_type;

    //Read parameters from command line
	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);
	
    //Default is PETSC_COMM_WORLD = FronTier_COMM_WORLD = MPI_COMM_WORLD
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	/*Construct Incompress Solver l_cartesian*/
	Incompress_Solver_Smooth_Basis *l_cartesian = nullptr;
	if (f_basic.dim == 3)
	    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);
    else
    {
        printf("dim must == 3\n");
        clean_up(EXIT_FAILURE);
    }

    in_name                 = f_basic.in_name;
    restart_state_name      = f_basic.restart_state_name;
    out_name                = f_basic.out_name;
    restart_name            = f_basic.restart_name;
    RestartRun              = f_basic.RestartRun;
    RestartStep             = f_basic.RestartStep;

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

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	if (debugging("trace"))
        printf("Passed FT_StartUp()\n");

    iFparams.dim = f_basic.dim;
	front.extra1 = (POINTER)&iFparams;
	front.extra3 = (POINTER)&rgb_params;
	
    read_iF_prob_type(in_name,&prob_type);
	read_iFparams(in_name,&iFparams);

    if (debugging("trace"))
        printf("Passed read_iFparams()\n");

	/* Initialize interface through level function */

	setInitialIntfc(&front,&level_func_pack,in_name,prob_type);

	if (debugging("trace"))
        printf("Passed setInitialIntfc()\n");

	if (!RestartRun)
	{
	    if (f_basic.dim == 3)
            level_func_pack.set_3d_bdry = YES;
	    
        FT_InitIntfc(&front,&level_func_pack);
        initRigidBody(&front);
        setRigidBodyMotionParams(&front,&rgb_params);
	    FT_PromptSetMixedTypeBoundary2d(in_name,&front);

	    if (debugging("trace"))
	    {
            printf("Passed FT_InitIntfc()\n");

		    char test_name[100];
            switch (f_basic.dim)
            {
            case 2:
                sprintf(test_name,"%s/init_intfc-%d.xg",out_name,pp_mynode());
                xgraph_2d_intfc(test_name,front.interf);
                break;
            case 3:
                sprintf(test_name,"%s/gv-init",out_name);
                gview_plot_interface(test_name,front.interf);
                break;
            }
	    }

	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    
        if (f_basic.dim < 3)
	    	FT_ClipIntfcToSubdomain(&front);
	    
        if (debugging("trace")) 
            printf("Passed read_iF_dirichlet_bdry_data()\n");
	}
	else
    {
        read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
    }

	/* Initialize velocity field function */

    front._compute_force_and_torque = ifluid_compute_force_and_torque;
	velo_func_pack.func_params = (POINTER)l_cartesian;
	velo_func_pack.func = l_cartesian_vel;
	velo_func_pack.point_propagate = ifluid_point_propagate;
	
    FT_InitFrontVeloFunc(&front,&velo_func_pack);
	
    if (debugging("trace"))
	    printf("Passed FT_InitFrontVeloFunc()\n");

	l_cartesian->initMesh();
	l_cartesian->initMovieVariables();
	l_cartesian->findStateAtCrossing = ifluid_find_state_at_crossing;
	
    if (debugging("sample_velocity"))
	    l_cartesian->initSampleVelocity(in_name);

	init_fluid_state_func(l_cartesian,prob_type);
	if (debugging("trace"))
	    printf("Passed l_cartesian.initMesh()\n");

	if (!RestartRun)
	    l_cartesian->setInitialCondition();
	else
	    l_cartesian->readFrontInteriorStates(restart_state_name);
	
    if (debugging("trace"))
        printf("Passed state initialization()\n");

	if (iFparams.surf_tension != 0.0)
        front._contact_node_propagate = contact_node_propagate;

	/* Propagate the front */

	ifluid_driver(&front, l_cartesian);

	PetscFinalize();
	clean_up(0);
}

static void ifluid_driver(Front *front,
        Incompress_Solver_Smooth_Basis *l_cartesian)
{
	Curve_redistribution_function(front) = full_redistribute;
	FT_ReadTimeControl(in_name,front);
	double CFL = Time_step_factor(front);

	if (!RestartRun)
	{
	    FT_RedistMesh(front);
	}

    if (!RestartRun)
    {
        FT_ResetTime(front);
        FrontPreAdvance(front);
        
        if (debugging("trace"))
            printf("Before FT_Propagate() front->dt = %f\n",front->dt);
        
        FT_Propagate(front);
    
        if (debugging("trace"))
            printf("Calling ifluid solve()\n");
        
        l_cartesian->solve(front->dt);
    
        if (debugging("trace"))
            printf("Passed ifluid solve()\n");
    
        FT_SetTimeStep(front);
        front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
    
        FT_Draw(front);
        FT_Save(front);
        l_cartesian->printFrontInteriorStates(out_name);
        
        FT_SetOutputCounter(front);
    }
    else
    {
        FT_SetOutputCounter(front);
    }

	FT_TimeControlFilter(front);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("step_size"))
        printf("Time step from start: %f\n",front->dt);
    
    // Write mesh file
    l_cartesian->writeMeshFile3d();

    char tv_name[100];
    sprintf(tv_name,"%s/time.txt",out_name);

    char velm_name[100], vortm_name[100];
    auto imax = l_cartesian->getMaxIJK();
    sprintf(velm_name,"%s/velmat-%d-%d-%d",
                out_name,imax[0],imax[1],imax[2]);
    sprintf(vortm_name,"%s/vortmat-%d-%d-%d",
                out_name,imax[0],imax[1],imax[2]);
    if (pp_numnodes() > 1)
	{
        sprintf(velm_name,"%s-nd%s",velm_name,
                right_flush(pp_mynode(),4));
        sprintf(vortm_name,"%s-nd%s",vortm_name,
                right_flush(pp_mynode(),4));
	}
    sprintf(velm_name,"%s.txt",velm_name);
    sprintf(vortm_name,"%s.txt",vortm_name);


    int tslice = 0;
    for (;;)
    {
        //For generating PINN training data
            //if (FT_IsSaveTime(front) || FT_IsDrawTime(front) || tslice == 0)
        if (FT_IsDrawTime(front) || tslice == 0)
        {
            auto vdata = l_cartesian->getVelData3d();
            
            if (pp_mynode() == 0)
            {
                FILE* tv_file = fopen(tv_name,"a");
                fprintf(tv_file,"%20.14f %20.14f\n",
                        vdata.time,vdata.dt);
                fclose(tv_file);
            }

            FILE* velm_file = fopen(velm_name,"a");
            FILE* vortm_file = fopen(vortm_name,"a");
            
            for (auto it : vdata.data)
            {
                auto ic = it.icoords;
                
                auto iv = it.vel;
                fprintf(velm_file,"%20.14f %20.14f %20.14f\n",iv[0],iv[1],iv[2]);

                auto vort = it.vort;
                fprintf(vortm_file,"%20.14f %20.14f %20.14f\n",vort[0],vort[1],vort[2]);
            }
        
            fclose(velm_file);
            fclose(vortm_file);
            tslice++;
        }
        
        /* Propagating interface for time step dt */

        if (debugging("CLOCK"))
            reset_clock();

        if (debugging("trace"))
                printf("Before FT_Propagate()\n");
        
        FrontPreAdvance(front);
        FT_Propagate(front);
        setContactNodeType(front);
    
        if (debugging("trace"))
            printf("Passed FT_Propagate()\n");

        if (debugging("trace"))
            printf("Calling ifluid solve()\n");
        
        l_cartesian->solve(front->dt);
        
        if (debugging("trace"))
            printf("Passed ifluid solve()\n");
        
        if (debugging("trace"))
        {
            (void) printf("After solve()\n");
            (void) print_storage("at end of time step","trace");
        }

        FT_AddTimeStepToCounter(front);        
        FT_SetTimeStep(front);
        
        if (debugging("step_size"))
                (void) printf("Time step from FT_SetTimeStep(): %20.14f\n",
                    front->dt);
        
        front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
        
        if (debugging("step_size"))
                (void) printf("Time step from l_cartesian->max_dt(): %20.14f\n",
                    front->dt);

        /* Output section */

            //l_cartesian->printEnstrophy();

        if (FT_IsSaveTime(front))
        {
            FT_Save(front);
            l_cartesian->printFrontInteriorStates(out_name);
        }
        
        //Supress Visualization Data
        /*
        if (FT_IsDrawTime(front))
        {
            FT_Draw(front);
        }
        */
    
        //recordBdryEnergyFlux(front,out_name);

        if (FT_TimeLimitReached(front))
        {
            FT_PrintTimeStamp(front);
            
            if (debugging("CAUCHY_ERROR"))
                l_cartesian->compareWithBaseSoln();
            
            break;
        }

        if (debugging("storage"))
        {
            char s[100];
            sprintf(s,"Storage at end of time step %d",front->step);
            print_storage(s,"trace");
        }
        
        FT_TimeControlFilter(front);
        FT_PrintTimeStamp(front);
        
        if (debugging("step_size"))
            (void) printf("Time step from FT_TimeControlFilter(): %f\n",
                    front->dt);
    }// end time loop

    if (debugging("trace"))
        printf("After time loop\n");

    
    //Append final tslice to file names
    if (pp_mynode() == 0)
    {
        char new_tv_name[200];
        sprintf(new_tv_name,"%s/time-%d.txt",out_name,tslice);
        std::rename(tv_name,new_tv_name);
    }

    char new_velm_name[200];
    sprintf(new_velm_name,"%s/velmat-%d-%d-%d-%d",
                out_name,imax[0],imax[1],imax[2],tslice);
    
    char new_vortm_name[200];
    sprintf(new_vortm_name,"%s/vortmat-%d-%d-%d-%d",
                out_name,imax[0],imax[1],imax[2],tslice);
    
    if (pp_numnodes() > 1)
	{
        sprintf(new_velm_name,"%s-nd%s",new_velm_name,
                right_flush(pp_mynode(),4));
        sprintf(new_vortm_name,"%s-nd%s",new_vortm_name,
                right_flush(pp_mynode(),4));
	}
    sprintf(new_velm_name,"%s.txt",new_velm_name);
    sprintf(new_vortm_name,"%s.txt",new_vortm_name);

    std::rename(vortm_name,new_vortm_name);
    std::rename(velm_name,new_velm_name);

}       /* end ifluid_driver */

static int l_cartesian_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	double *coords = Coords(p);
	((Incompress_Solver_Smooth_Basis*)params)->getVelocity(coords, vel);
	return YES;
}	/* end l_cartesian_vel */
