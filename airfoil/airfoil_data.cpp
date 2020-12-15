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
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc before FT_StartUp()      
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    if (debugging("trace")) printf("Passed PetscInitialize()\n");

	Incompress_Solver_Smooth_Basis *l_cartesian = nullptr;
    if (f_basic.dim == 3)
        l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);
    else
    {
        printf("dim must == 3\n");
        clean_up(EXIT_FAILURE);
    }
	
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

    iFparams.dim = f_basic.dim;
    front.extra1 = (POINTER)&iFparams;
    front.extra2 = (POINTER)&af_params;
    read_iFparams(in_name,&iFparams);
    
    if (debugging("trace")) 
        (void) printf("Passed read_iFparams()\n");

    //2d initialization using old method 
    if (FT_Dimension() == 2)
        setInitialIntfcAF(&front,&level_func_pack,in_name);
    else
        level_func_pack.pos_component = LIQUID_COMP2;
       
    //TODO: to be removed along with initRigidBody_OLD()
	    //setInitialIntfcAF(&front,&level_func_pack,in_name);
	if (!RestartRun)
	{
	    FT_InitIntfc(&front,&level_func_pack);
	    if (f_basic.dim == 3)
            initIsolated3dCurves(&front);
	    
        initRigidBody(&front);
	    setRigidBodyMotionParams(&front,&rgb_params);
	    
        if (f_basic.dim == 3 && debugging("trace"))
	    {
            char gvdir[100];
            sprintf(gvdir,"%s/gv-init",out_name);
            gview_plot_interface(gvdir,front.interf);
	    }

	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim < 3)
            FT_ClipIntfcToSubdomain(&front);
	}
	else
	{
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	}

	initMovieStress(in_name,&front);

	/* Time control */
	FT_ReadTimeControl(in_name,&front);
        
	if (!RestartRun)
	{
	    optimizeElasticMesh(&front);
	    set_equilibrium_mesh(&front);
	    FT_SetGlobalIndex(&front);
	}

	/* Initialize velocity field function */
	setMotionParams(&front);

	front._compute_force_and_torque = ifluid_compute_force_and_torque;
	l_cartesian->findStateAtCrossing = af_find_state_at_crossing;
	l_cartesian->getInitialState = zero_state;
	l_cartesian->initMesh();
        
    if (RestartRun)
	{
	    if (ReSetTime) 
	    {
		    /* forbidden if restart with inherited states */
	    	readAfExtraData(&front,restart_state_name);
	    	modifyInitialization(&front);
	    	read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
		    l_cartesian->initMesh();
            l_cartesian->setInitialCondition();
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

	if (debugging("sample_velocity"))
        l_cartesian->initSampleVelocity(in_name);
    
    static_mesh(front.interf) = YES;

    l_cartesian->initMovieVariables();
    initMovieStress(in_name,&front);

    //TODO: NEED TO ZERO OUT OTHER FIELDS IN resetFronVelocity()?
	if (!RestartRun || ReSetTime)
	    resetFrontVelocity(&front);

    if (debugging("trace"))
        (void) printf("Passed state initialization()\n");

	/* Propagate the front */

	airfoil_driver(&front,l_cartesian);

	clean_up(0);
}

void airfoil_driver(Front *front,
        Incompress_Solver_Smooth_Basis *l_cartesian)
{
    double CFL;
    int  dim = front->rect_grid->dim;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

    CFL = Time_step_factor(front);
	Tracking_algorithm(front) = STRUCTURE_TRACKING;

	(void) printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun || ReSetTime)
	{
	    FT_ResetTime(front);

	    if (dim == 2)
	    {
	    	xgraph_front(front,out_name);
	    }

	    if (debugging("trace"))
            (void) printf("Calling FT_Save()\n");

	    setStressColor(front);
	    FT_Save(front);

       if (debugging("trace"))
           (void) printf("Calling printFrontInteriorStates()\n");
       
       l_cartesian->printFrontInteriorStates(out_name);
       printAfExtraData(front,out_name);

       if (debugging("trace"))
           (void) printf("Calling FT_Draw()\n");
       FT_Draw(front);

	    FrontPreAdvance(front);
	    FT_Propagate(front);
        FT_RelinkGlobalIndex(front);

	    if (!af_params->no_fluid)
	    {
            if (debugging("trace")) printf("Calling ifluid solve()\n");
            l_cartesian->solve(front->dt);
            if (debugging("trace")) printf("Passed ifluid solve()\n");
	    }
	    print_airfoil_stat(front,out_name);

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
	    FT_SetOutputCounter(front);
	}
	
    FT_TimeControlFilter(front);
	    //FT_PrintTimeStamp(front);

    
    //For generating PINN training data
    ///////////////////////////////////////////////////////
    l_cartesian->writeMeshFile3d();

    char tv_name[100];
    sprintf(tv_name,"%s/time.txt",out_name);

    auto imax = l_cartesian->getMaxIJK();

    char velm_name[100];
    sprintf(velm_name,"%s/velmat-%d-%d-%d.txt",
                out_name,imax[0],imax[1],imax[2]);
    
    char vortm_name[100];
    sprintf(vortm_name,"%s/vortmat-%d-%d-%d.txt",
                out_name,imax[0],imax[1],imax[2]);
    
    char intfc_name[100];
    sprintf(intfc_name,"%s/posintfc.txt",out_name);

    char intfc_disp_name[100];
    sprintf(intfc_disp_name,"%s/dispintfc.txt",out_name);
    
    char veli_name[100];
    sprintf(veli_name,"%s/velintfc.txt",out_name);
    
    /*
    char vorti_name[100];
    sprintf(vorti_name,"%s/vortintfc.txt",out_name);
    */
    
    
    IDATA3d idata0;
    int intfc_data_size;
    
    int tslice = 0;
    ///////////////////////////////////////////////////////

    for (;;)
    {
        //For generating PINN training data
        ///////////////////////////////////////////////////////
        if (FT_IsDrawTime(front) || tslice == 0)
        {
            auto vdata = l_cartesian->getVelData3d();

            FILE* tv_file = fopen(tv_name,"a");
            fprintf(tv_file,"%20.14f %20.14f\n",vdata.time,vdata.dt);
            fclose(tv_file);

            FILE* velm_file = fopen(velm_name,"a");
            FILE* vortm_file = fopen(vortm_name,"a");
            
            for (auto it : vdata.data)
            {
                //auto ic = it.icoords;
                
                auto iv = it.vel;
                fprintf(velm_file,"%20.14f %20.14f %20.14f\n",
                        iv[0],iv[1],iv[2]);

                /*
                auto vort = it.vort;
                fprintf(vortm_file,"%20.14f %20.14f %20.14f\n",
                        vort[0],vort[1],vort[2]);
                */
            }
        
            fclose(velm_file);
            fclose(vortm_file);
            
            FILE* intfc_file = fopen(intfc_name,"a");
            FILE* veli_file = fopen(veli_name,"a");
            FILE* intfc_disp_file = fopen(intfc_disp_name,"a");
            //FILE* vorti_file = fopen(vorti_name,"a");

            if (tslice == 0)
            {
                idata0 = l_cartesian->getIntfcData3d();
                intfc_data_size = idata0.data.size();
            }

            auto idata = l_cartesian->getIntfcData3d();

            //TODO: This is just a sanity check for now.
            //      Fabric runs should never change number
            //      of mesh points -- To Be Removed. 
            if (idata.data.size() != intfc_data_size)
            {
                printf("ERROR: num interface points has changed\n");
                printf("initial num intfc points: %d\n",intfc_data_size); 
                printf("current num intfc points: %d\n",idata.data.size()); 
                LOC(); clean_up(EXIT_FAILURE);
            }

            auto ientries = idata.data;
            auto ientries0 = idata0.data;

            //for (auto it : idata.data)
            for (int i = 0; i < intfc_data_size; ++i)
            {
                //auto ic = it.coords;
                auto ic = ientries[i].coords;
                fprintf(intfc_file,"%20.14f %20.14f %20.14f\n",
                        ic[0],ic[1],ic[2]);
            
                auto ic0 = ientries0[i].coords;
                fprintf(intfc_disp_file,"%20.14f %20.14f %20.14f\n",
                        ic[0]-ic0[0], ic[1]-ic0[1], ic[2]-ic0[2]);
                
                //auto iv = it.vel;
                auto iv = ientries[i].vel;
                fprintf(veli_file,"%20.14f %20.14f %20.14f\n",
                        iv[0],iv[1],iv[2]);
                
                /*
                //auto vort = it.vort;
                auto vort = ientries[i].vort;
                fprintf(vorti_file,"%20.14f %20.14f %20.14f\n",
                        vort[0],vort[1],vort[2]);
                */
            }

            fclose(intfc_file);
            fclose(intfc_disp_file);
            fclose(veli_file);
            //fclose(vorti_file);

            tslice++;
        }
        ///////////////////////////////////////////////////////
        
        /* Propagating interface for time step dt */
	    if (debugging("CLOCK"))
            reset_clock();

	    start_clock("time_step");
	    if (debugging("trace"))
            (void) printf("Before FT_Propagate()\n");

	    if (!af_params->no_fluid)
	    {
	    	coating_mono_hyper_surf(front);
	    	l_cartesian->applicationSetComponent();
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
	    	l_cartesian->applicationSetStates();
	    }
        
	    if (!af_params->no_fluid)
	    {
            if (debugging("trace")) printf("Calling ifluid solve()\n");
            l_cartesian->solve(front->dt);
            if (debugging("trace")) printf("Passed ifluid solve()\n");
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
            printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
        
        if (!af_params->no_fluid)
            front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
        
        if (debugging("step_size"))
            printf("Time step from l_cartesian->max_dt(): %f\n",front->dt);

	    /* Output section */

	    print_airfoil_stat(front,out_name);

        if (FT_IsSaveTime(front))
	    {
            setStressColor(front);
            FT_Save(front);
            l_cartesian->printFrontInteriorStates(out_name);
	    	printAfExtraData(front,out_name);
	    }
        if (debugging("trace"))
            (void) printf("After print output()\n");
        
        //TODO: suppress visualization data when, PINN training data working
        if (FT_IsDrawTime(front))
	    {
            FT_Draw(front);
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
	    stop_clock("time_step");
    }


    //For generating PINN training data
    ///////////////////////////////////////////////////////
    
    //Append final tslice to file names
    char new_tv_name[100];
    sprintf(new_tv_name,"%s/time-%d.txt",out_name,tslice);
    std::rename(tv_name,new_tv_name);

    char new_velm_name[100];
    sprintf(new_velm_name,"%s/velmat-%d-%d-%d-%d.txt",
                out_name,imax[0],imax[1],imax[2],tslice);
    std::rename(velm_name,new_velm_name);
    
    char new_vortm_name[100];
    sprintf(new_vortm_name,"%s/vortmat-%d-%d-%d-%d.txt",
                out_name,imax[0],imax[1],imax[2],tslice);
    std::rename(vortm_name,new_vortm_name);

    char new_intfc_name[100];
    sprintf(new_intfc_name,"%s/posintfc-%d-%d.txt",
                out_name,intfc_data_size,tslice);
    std::rename(intfc_name,new_intfc_name);
    
    char new_intfc_disp_name[100];
    sprintf(new_intfc_disp_name,"%s/dispintfc-%d-%d.txt",
            out_name,intfc_data_size,tslice);
    std::rename(intfc_disp_name,new_intfc_disp_name);
    
    char new_veli_name[100];
    sprintf(new_veli_name,"%s/velintfc-%d-%d.txt",
                out_name,intfc_data_size,tslice);
    std::rename(veli_name,new_veli_name);
    
    /*
    char new_vorti_name[100];
    sprintf(new_vorti_name,"%s/vortintfc-%d-%d.txt",
                out_name,intfc_data_size,tslice);
    std::rename(vorti_name,new_vorti_name);
    */
    ///////////////////////////////////////////////////////
    
    FT_FreeMainIntfc(front);
}       /* end airfoil_driver */


void zero_state(
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

