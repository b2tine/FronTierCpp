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


#include "cFluid.h"

#include <functional>
#include <chrono>
#include <random>


/*  Function Declarations */
static void promptForDirichletBdryState(FILE*,Front*,HYPER_SURF**,int,int);

static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);

static void rgbody_point_propagate(Front*,POINTER,POINT*,POINT*,
        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void rgbody_point_propagate_in_vacuum(Front*,POINTER,POINT*,POINT*,
        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void rgbody_point_propagate_in_fluid(Front*,POINTER,POINT*,POINT*,
        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);

static void cfluid_compute_force_and_torque2d(Front*,HYPER_SURF*,double,
                        double*,double*);
static void cfluid_compute_force_and_torque3d(Front*,HYPER_SURF*,double,
                        double*,double*);
static boolean force_on_hse(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,double*,
                                        double*,double*,boolean);
static boolean force_on_hse2d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static boolean force_on_hse3d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static double intrp_between(double,double,double,double,double);

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};

static void set_state_max_speed(Front*,STATE*,double*);
static void get_variable_bdry_params(int,FILE*,POINTER*);
static void cF_variableBoundaryState2d(double*,HYPER_SURF*,Front*,
					POINTER,POINTER);
static void cF_variableBoundaryState3d(double*,HYPER_SURF*,Front*,
					POINTER,POINTER);
/* test of open boundary */
static void pipe_end_func(Front*,POINTER,int*,COMPONENT,
				int,int,int*,Locstate);

static void get_turbulent_inlet_bdry_params(int dim, FILE *infile, POINTER *func_params);


extern void read_dirichlet_bdry_data(
	char *inname,
	Front *front)
{
	char msg[100];
	int i,j,k,nhs,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	HYPER_SURF *hs,**hss;
	INTERFACE *intfc = front->interf;
	int i_hs = 0;

	for (i = 0; i < dim; ++i)
	for (j = 0; j < 2; ++j)
	{
	    if (rect_boundary_type(intfc,i,j) == DIRICHLET_BOUNDARY)
	    {
            hs = FT_RectBoundaryHypSurf(intfc,DIRICHLET_BOUNDARY,i,j);
            if (hs == NULL)
            {
                printf("ERROR: cannot find Dirichlet boundary"
                   " in dimension %d direction %d\n",i,j);
                clean_up(ERROR);
            }
            if (j == 0)
                sprintf(msg,"For lower boundary in %d-th dimension",i);
            else
                sprintf(msg,"For upper boundary in %d-th dimension",i);
            CursorAfterString(infile,msg);
            (void) printf("\n");
            
            promptForDirichletBdryState(infile,front,&hs,1,i_hs);
            i_hs++;
	    }
	    else if (rect_boundary_type(intfc,i,j) == MIXED_TYPE_BOUNDARY)
	    {
            hss = FT_MixedBoundaryHypSurfs(intfc,i,j,DIRICHLET_BOUNDARY,&nhs);
            printf("Number of Dirichlet boundaries on dir %d side %d: %d\n",
                        i,j,nhs);
            if (dim == 2)
            {
                for (k = 0; k < nhs; ++k)
                {
                CURVE *c = Curve_of_hs(hss[k]);
                    (void) printf("Curve %d start and end at: ",k+1);
                    (void) printf("(%f %f)->(%f %f)\n",
                      Coords(c->start->posn)[0],
                      Coords(c->start->posn)[1],
                      Coords(c->end->posn)[0],
                      Coords(c->end->posn)[1]);
                promptForDirichletBdryState(infile,front,hss+k,1,i_hs);
                i_hs++;
                }
            }
	    }
	}

	hss = FT_InteriorHypSurfs(intfc,DIRICHLET_BOUNDARY,&nhs);
	if (nhs == 0 || hss == NULL) return;

	sprintf(msg,"For interior Dirichlet boundary:");
	CursorAfterString(infile,msg);
	(void) printf("\n");
	promptForDirichletBdryState(infile,front,hss,nhs,i_hs);
	i_hs++;
}	/* end read_dirichlet_bdry_data */

static void promptForDirichletBdryState(
	FILE *infile,
	Front *front,
	HYPER_SURF **hs,
	int nhs,
	int i_hs)
{
	static STATE *state;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	char s[100];
	COMPONENT comp;
	int dim = front->rect_grid->dim;
	POINTER func_params;

	FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
	state->dim = dim;

	CursorAfterString(infile,"Enter type of Dirichlet boundary:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	switch (s[0])
    {
        case 'c':			// Constant state
        case 'C':
        {
            comp = gas_comp(positive_component(hs[0])) ? 
                    positive_component(hs[0]) : negative_component(hs[0]);
            
            state->eos = &(eqn_params->eos[comp]);

            bool set_inlet_mach_number = false;
            if (CursorAfterStringOpt(infile,"Enter yes to set mach number at inlet:"))
            {
                char string[25];
                fscanf(infile,"%s",string);
                (void) printf("%s\n",string);
                if (string[0] == 'y' || string[0] == 'Y')
                {
                    set_inlet_mach_number = true;
                }
            }

            if (set_inlet_mach_number) //TODO: NOT CORRECT
            {
                /*
                CursorAfterString(infile,"Enter pressure:");
                fscanf(infile,"%lf",&state->pres);
                (void) printf("%f\n",state->pres);
                */

                CursorAfterString(infile,"Enter Mach number of shock:");
                fscanf(infile,"%lf",&eqn_params->Mach_number);
                (void) printf("%f\n",eqn_params->Mach_number);

                CursorAfterString(infile,"Enter idir of shock:");
                fscanf(infile,"%d",&eqn_params->idir);
                (void) printf("%d\n",eqn_params->idir);

                CursorAfterString(infile,"Enter direction of shock:");
                fscanf(infile,"%d",&eqn_params->shock_side);
                (void) printf("%d\n",eqn_params->shock_side);

                getChannelInletState(state,eqn_params,comp);

                printf("\nShock Speed: %f\n\n",eqn_params->shock_speed);    
            }
            else
            {
                CursorAfterString(infile,"Enter velocity:");
                for (int k = 0; k < dim; ++k)
                {
                    fscanf(infile,"%lf",&state->vel[k]);
                    (void) printf("%f ",state->vel[k]);
                }
                (void) printf("\n");

                CursorAfterString(infile,"Enter pressure:");
                fscanf(infile,"%lf",&state->pres);
                (void) printf("%f\n",state->pres);

                //Can specify temperature or density at the inlet.
                // Computes the unspecified value using the eqn of state.
                // Uses temperature as input if both values are specified.
                if (CursorAfterStringOpt(infile,"Enter temperature:"))
                {
                    fscanf(infile,"%lf",&state->temp);
                    printf("%f\n",state->temp);
                    state->dens = EosDensity(state);
                }
                else
                {
                    CursorAfterString(infile,"Enter density:");
                    fscanf(infile,"%lf",&state->dens);
                    printf("%f\n",state->dens);
                    state->temp = EosTemperature(state);
                }

                for (int k = 0; k < dim; ++k)
                {
                    state->momn[k] = state->dens*state->vel[k];
                }

                state->engy = EosEnergy(state);
                state->mu = EosViscosity(state);
            }
            
            ////////////////////////////////////
            state->k_turb = 0.0;
            ////////////////////////////////////


            //////////////////////////////////////////////////////////////////////////////////
            //FREESTREAM VALUES -- assumed to be that of inlet state
            //
            //TODO: Is this a correct/valid assumption?
            eqn_params->dens_freestream = state->dens;
            
            eqn_params->c_freestream = EosSoundSpeed(state);
            
            eqn_params->Mach_freestream = Magd(state->vel,dim)/eqn_params->c_freestream;
            
            eqn_params->pres_freestream = 
                eqn_params->dens_freestream * sqr(eqn_params->c_freestream) / (state->eos)->gamma;

            eqn_params->dir_freestream = (dim == 2) ? 0 : 2;
            if (CursorAfterStringOpt(infile,"Enter freestream direction:"))
            {
                fscanf(infile,"%d",&eqn_params->dir_freestream);
                printf("%d\n",eqn_params->dir_freestream);
            }

            eqn_params->vel_freestream = state->vel[eqn_params->dir_freestream];

            eqn_params->alpha = 0.0;
            eqn_params->beta = 0.0;
            
            //TODO: input options for alpha and beta (angle of attack and sideslip angle)
            
            //TODO: precompute U_dimless during initialization like the other values above.

            //TODO: make sure values are correct
            /*
            printf("\n");
            printf("Free-stream Density: %f\n",eqn_params->dens_freestream);
            printf("Free-stream Pressure: %f\n",eqn_params->pres_freestream);
            printf("Free-stream Sound Speed: %f\n",eqn_params->c_freestream);
            printf("Free-stream Mach Number: %f\n",eqn_params->Mach_freestream);
            printf("\n");
            */

            //////////////////////////////////////////////////////////////////////////////////

            
            printf("\n");
            printf("Inlet Density: %f\n",state->dens);
            printf("Inlet Velocity:");
            for (int i = 0; i < dim; ++i)
            {
                printf(" %f", state->vel[i]);
            }
            printf("\n");
            printf("Inlet Energy: %f\n",state->engy);
            printf("Inlet Pressure: %f\n",state->pres);
            printf("Inlet Temperature: %f\n",state->temp);
            printf("Inlet Viscosity: %g\n",state->mu);
            printf("\n");


            //void (*state_func)(double*,HYPER_SURF*,Front*,POINTER,POINTER) = nullptr;
            //POINTER state_func_params = nullptr;

            bool add_whitenoise = false;
            if (CursorAfterStringOpt(infile,"Enter yes to add white noise to inlet:"))
            {
                char string[25];
                fscanf(infile,"%s",string);
                (void) printf("%s\n",string);
                if (string[0] == 'y' || string[0] == 'Y')
                {
                    add_whitenoise = true;
                }
            }

            if (add_whitenoise)
            {
                //WHITE_NOISE_PARAMS* bdry_params = getWhiteNoiseParams(dim,infile);
                static WHITE_NOISE_PARAMS bdry_params;
                bdry_params.dim = dim;
                bdry_params.amplitude = 175.0; //TODO: read from input file
                std::copy(state,state+1,&(bdry_params.base_state));
                
                //state_func = cF_constantWithWhiteNoise;
                //state_func_params = (POINTER)&bdry_params;
            
                FT_InsertDirichletBoundary(front,cF_constantWithWhiteNoise,
                        "cF_constantWithWhiteNoise",(POINTER)&bdry_params,
                        nullptr,*hs,i_hs);
                /*FT_InsertDirichletBoundary(front,state_func,"cF_constantWithWhiteNoise",
                        state_func_params,(POINTER)state,*hs,i_hs);*/
            }
            else
            {
                FT_InsertDirichletBoundary(front,NULL,NULL,NULL,(POINTER)state,*hs,i_hs);
            }

                //FT_InsertDirichletBoundary(front,NULL,NULL,NULL,(POINTER)state,*hs,i_hs);

            for (int i = 1; i < nhs; ++i)
            {
                bstate_index(hs[i]) = bstate_index(hs[0]);
            }
        }
        break;
        
        /*
        //TODO:
        case 't':			// Turbulent Inlet state
        case 'T':
            get_turbulent_inlet_bdry_params(dim,infile,&func_params);
            FT_InsertDirichletBoundary(front,cF_turbulentInletBoundaryState,
                "cF_turbulentInletBoundaryState",func_params,NULL,hs[0],i_hs);
            for (int i = 1; i < nhs; ++i)
            {
                bstate_index(hs[i]) = bstate_index(hs[0]);
            }
            break;
        */

        case 's':
        case 'S':			// Supersonic outflow
        {
            if (s[1] == 'U' || s[1] == 'u')
            {
                FT_InsertDirichletBoundary(front,cF_supersonicOutflowState,
                        "cF_supersonicOutflowState",NULL,NULL,*hs,i_hs);
                for (int i = 1; i < nhs; ++i)
                {
                    bstate_index(hs[i]) = bstate_index(hs[0]);
                }
            }
            else if (s[1] == 'Y' || s[1] == 'y')
            {
                FT_InsertDirichletBoundary(front,cF_symmetryBoundaryState,
                        "cF_symmetryBoundaryState",NULL,NULL,*hs,i_hs);
                for (int i = 1; i < nhs; ++i)
                {
                    bstate_index(hs[i]) = bstate_index(hs[0]);
                }
            }
            else
            {
                printf("ERROR\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
        }
        break;

        case 'f':
        case 'F':
        {
	        switch (s[1])
            {
                case 'a':
                case 'A':			// Far-field state
                {
                    FT_InsertDirichletBoundary(front,cF_farfieldBoundaryState,
                            "cF_farfieldBoundaryState",NULL,NULL,*hs,i_hs);
                    for (int i = 1; i < nhs; ++i)
                    {
                        bstate_index(hs[i]) = bstate_index(hs[0]);
                    }
                }
                break;

                case 'l':
                case 'L':			// Flow through state
                {
                    FT_InsertDirichletBoundary(front,cF_flowThroughBoundaryState,
                            "cF_flowThroughBoundaryState",NULL,NULL,*hs,i_hs);
                    for (int i = 1; i < nhs; ++i)
                    {
                        bstate_index(hs[i]) = bstate_index(hs[0]);
                    }
                }
                break;
            }
        }
        break;
        
        case 'v':			// Variable state
        case 'V':
        {
            get_variable_bdry_params(dim,infile,&func_params);
            FT_InsertDirichletBoundary(front,cF_variableBoundaryState,
                "cF_variableBoundaryState",func_params,NULL,hs[0],i_hs);
            for (int i = 1; i < nhs; ++i)
            {
                bstate_index(hs[i]) = bstate_index(hs[0]);
            }
        }
        break;
    }
} 	/* end  promptForDirichletBdryState */

extern void cF_constantWithWhiteNoise(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
    FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
    EQN_PARAMS *eqn_params = ft_params->eqn_params;
    COMPONENT comp = ft_params->comp;
    POINT *oldp = ft_params->oldp;
    HYPER_SURF *oldhs = oldp->hs;
    HYPER_SURF_ELEMENT *oldhse = oldp->hse;

    WHITE_NOISE_PARAMS* bdry_params = (WHITE_NOISE_PARAMS*)boundary_state_function_params(hs);
    int dim = bdry_params->dim;
    double amp = bdry_params->amplitude;

    STATE *newst = (STATE*)state;
    newst->eos = &eqn_params->eos[comp]; 
	
    for (int i = 0; i < MAXD; ++i)
    {
        newst->vel[i] = 0.0;
        newst->momn[i] = 0.0;
    }
    
    STATE* bstate = &(bdry_params->base_state);
    //if (boundary_state(oldhs) != NULL)
    if (bstate != nullptr)
	{
	    //STATE *bstate = (STATE*)boundary_state(oldhs);

	    newst->dens = bstate->dens;
        newst->pres = bstate->pres;
        newst->temp = bstate->temp;
        newst->mu = bstate->mu;
        
        ////////////////////////////////////////////////////////////////////
        double white_noise_vel[MAXD] = {0.0};
        for (int i = 0; i < dim; ++i)
        {
            white_noise_vel[i] = (bstate->vel[i] == 0) ? 0.0 : amp;
        }

        auto seed = static_cast<long unsigned int>(
                std::chrono::system_clock::now().time_since_epoch().count());
        
        auto gen_sign =
        std::bind(std::uniform_int_distribution<int>{-1,1},
                                    std::default_random_engine{seed});

        int sign = gen_sign();

        for (int i = 0; i < dim; ++i)
        {
            newst->vel[i] = bstate->vel[i] + sign*white_noise_vel[i];
        }
        ////////////////////////////////////////////////////////////////////

        for (int i = 0; i < dim; ++i)
        {
	    	newst->momn[i] = newst->dens*newst->vel[i];
        }
	    newst->engy = EosEnergy(newst);
	    
        //TODO: Should vort/vorticity be non-zero for turbulent inlet bdry?
        newst->vort = 0.0;

        set_state_max_speed(front,newst,p0);

	    if (debugging("const_whitenoise_bdry"))
	    {
            printf("\n");
            printf("Coords: %f %f %f\n", p0[0], p0[1], p0[2]);
            printf("Velocity: %f %f %f\n", newst->vel[0], newst->vel[1], newst->vel[2]);
            //print_general_vector("Coords: ",p0,dim,"\n");
            //print_general_vector("Velocity: ",newst->vel,dim,"\n");
            //print_general_vector("Momentum: ",newst->momn,dim,"\n");
            printf("Density: %f\n",newst->dens);
            printf("Energy: %f\n",newst->engy);
            printf("Pressure: %f\n",newst->pres);
            printf("Temperature: %f\n",newst->temp);
            printf("Viscosity: %f\n",newst->mu);
	    }
    }
    else
    {
        printf("\n\nERROR cF_constantWithWhiteNoise(): no boundary state found!\n");
        LOC(); clean_up(EXIT_FAILURE);
    }
}

extern void cF_variableBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    cF_variableBoundaryState2d(p0,hs,front,params,state);
	    return;
	case 3:
	    cF_variableBoundaryState3d(p0,hs,front,params,state);
	    return;
	}
}	/* end cF_variableBoundaryState */

extern void cF_variableBoundaryState2d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	VAR_BDRY_PARAMS *bdry_params;
	STATE *newst = (STATE*) state;
	int i=0, dim, nbr_pist;
	double *angles,half_angular_width,*center,jet_duration_time;
	double radius,theta,vec[MAXD];
	boolean within_piston = NO;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;

	bdry_params = (VAR_BDRY_PARAMS*)boundary_state_function_params(hs);
	jet_duration_time = bdry_params->jet_duration_time;
	if (front->time > jet_duration_time)
	{
	    cF_flowThroughBoundaryState(p0,hs,front,params,state);
	    return;
	}

	dim = bdry_params->dim;
	center = bdry_params->center;
	nbr_pist = bdry_params->number_pistons;
	half_angular_width = bdry_params->half_angular_width;
	angles = bdry_params->angles_pistons;

	radius = 0.0;
	for (i = 0; i < dim; ++i) 
	{
	    vec[i] = p0[i] - center[i];
	    radius += sqr(vec[i]);
	}
	radius = sqrt(radius);
	for (i = 0; i < dim; ++i) 
	    vec[i] /= -radius;
	theta = asin(fabs(p0[1] - center[1])/radius);
	if (p0[0]-center[0] < 0 && p0[1]-center[1] > 0)
            theta = PI - theta;
	else if (p0[0]-center[0] < 0 && p0[1]-center[1] < 0)
            theta = PI + theta;
	else if (p0[0]-center[0] > 0 && p0[1]-center[1] < 0)
            theta = 2*PI - theta;
	for (i = 0; i < nbr_pist; ++i)
	{
	    if (theta > angles[i] - half_angular_width &&
		theta < angles[i] + half_angular_width)
	    {
		within_piston = YES;
	    }
	}
	if (within_piston)
	{
	    POINT *oldp = ft_params->oldp;
	    HYPER_SURF *oldhs = oldp->hs;
	    HYPER_SURF_ELEMENT *oldhse = oldp->hse;
	    STATE *sl,*sr;
	    COMPONENT comp;
	    slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	    if (gas_comp(negative_component(oldhs)))
	    {
	        newst = (STATE*)state;
	        comp = negative_component(oldhs);
	        newst->eos = sl->eos;
	    }
	    else if (gas_comp(positive_component(oldhs)))
	    {
	        newst = (STATE*)state;
	        comp = positive_component(oldhs);
	        newst->eos = sr->eos;
	    }

	    newst->dens = bdry_params->bdry_dens;
	    newst->pres = bdry_params->bdry_pres;
	    for (i = 0; i < dim; ++i)
	    {
		newst->vel[i] = bdry_params->bdry_vel*vec[i];
		newst->momn[i] = (newst->dens)*(newst->vel[i]);
	    }
	    newst->engy = EosEnergy(newst);
	    set_state_max_speed(front,newst,p0);
	}
	else
	{
	    cF_flowThroughBoundaryState(p0,hs,front,params,state);
	}
}	/* end cF_variableBoundaryState2d */

void cF_variableBoundaryState3d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
    ////////////////////////////////////////////////////////////////////////
    printf("\nERROR cF_variableBoundaryState3d() not implemented yet!\n");
    LOC(); clean_up(EXIT_FAILURE);
    ////////////////////////////////////////////////////////////////////////
}	/* end cF_variableBoundaryState3d */

extern void cF_farfieldBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
    int dim = front->rect_grid->dim;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;

	STATE* newst = (STATE*)state;
    newst->eos = &eqn_params->eos[comp]; 

    //TODO: read these from input file; hardcoded now to prototype
    //
    //      SEE EQN_PARAMS for declarations
    
    int dir_free = eqn_params->dir_freestream;
    double alpha = eqn_params->alpha;
    double beta = eqn_params->beta;

    double M_free = eqn_params->Mach_freestream;
    double c_free = eqn_params->c_freestream;
    double dens_free = eqn_params->dens_freestream;
    double pres_free = eqn_params->pres_freestream;
    
    /*
    double M_free = 2.0;
    double alpha = 0.0; //angle of attack
    double beta = 0.0;  //yaw

    double c_freestream = 331.61;
    double dens_freestream = 1.83;
    
    double pres_freestream = dens_freestream * sqr(c_freestream) / (newst->eos)->gamma;
    */

    newst->dens = dens_free;
    newst->pres = pres_free;

    double U_dimless[3] = {0.0};

    //TODO: Can also precompute U_dimless during initialization like the other values
    
    /*
    U_dimless[0] = M_freestream * std::cos(alpha) * std::cos(beta);
    U_dimless[1] = -1.0*M_freestream * std::sin(beta);
    U_dimless[2] = M_freestream * std::sin(alpha) * std::cos(beta);
    */

    //TODO: Double check these for accuracy, and simplify if possible.
    switch (dir_free)
    {
        case 0: //alpha in x-z plane, beta in x-y plane
            U_dimless[0] = M_free* std::cos(alpha) * std::cos(beta);
            U_dimless[1] = -1.0*M_free* std::sin(beta);
            U_dimless[2] = M_free* std::sin(alpha) * std::cos(beta);
            break;
        case 1: //alpha in y-x plane, beta in y-z plane
            U_dimless[0] = M_free* std::sin(alpha) * std::cos(beta);
            U_dimless[1] = M_free* std::cos(alpha) * std::cos(beta);
            U_dimless[2] = -1.0*M_free* std::sin(beta);
            break;
        case 2: //alpha in z-y plane, beta in z-x plane
            U_dimless[0] = -1.0*M_free* std::sin(beta);
            U_dimless[1] = M_free* std::sin(alpha) * std::cos(beta);
            U_dimless[2] = M_free* std::cos(alpha) * std::cos(beta);
            break;
    }
    

    for (int i = 0; i < dim; ++i)
    {
        newst->vel[i] = U_dimless[i] * c_free;
        newst->momn[i] = newst->dens * newst->vel[i];
    }

    newst->engy = EosEnergy(newst);
    newst->temp = EosTemperature(newst);
    newst->mu = EosViscosity(newst);
    
    set_state_max_speed(front,newst,p0);
    
    if (debugging("far_field"))
	{
	    printf("far-field boundary state:\n");
        printf("Coords:");
        for (int i = 0; i < dim; ++i)
        {
            printf(" %f",Coords(oldp)[i]);
        }
        printf("\n\n");

	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
	    printf("Density: %f\n",newst->dens);
	    printf("Energy: %f\n",newst->engy);
	    printf("Temperature: %f\n",newst->temp);
	    printf("Viscosity: %f\n",newst->mu);
	}
}

/*
//OLD VERSION
extern void cF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;
	
    double dir[MAXD];
	double u[3];		// velocity in the sweeping direction
	double v[3][MAXD];	// velocity in the orthogonal direction
	double vort[3];		// vorticity stencil
	double pres[3];		// pressure stencil
	double dens[3];		// pressure stencil
	double f_u;		    // u flux in the sweeping direction
	double f_v[MAXD];	// v flux in the orthogonal direction
	double f_vort;		// vort flux
	double f_pres;		// pressure flux
	double f_dens;		// density flux
	
    double dn, dt = front->dt;
	STATE  *s0,*sl,*sr,**sts;
	int i,j,dim = front->rect_grid->dim;
	
    int nrad = 2;
	
    STATE *newst = (STATE*)state;
    newst->eos = &eqn_params->eos[comp]; 


	if (debugging("flow_through"))
	    printf("Entering cF_flowThroughBoundaryState()\n");
	
	static STATE *s1;
    if (s1 == NULL)
	{
	    FT_ScalarMemoryAlloc((POINTER*)&s1,sizeof(STATE));
    }
	 

	tsten = FrontGetTanStencils(front,oldp,nrad);
	if (tsten != NULL)
	{
	    for (i = 0; i < dim; ++i)
	    	dir[i] = tsten[0]->dir[i];
	    dn = FT_GridSizeInDir(dir,front);

	    if (comp == negative_component(hs))  
	    {
	        sts = (STATE**)tsten[0]->leftst;
		s0 = sts[0];
	    }
	    else 
	    {
	        sts = (STATE**)tsten[0]->rightst;
		s0 = sts[0];
	    }

	    if (debugging("flow_through"))
	    {
	    	(void) printf("Ambient component: %d\n",comp);
	    	(void) printf("hs = %p  oldp->hs = %p\n",
					(POINTER)hs,(POINTER)oldp->hs);
	    	(void) printf("Time step = %f  Tangential grid size = %f\n",
					dt,dn);
	    	(void) printf("Tangential direction: ");
	    	for (j = 0; j < dim; ++j)
		    (void) printf("%f ",tsten[0]->dir[j]);
	    	(void) printf("\n");
	    	(void) printf("Tan_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    	(void) printf("Left points:\n");
	    	for (i = 0; i < nrad; ++i)
	    	{
		    for (j = 0; j < dim; ++j)
	    	    	(void) printf("%f ",Coords(tsten[0]->p[-i])[j]);
		    (void) printf("\n");
	    	}
	    	(void) printf("Right points:\n");
	    	for (i = 0; i < nrad; ++i)
	    	{
		    for (j = 0; j < dim; ++j)
	    	    	(void) printf("%f ",Coords(tsten[0]->p[i])[j]);
		    (void) printf("\n");
	    	}
	    }

	    for (j = 0; j < 3; ++j)
	    	u[j] = 0.0;
	    
        for (j = 0; j < 3; ++j)
	    {
	    	vort[j] = sts[j-1]->vort;
	    	pres[j] = sts[j-1]->pres;
	    	dens[j] = sts[j-1]->dens;
	    	for (i = 0; i < dim; ++i)
	    	{
		    u[j] += sts[j-1]->vel[i]*dir[i];
		    v[j][i] = sts[j-1]->vel[i]*(1.0 - dir[i]);
	    	}
	    }

	    f_u = burger_flux(u[0],u[1],u[2]);
	    for (i = 0; i < dim; ++i)
	    	f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	    f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	    f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	    f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	    for (i = 0; i < dim; ++i)
	    	newst->vel[i] = sts[0]->vel[i] - dt/dn*(f_u*dir[i] + f_v[i]);
	    newst->vort = sts[0]->vort - dt/dn*f_vort;
	    newst->pres = sts[0]->pres - dt/dn*f_pres;
	    newst->dens = sts[0]->dens - dt/dn*f_dens;
	}
	else
	{
	    slsr(oldp,oldp->hse,oldp->hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (comp == negative_component(hs))  
		    s0 = sl;
	    else
		    s0 = sr;
	}
	
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

	if (debugging("flow_through"))
	{
	    printf("Time step = %f  Normal grid size = %f\n",dt,dn);
	    printf("Normal direction: ");
	    for (j = 0; j < dim; ++j)
		printf("%f ",nsten->nor[j]);
	    printf("\n");
	    printf("Nor_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    printf("Nor_stencil:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",nsten->pts[i][j]);
		printf("\n");
	    }
	}

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    vort[j] = s0->vort;
	    pres[j] = s0->pres;
	    dens[j] = s0->dens;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += s0->vel[i]*dir[i];
		v[j][i] = s0->vel[i]*(1.0 - dir[i]);
	    }
	}

	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vel[i],getStateVel[i],&vtmp,&s0->vel[i]);
	    s1->vel[i] = vtmp;
	}
	
    if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vort,getStateVort,&vort[2],&s0->vort);
	    s1->vort = vort[2];
	}
	
    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->pres,
                            getStatePres,&pres[2],&s0->pres);
	
    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->dens,
                            getStateDens,&dens[2],&s0->dens);
	
    s1->pres = pres[2];
	s1->dens = dens[2];
	
    for (i = 0; i < dim; ++i)
	{
	    u[2] += s1->vel[i]*dir[i];
	    v[2][i] = s1->vel[i] - s1->vel[i]*dir[i];
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]);
	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;
	newst->dens += - dt/dn*f_dens;
	
	for (i = 0; i < dim; ++i)
        newst->momn[i] = newst->dens * newst->vel[i];

    newst->temp = EosTemperature(newst);
    newst->engy = EosEnergy(newst);
    newst->mu = EosViscosity(newst);
    
    set_state_max_speed(front,newst,p0);
	
    if (debugging("flow_through"))
	{
	    printf("flow through boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
        printf("Vorticity: %f\n",newst->vort);
	    printf("Pressure: %f\n",newst->pres);
	    printf("Density: %f\n",newst->dens);
	    printf("Energy: %f\n",newst->engy);
	    printf("Temperature: %f\n",newst->temp);
	    printf("Viscosity: %f\n",newst->mu);
	}
}*/       /* end cF_flowThroughBoundaryState */


//NEW VERSION
extern void cF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
    switch (front->rect_grid->dim)
    {
    case 2:
        return cF_flowThroughBoundaryState2d(p0,hs,front,params,state);
    case 3:
        return cF_flowThroughBoundaryState3d(p0,hs,front,params,state);
    default:
        printf("\nERROR cF_flowThroughBoundaryState() : \
                unsupported spatial dimension! \n\tdim = %d\n",
                front->rect_grid->dim);
        LOC(); clean_up(EXIT_FAILURE);
    }
}

extern void cF_flowThroughBoundaryState2d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;
	
    double dir[MAXD];
	double u[3];		// velocity in the sweeping direction
	double v[3][MAXD];	// velocity in the orthogonal direction
	double vort[3];		// vorticity stencil
	double pres[3];		// pressure stencil
	
    double dens[3];		// pressure stencil
	double f_u;		    // u flux in the sweeping direction
	double f_v[MAXD];	// v flux in the orthogonal direction
	double f_vort;		// vort flux
	double f_pres;		// pressure flux
	double f_dens;		// density flux
	
    int dim = front->rect_grid->dim;
    double dn, dt = front->dt;
	int i,j;
	
    STATE* oldst;
    STATE* newst = (STATE*)state;
	STATE** sts;

    int nrad = 2;

	if (debugging("flow_through"))
	    printf("Entering cF_flowThroughBoundaryState2d()\n");
	
	POINTER sl, sr;
    FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
    
    if (comp == negative_component(hs))  
        oldst = (STATE*)sl;
    else
        oldst = (STATE*)sr;

    newst->eos = &eqn_params->eos[comp]; 
    for (i = 0; i < dim; ++i)
    {
        newst->vel[i] = oldst->vel[i];
    }
    newst->vort = oldst->vort;
    newst->pres = oldst->pres;
    newst->dens = oldst->dens;

    //Normal
	Nor_stencil* nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);

	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

    if (debugging("flow_through"))
    {
        (void) printf("Normal grid size = %f\n",dn);
        (void) print_Nor_stencil(front,nsten);
    }

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;

	for (j = 0; j < 2; ++j)
	{
        //normal component of velocity
	    for (i = 0; i < dim; ++i)
	    {
            u[j] += oldst->vel[i]*dir[i];
	    }

        //orthogonal direction velocity
	    for (i = 0; i < dim; ++i)
	    {
            v[j][i] = oldst->vel[i] - u[j]*dir[i];
	    }

	    vort[j] = oldst->vort;
	    pres[j] = oldst->pres;
	    dens[j] = oldst->dens;
	}

    STATE s1;
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
	    s1.vel[i] = vtmp;
	}

    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->vort,
            getStateVort,&vort[2],&oldst->vort);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->pres,
            getStatePres,&pres[2],&oldst->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->dens,
            getStateDens,&dens[2],&oldst->dens);

    s1.vort = vort[2];
	s1.pres = pres[2];
	s1.dens = dens[2];
	
    for (i = 0; i < dim; ++i)
	{
	    u[2] += s1.vel[i]*dir[i];
	}

    for (i = 0; i < dim; ++i)
	{
        v[2][i] = s1.vel[i] - u[2]*dir[i];
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
    {
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
    }
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	for (i = 0; i < dim; ++i)
    {
	    newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
    }
	newst->vort -= dt/dn*f_vort;
	newst->pres -= dt/dn*f_pres;
	newst->dens -= dt/dn*f_dens;

    if (debugging("flow_through"))
    {
        (void) print_Nor_stencil(front,nsten);
        (void) printf("new velocity after normal prop: %f %f %f\n",
            newst->vel[0],newst->vel[1],newst->vel[2]);
    }


    ////////////////////////////////////////////////////////////
    //Tangential
	Tan_stencil** tsten = FrontGetTanStencils(front,oldp,nrad);

	if (tsten != nullptr)
	{
	    for (i = 0; i < dim; ++i)
	    	dir[i] = tsten[0]->dir[i];
	    dn = FT_GridSizeInDir(dir,front);

        if (debugging("flow_through"))
        {
            (void) printf("Ambient component: %d\n",comp);
            (void) printf("Tangential grid size = %f\n",dn);
            (void) print_Tan_stencil(front,tsten[0]);
        }

        for (j = 0; j < 3; ++j)
            u[j] = 0.0;

	    if (comp == negative_component(hs))  
	        sts = (STATE**)tsten[0]->leftst;
	    else 
	        sts = (STATE**)tsten[0]->rightst;

        for (j = 0; j < 3; ++j)
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
		        u[j] += sts[j-1]->vel[i]*dir[i];
	    	}

	    	for (i = 0; i < dim; ++i)
	    	{
		        v[j][i] = sts[j-1]->vel[i] - u[j]*dir[i];
	    	}

	    	vort[j] = sts[j-1]->vort;
	    	pres[j] = sts[j-1]->pres;
	    	dens[j] = sts[j-1]->dens;
	    }

	    f_u = burger_flux(u[0],u[1],u[2]);
	    for (i = 0; i < dim; ++i)
        {
	    	f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
        }
	    f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	    f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	    f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	    for (i = 0; i < dim; ++i)
        {
	    	newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
        }

	    newst->vort -= dt/dn*f_vort;
	    newst->pres -= dt/dn*f_pres;
	    newst->dens -= dt/dn*f_dens;
	}
    ////////////////////////////////////////////////////////////

    for (i = 0; i < dim; ++i)
    {
        newst->momn[i] = newst->dens * newst->vel[i];
    }
    
    newst->temp = EosTemperature(newst);
    newst->engy = EosEnergy(newst);
    newst->mu = EosViscosity(newst);
    
    set_state_max_speed(front,newst,p0);
	
    if (debugging("flow_through"))
	{
	    printf("flow through boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
        printf("Vorticity: %f\n",newst->vort);
	    printf("Pressure: %f\n",newst->pres);
	    printf("Density: %f\n",newst->dens);
	    printf("Energy: %f\n",newst->engy);
	    printf("Temperature: %f\n",newst->temp);
	    printf("Viscosity: %f\n",newst->mu);
	}
}       /* end cF_flowThroughBoundaryState2d */

extern void cF_flowThroughBoundaryState3d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;
	
    double dir[MAXD];
	double u[3];		// velocity in the sweeping direction
	double v[3][MAXD];	// velocity in the orthogonal direction
	double vort[3];		// vorticity stencil
	double pres[3];		// pressure stencil
	double dens[3];		// pressure stencil
	
    double f_u;		    // u flux in the sweeping direction
	double f_v[MAXD];	// v flux in the orthogonal direction
	double f_vort;		// vort flux
	double f_pres;		// pressure flux
	double f_dens;		// density flux
	
	int i,j,dim = front->rect_grid->dim;
    double dn, dt = front->dt;
	
    STATE* oldst;
    STATE* newst = (STATE*)state;
	STATE** sts;
	
    int nrad = 2;

	if (debugging("flow_through"))
	    printf("Entering cF_flowThroughBoundaryState3d()\n");
	
	POINTER sl, sr;
    FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
    
    if (comp == negative_component(hs))  
        oldst = (STATE*)sl;
    else
        oldst = (STATE*)sr;

    newst->eos = &eqn_params->eos[comp]; 
    for (i = 0; i < dim; ++i)
    {
        newst->vel[i] = oldst->vel[i];
    }
    newst->pres = oldst->pres;
    newst->dens = oldst->dens;

    //Normal
	Nor_stencil* nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);

	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

    if (debugging("flow_through"))
    {
        (void) printf("Normal grid size = %f\n",dn);
        (void) print_Nor_stencil(front,nsten);
    }

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;

	for (j = 0; j < 2; ++j)
	{
	    for (i = 0; i < dim; ++i)
	    {
            u[j] += oldst->vel[i]*dir[i];
	    }

	    for (i = 0; i < dim; ++i)
	    {
            v[j][i] = oldst->vel[i] - u[j]*dir[i];
	    }

	    pres[j] = oldst->pres;
	    dens[j] = oldst->dens;
	}

    STATE s1;
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
	    s1.vel[i] = vtmp;
	}

	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->pres,
            getStatePres,&pres[2],&oldst->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->dens,
            getStateDens,&dens[2],&oldst->dens);

	s1.pres = pres[2];
	s1.dens = dens[2];
	
    for (i = 0; i < dim; ++i)
	{
	    u[2] += s1.vel[i]*dir[i];
	}

    for (i = 0; i < dim; ++i)
	{
        v[2][i] = s1.vel[i] - u[2]*dir[i];
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
    {
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
    }
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

	for (i = 0; i < dim; ++i)
    {
	    newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
    }
	newst->pres -= dt/dn*f_pres;
	newst->dens -= dt/dn*f_dens;

    if (debugging("flow_through"))
    {
        (void) print_Nor_stencil(front,nsten);
        (void) printf("new velocity after normal prop: %f %f %f\n",
            newst->vel[0],newst->vel[1],newst->vel[2]);
    }
    
    ////////////////////////////////////////////////////////////
    //Tangential
	Tan_stencil** tsten = FrontGetTanStencils(front,oldp,nrad);

	if (tsten != nullptr)
	{
        for (int k = 0; k < dim-1; ++k)
        {
            for (i = 0; i < dim; ++i)
                dir[i] = tsten[k]->dir[i];
            dn = FT_GridSizeInDir(dir,front);

            if (debugging("flow_through"))
            {
                (void) printf("Ambient component: %d\n",comp);
                (void) printf("Tangential grid size = %f\n",dn);
                (void) printf("For direction %d\n",k);
                (void) print_Tan_stencil(front,tsten[k]);
            }

            for (j = 0; j < 3; ++j)
                u[j] = 0.0;

            if (comp == negative_component(hs))  
                sts = (STATE**)tsten[k]->leftst;
            else 
                sts = (STATE**)tsten[k]->rightst;

            for (j = 0; j < 3; ++j)
            {
                for (i = 0; i < dim; ++i)
                {
                    u[j] += sts[j-1]->vel[i]*dir[i];
                }

                for (i = 0; i < dim; ++i)
                {
                    v[j][i] = sts[j-1]->vel[i] - u[j]*dir[i];
                }

                pres[j] = sts[j-1]->pres;
                dens[j] = sts[j-1]->dens;
            }

            f_u = burger_flux(u[0],u[1],u[2]);
            for (i = 0; i < dim; ++i)
            {
                f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
            }
            f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
            f_dens = linear_flux(u[1],dens[0],dens[1],dens[2]);

            for (i = 0; i < dim; ++i)
            {
                newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
            }

            newst->pres -= dt/dn*f_pres;
            newst->dens -= dt/dn*f_dens;
        }
	}
    ////////////////////////////////////////////////////////////
    
    for (i = 0; i < dim; ++i)
    {
        newst->momn[i] = newst->dens * newst->vel[i];
    }

    newst->temp = EosTemperature(newst);
    newst->engy = EosEnergy(newst);
    newst->mu = EosViscosity(newst);
	
    set_state_max_speed(front,newst,p0);
	
    if (debugging("flow_through"))
	{
	    printf("flow through boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
	    printf("Density: %f\n",newst->dens);
	    printf("Energy: %f\n",newst->engy);
	    printf("Temperature: %f\n",newst->temp);
	    printf("Viscosity: %f\n",newst->mu);
	}
}       /* end cF_flowThroughBoundaryState3d */

//TODO: Should only be used if local outflow mach number >= 1
extern void cF_supersonicOutflowState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;

	POINTER sl, sr;
    FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
    
    STATE* oldst;
    if (comp == negative_component(hs))  
        oldst = (STATE*)sl;
    else
        oldst = (STATE*)sr;

    STATE* newst = (STATE*)state;
    newst->eos = &eqn_params->eos[comp]; 
    
    int nrad = 2;
	Nor_stencil* nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
    double* nor = nsten->nor;

	int dim = front->rect_grid->dim;

	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->pres,
            getStatePres,&newst->pres,&oldst->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->dens,
            getStateDens,&newst->dens,&oldst->dens);
    
    double vel[MAXD] = {0.0};
	for (int i = 0; i < dim; ++i)
	{
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vel[i],getStateVel[i],&newst->vel[i],&oldst->vel[i]);
	    //FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
		//	eqn_params->vel[i],getStateVel[i],&vel[i],&oldst->vel[i]);
	}

    /*
    double vn = 0.0;
	for (int i = 0; i < dim; ++i)
	{
        vn += vel[i]*nor[i];
    }

	for (int i = 0; i < dim; ++i)
	{
        newst->vel[i] = vn*nor[i];
    }
    */

	for (int i = 0; i < dim; ++i)
	{
        newst->momn[i] = newst->dens*newst->vel[i];
    }

    newst->temp = EosTemperature(newst);
    newst->engy = EosEnergy(newst);
    newst->mu = EosViscosity(newst);
	
    set_state_max_speed(front,newst,p0);
	
    if (debugging("outflow"))
	{
	    printf("supersonic outflow boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
	    printf("Density: %f\n",newst->dens);
	    printf("Energy: %f\n",newst->engy);
	    printf("Temperature: %f\n",newst->temp);
	    printf("Viscosity: %f\n",newst->mu);
	}
}

extern void cF_symmetryBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	EQN_PARAMS *eqn_params = ft_params->eqn_params;

	POINTER sl, sr;
    FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
    
    STATE* oldst;
    if (comp == negative_component(hs))  
        oldst = (STATE*)sl;
    else
        oldst = (STATE*)sr;

    STATE* newst = (STATE*)state;
    newst->eos = &eqn_params->eos[comp]; 
    
    int nrad = 2;
	Nor_stencil* nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
    double* nor = nsten->nor;

	int dim = front->rect_grid->dim;

	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->pres,
            getStatePres,&newst->pres,&oldst->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->dens,
            getStateDens,&newst->dens,&oldst->dens);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->engy,
            getStateEngy,&newst->engy,&oldst->engy);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->temp,
            getStateTemp,&newst->temp,&oldst->temp);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],eqn_params->mu,
            getStateMu,&newst->mu,&oldst->mu);
    
    double vel[MAXD] = {0.0};
	for (int i = 0; i < dim; ++i)
	{
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			eqn_params->vel[i],getStateVel[i],&vel[i],&oldst->vel[i]);
	}

    double vn = 0.0;
	for (int i = 0; i < dim; ++i)
	{
        vn += vel[i]*nor[i];
    }

	for (int i = 0; i < dim; ++i)
	{
        newst->vel[i] = vel[i] - 2.0*vn*nor[i];
    }

	for (int i = 0; i < dim; ++i)
	{
        newst->momn[i] = newst->dens*newst->vel[i];
    }

    set_state_max_speed(front,newst,p0);
	
    if (debugging("symmetry_bdry"))
	{
	    printf("symmetry boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
	    printf("Density: %f\n",newst->dens);
	    printf("Energy: %f\n",newst->engy);
	    printf("Temperature: %f\n",newst->temp);
	    printf("Viscosity: %f\n",newst->mu);
	}
}

extern void cFluid_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	switch(wave_type(oldhs))
	{
        case MOVABLE_BODY_BOUNDARY:
            return rgbody_point_propagate(front,wave,oldp,newp,
                    oldhse,oldhs,dt,V);
        //case ELASTIC_BOUNDARY: //TODO: MAKE SURE THAT COMMENTING THIS OUT DOESN"T CREATE PROBLEMS
        case NEUMANN_BOUNDARY:
            return neumann_point_propagate(front,wave,oldp,newp,
                    oldhse,oldhs,dt,V);
        case DIRICHLET_BOUNDARY:
            return dirichlet_point_propagate(front,wave,oldp,newp,
                    oldhse,oldhs,dt,V);
        case SUBDOMAIN_BOUNDARY:
        case PASSIVE_BOUNDARY:
            return;
        default:
            return contact_point_propagate(front,wave,oldp,newp,
                    oldhse,oldhs,dt,V);
	}
}       /* cFluid_point_propagate */

static  void neumann_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	
    double *m_pres = eqn_params->pres;
	double *m_dens = eqn_params->dens;
	double *m_engy = eqn_params->engy;
	double *m_temp = eqn_params->temp;
	double *m_mu = eqn_params->mu;

    double *m_kturb = eqn_params->k_turb;

	double nor[MAXD],tan[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}
	FT_NormalAtPoint(oldp,front,nor,comp);
	newst->dim = dim;

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;
    
    //TODO: tangent vector not used,
    //      and this is in general not correct
	tan[0] = -nor[1];
    tan[1] = nor[0];

    //TODO:Can get rid of this -- handled by rgbody_point_propagate()
	if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
	{
        double omega_dt,crds_com[MAXD];
        omega_dt = angular_velo(oldhs)*dt;
        for (i = 0; i < dim; ++i)
        {
            vel[i] = center_of_mass_velo(oldhs)[i];
            crds_com[i] = Coords(oldp)[i] + dt*vel[i] - center_of_mass(oldhs)[i];
        }
    
        //if (dim == 2) //This is only for rotations in xy plane
        vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
                 angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
        vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
                 angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
	}
	else
	{
        for (i = 0; i < dim; ++i) vel[i] = 0.0;
	}

	for (i = 0; i < dim; ++i)
	{
        Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
	    newst->vel[i] = vel[i];
        FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
	}


    FT_IntrpStateVarAtCoords(front,comp,p1,m_pres,
			getStatePres,&newst->pres,&oldst->pres);
    
    FT_IntrpStateVarAtCoords(front,comp,p1,m_dens,
			getStateDens,&newst->dens,&oldst->dens);

    //TODO: should a wall model/relation be used for temperature here?
    //      e.g. coco-brusemann etc.
    //
	FT_IntrpStateVarAtCoords(front,comp,p1,m_temp,
			getStateTemp,&newst->temp,&oldst->temp);

	FT_IntrpStateVarAtCoords(front,comp,p1,m_mu,
			getStateMu,&newst->mu,&oldst->mu);
	
	FT_IntrpStateVarAtCoords(front,comp,p1,m_kturb,
			getStateKTurb,&newst->k_turb,&oldst->k_turb);
    
    for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	    newst->momn[i] = newst->dens*vel[i];
	}
    newst->eos = &(eqn_params->eos[comp]);
	
    newst->engy = EosEnergy(newst);
    //newst->temp = EosTemperature(newst);
    //newst->mu = EosViscosity(newst);

    s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
	set_state_max_speed(front,newst,Coords(newp));
}	/* end neumann_point_propagate */

static void dirichlet_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
    int dim = front->rect_grid->dim;
	STATE *sl,*sr,*newst = NULL;
	STATE *bstate;
	FLOW_THROUGH_PARAMS ft_params;
	COMPONENT comp;

	if (debugging("trace"))
	{
	    printf("Entering dirichlet_point_propagate()\n");
	}

	slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    newst = (STATE*)left_state(newp);
	    comp = negative_component(oldhs);
        newst->eos = &(eqn_params->eos[comp]);
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    newst = (STATE*)right_state(newp);
	    comp = positive_component(oldhs);
        newst->eos = &(eqn_params->eos[comp]);
	}
	if (newst == NULL) return;	// node point

	if (boundary_state_function(oldhs))
	{
	    oldp->hse = oldhse;
	    oldp->hs = oldhs;
	    ft_params.oldp = oldp;
	    ft_params.eqn_params = eqn_params;
	    ft_params.comp = comp;
	    (*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)&ft_params,(POINTER)newst);	
	}
    else if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);

	    newst->dens = bstate->dens;
        newst->pres = bstate->pres;
        
        for (int i = 0; i < dim; ++i)
        {
	    	newst->vel[i] = bstate->vel[i];
        }

        for (int i = 0; i < dim; ++i)
        {
	    	newst->momn[i] = newst->dens*newst->vel[i];
        }

	    newst->engy = bstate->engy;
	    newst->temp = bstate->temp;
	    newst->mu = bstate->mu;
	    
        //TODO: Should vort/vorticity be non-zero for turbulent inlet bdry?
        //      Needs to depend on the velocity of the boundary state,
        //      which in general we do not want to restrict to being
        //      a uniform/irrotational quantity.
        
        //////////////////////////////////////////////////////////////////
        newst->vort = 0.0;
        //////////////////////////////////////////////////////////////////
	}

    set_state_max_speed(front,newst,Coords(oldp));


    if (debugging("dirichlet_bdry"))
    {
        printf("\nDirichlet boundary state:\n");
        print_general_vector("Coords: ",Coords(oldp),dim,"\n\n");
        print_general_vector("Velocity: ",newst->vel,dim,"\n");
        printf("\tDensity: %f\n",newst->dens);
        printf("\tEnergy: %f\n",newst->engy);
        printf("\tPressure: %f\n",newst->pres);
        printf("\tTemperature: %f\n",newst->temp);
        printf("\tViscosity: %f\n",newst->mu);
        printf("\tVorticity: %f\n",newst->vort);
    }

	if (debugging("trace"))
	    printf("Leaving dirichlet_point_propagate()\n");
}	/* end dirichlet_point_propagate */

static  void contact_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
    int i, dim = front->rect_grid->dim;
	
    double **m_mom = eqn_params->mom;
	double *m_dens = eqn_params->dens;
	double *m_engy = eqn_params->engy;
	double *m_kturb = eqn_params->k_turb;

	double *p0;
	STATE *oldst,*newst;
	POINTER sl,sr;
	EOS_PARAMS *eos = eqn_params->eos;
	double default_var;

	switch (eqn_params->point_prop_scheme)
	{
	case FIRST_ORDER:
	    first_order_point_propagate(front,wave,oldp,newp,
				oldhse,oldhs,dt,V);
	    break;
	case SECOND_ORDER:
	    second_order_point_propagate(front,wave,oldp,newp,
				oldhse,oldhs,dt,V);
	    break;
	case FOURTH_ORDER:
	    fourth_order_point_propagate(front,wave,oldp,newp,
				oldhse,oldhs,dt,V);
	    break;
	}
	if (debugging("point_propagate"))
	{
	    double hmin,dist;
	    hmin = front->rect_grid->h[0];
	    dist = 0.0;
	    for (i = 0; i < dim; ++i)
		dist += sqr(Coords(newp)[i] - Coords(oldp)[i]);
	    dist = sqrt(dist);
	    if (dist > 0.5*hmin || isnan(dist))
	    {
		printf("WARNING: propagating over half grid size!\n");
		printf("dist/hmin = %f\n",dist/hmin);
		if (dist > hmin || isnan(dist))
		    clean_up(ERROR);
	    }
	}

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	p0 = Coords(newp);
	newst = (STATE*)left_state(newp);
	oldst = (STATE*)sl;
	*newst = *oldst;
	newst->dim = dim;
	newst->eos = &eos[negative_component(oldhs)];
	
    FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_dens,2,&default_var);
	FT_IntrpStateVarAtCoords(front,negative_component(oldhs),p0,
		m_dens,getStateDens,&newst->dens,&default_var);
	
    FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_engy,2,&default_var);
	FT_IntrpStateVarAtCoords(front,negative_component(oldhs),p0,
		m_engy,getStateEngy,&newst->engy,&default_var);
	
    FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_kturb,2,&default_var);
	FT_IntrpStateVarAtCoords(front,negative_component(oldhs),p0,
		m_kturb,getStateKTurb,&newst->k_turb,&default_var);
	
    for (i = 0; i < dim; ++i)
	{
	    FT_NearestRectGridVarInRange(front,negative_component(oldhs),p0,
			m_mom[i],2,&default_var);
	    FT_IntrpStateVarAtCoords(front,negative_component(oldhs),p0,
		    m_mom[i],getStateMom[i],&newst->momn[i],&default_var);
	    newst->vel[i] = newst->momn[i]/newst->dens;
	}
	newst->pres = EosPressure(newst);
	set_state_max_speed(front,newst,Coords(oldp));

	newst = (STATE*)right_state(newp);
	oldst = (STATE*)sr;
	*newst = *oldst;
	newst->dim = dim;
	newst->eos = &eos[positive_component(oldhs)];
	
    //TODO: Temperature, Viscosity

    FT_NearestRectGridVarInRange(front,positive_component(oldhs),p0,
			m_dens,2,&default_var);
	FT_IntrpStateVarAtCoords(front,positive_component(oldhs),p0,
		m_dens,getStateDens,&newst->dens,&default_var);

	FT_NearestRectGridVarInRange(front,positive_component(oldhs),p0,
			m_engy,2,&default_var);
	FT_IntrpStateVarAtCoords(front,positive_component(oldhs),p0,
		m_engy,getStateEngy,&newst->engy,&default_var);
	
	FT_NearestRectGridVarInRange(front,positive_component(oldhs),p0,
			m_kturb,2,&default_var);
	FT_IntrpStateVarAtCoords(front,positive_component(oldhs),p0,
		m_kturb,getStateKTurb,&newst->k_turb,&default_var);
	
    for (i = 0; i < dim; ++i)
	{
	    FT_NearestRectGridVarInRange(front,positive_component(oldhs),p0,
			m_mom[i],2,&default_var);
	    FT_IntrpStateVarAtCoords(front,positive_component(oldhs),p0,
		    m_mom[i],getStateMom[i],&newst->momn[i],&default_var);
	    newst->vel[i] = newst->momn[i]/newst->dens;
	}
	newst->pres = EosPressure(newst);
	set_state_max_speed(front,newst,Coords(oldp));
}	/* end contact_point_propagate */

static void rgbody_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
    rgbody_point_propagate_in_fluid(front,wave,
            oldp,newp,oldhse,oldhs,dt,V);
    
    /*
    //TODO: can we get rid of rgbody_point_propagate_in_vacuum()???
	
    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    if (af_params->no_fluid == YES)
        rgbody_point_propagate_in_vacuum(front,wave,
                oldp,newp,oldhse,oldhs,dt,V);
    else
        rgbody_point_propagate_in_fluid(front,wave,
                oldp,newp,oldhse,oldhs,dt,V);
    */

}	/* end rgbody_point_propagate */

static void rgbody_point_propagate_in_fluid(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
    double vel[MAXD];
    int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	
	double *m_pres = eqn_params->pres;
	double *m_dens = eqn_params->dens;
	double *m_engy = eqn_params->engy;
    double *m_temp = eqn_params->temp;
    double *m_mu = eqn_params->mu;
    double *m_kturb = eqn_params->k_turb;
	
    double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}

    FT_NormalAtPoint(oldp,front,nor,comp);
	dn = grid_size_in_direction(nor,h,dim);

	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;

        if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
        {
            if(!debugging("collision_off"))
            {
                for (i = 0; i < dim; ++i)
                    newst->x_old[i] = Coords(oldp)[i];
            }
            
            double omega_dt,crds_com[MAXD];
            omega_dt = angular_velo(oldhs)*dt;

            //TODO: test/verify
            for (i = 0; i < dim; ++i)
    	    {
                vel[i] = center_of_mass_velo(oldhs)[i];
                crds_com[i] = Coords(oldp)[i] + dt*(vel[i] + oldst->vel[i])
				*0.5 - rotation_center(oldhs)[i];
	    }
            if (dim == 2)
            {
		vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
			angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
		vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
			angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
                for (i = 0; i < dim; ++i)
                {
                    Coords(newp)[i] = Coords(oldp)[i] + dt*(vel[i] + 
					oldst->vel[i])*0.5;
                    newst->vel[i] = vel[i];
                    FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,
						Coords(newp),front);
                }
	    }
	    else if (dim == 3)
	    {
		vel[0] += -p_angular_velo(oldhs)[2] * crds_com[1]
                          +p_angular_velo(oldhs)[1] * crds_com[2];
                vel[1] +=  p_angular_velo(oldhs)[2] * crds_com[0]
                          -p_angular_velo(oldhs)[0] * crds_com[2];
                vel[2] += -p_angular_velo(oldhs)[1] * crds_com[0]
                          +p_angular_velo(oldhs)[0] * crds_com[1];
		// propagate by euler parameters
		if (motion_type(oldhs) == ROTATION ||
		    motion_type(oldhs) == PRESET_ROTATION)
		{
                    double A[3][3],AI[3][3];
                    double ep[4];
                    int j,k;
                    double initial[MAXD];
                    for (i = 0; i< 4; i++)
                        ep[i] = old_euler_params(oldhs)[i];
                    AI[0][0] =   ep[0]*ep[0] + ep[1]*ep[1]
                               - ep[2]*ep[2] - ep[3]*ep[3];
                    AI[0][1] = 2.0 * (ep[1]*ep[2] + ep[0]*ep[3]);
                    AI[0][2] = 2.0 * (ep[1]*ep[3] - ep[0]*ep[2]);
                    AI[1][0] = 2.0 * (ep[1]*ep[2] - ep[0]*ep[3]);
                    AI[1][1] =   ep[0]*ep[0] - ep[1]*ep[1]
                               + ep[2]*ep[2] - ep[3]*ep[3];
                    AI[1][2] = 2.0 * (ep[2]*ep[3] + ep[0]*ep[1]);
                    AI[2][0] = 2.0 * (ep[1]*ep[3] + ep[0]*ep[2]);
                    AI[2][1] = 2.0 * (ep[2]*ep[3] - ep[0]*ep[1]);
                    AI[2][2] =   ep[0]*ep[0] - ep[1]*ep[1]
                               - ep[2]*ep[2] + ep[3]*ep[3];
                    for (j = 0; j < 3; j++)
                    {
                        initial[j] = 0.0;
                        for (k = 0; k < 3; k++)
                            initial[j] += AI[j][k]*crds_com[k];
                    }
                    for (i = 0; i< 4; i++)
                        ep[i] = euler_params(oldhs)[i];
                    A[0][0] =   ep[0]*ep[0] + ep[1]*ep[1]
                              - ep[2]*ep[2] - ep[3]*ep[3];
                    A[0][1] = 2.0 * (ep[1]*ep[2] - ep[0]*ep[3]);
                    A[0][2] = 2.0 * (ep[1]*ep[3] + ep[0]*ep[2]);
                    A[1][0] = 2.0 * (ep[1]*ep[2] + ep[0]*ep[3]);
                    A[1][1] =   ep[0]*ep[0] - ep[1]*ep[1]
                              + ep[2]*ep[2] - ep[3]*ep[3];
                    A[1][2] = 2.0 * (ep[2]*ep[3] - ep[0]*ep[1]);
                    A[2][0] = 2.0 * (ep[1]*ep[3] - ep[0]*ep[2]);
                    A[2][1] = 2.0 * (ep[2]*ep[3] + ep[0]*ep[1]);
                    A[2][2] =   ep[0]*ep[0] - ep[1]*ep[1]
                              - ep[2]*ep[2] + ep[3]*ep[3];
                    for (j = 0; j < 3; j++)
                    {
                        Coords(newp)[j] = rotation_center(oldhs)[j];
                        for (k = 0; k < 3; k++)
                            Coords(newp)[j] += A[j][k]*initial[k];
                    }
		}
		else
		    for (i = 0; i < dim; ++i)
                        Coords(newp)[i] = Coords(oldp)[i] + 
					dt*(vel[i] + oldst->vel[i])*0.5;
		for (i = 0; i < dim; ++i)
                {
                    newst->vel[i] = vel[i];
                    FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,
                                        Coords(newp),front);
		}
	    }
        }
        else
        {
            fourth_order_point_propagate(front,NULL,oldp,newp,oldhse,
                                    oldhs,dt,vel);
        }
	
	
	FT_IntrpStateVarAtCoords(front,comp,p1,m_kturb,
			getStateKTurb,&newst->k_turb,&oldst->k_turb);
    

    //TODO: Need to use wall model for temp?
    //
    FT_IntrpStateVarAtCoords(front,comp,p1,m_temp,
        getStateTemp,&newst->temp,&oldst->temp);
    
	FT_IntrpStateVarAtCoords(front,comp,p1,m_mu,
			getStateMu,&newst->mu,&oldst->mu);
	

	FT_IntrpStateVarAtCoords(front,comp,p1,m_pres,
			getStatePres,&newst->pres,&oldst->pres);
    
    FT_IntrpStateVarAtCoords(front,comp,p1,m_dens,
			getStateDens,&newst->dens,&oldst->dens);

    for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	    newst->momn[i] = newst->dens*vel[i];
	}

    newst->eos = oldst->eos; //TODO: Will this have valid memory?
    newst->engy = EosEnergy(newst);
    //newst->temp = EosTemperature(newst);
    //newst->mu = EosViscosity(newst);
	
    
    double s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
	set_state_max_speed(front,newst,Coords(newp));
        
    if(!debugging("collision_off"))
    {
        /* copy newst to the other STATE; used in collision solver */
        if (gas_comp(negative_component(oldhs)))
            std::copy(newst, newst+1, (STATE*)right_state(newp));
        else if (gas_comp(positive_component(oldhs)))
            std::copy(newst, newst+1, (STATE*)left_state(newp));
    }
}	/* end rgbody_point_propagate_in_fluid */

//TODO: is this necessary?
static void rgbody_point_propagate_in_vacuum(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (gas_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (gas_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}

        if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
        {
            if(!debugging("collision_off"))
            {
                for (i = 0; i < dim; ++i)
                    newst->x_old[i] = Coords(oldp)[i];
            }
            
            double omega_dt,crds_com[MAXD];
            omega_dt = angular_velo(oldhs)*dt;

            for (i = 0; i < dim; ++i)
    	    {
                vel[i] = center_of_mass_velo(oldhs)[i];
                crds_com[i] = Coords(oldp)[i] + dt*(vel[i] + oldst->vel[i])
				*0.5 - rotation_center(oldhs)[i];
	    }
            if (dim == 2)
            {
		vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
			angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
		vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
			angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
                for (i = 0; i < dim; ++i)
                {
                    Coords(newp)[i] = Coords(oldp)[i] + dt*(vel[i] + 
					oldst->vel[i])*0.5;
                    newst->vel[i] = vel[i];
                    FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,
						Coords(newp),front);
                }
	    }
	    else if (dim == 3)
	    {
		vel[0] += -p_angular_velo(oldhs)[2] * crds_com[1]
                          +p_angular_velo(oldhs)[1] * crds_com[2];
                vel[1] +=  p_angular_velo(oldhs)[2] * crds_com[0]
                          -p_angular_velo(oldhs)[0] * crds_com[2];
                vel[2] += -p_angular_velo(oldhs)[1] * crds_com[0]
                          +p_angular_velo(oldhs)[0] * crds_com[1];
		// propagate by euler parameters
		if (motion_type(oldhs) == ROTATION ||
		    motion_type(oldhs) == PRESET_ROTATION)
		{
                    double A[3][3],AI[3][3];
                    double ep[4];
                    int j,k;
                    double initial[MAXD];
                    for (i = 0; i< 4; i++)
                        ep[i] = old_euler_params(oldhs)[i];
                    AI[0][0] =   ep[0]*ep[0] + ep[1]*ep[1]
                               - ep[2]*ep[2] - ep[3]*ep[3];
                    AI[0][1] = 2.0 * (ep[1]*ep[2] + ep[0]*ep[3]);
                    AI[0][2] = 2.0 * (ep[1]*ep[3] - ep[0]*ep[2]);
                    AI[1][0] = 2.0 * (ep[1]*ep[2] - ep[0]*ep[3]);
                    AI[1][1] =   ep[0]*ep[0] - ep[1]*ep[1]
                               + ep[2]*ep[2] - ep[3]*ep[3];
                    AI[1][2] = 2.0 * (ep[2]*ep[3] + ep[0]*ep[1]);
                    AI[2][0] = 2.0 * (ep[1]*ep[3] + ep[0]*ep[2]);
                    AI[2][1] = 2.0 * (ep[2]*ep[3] - ep[0]*ep[1]);
                    AI[2][2] =   ep[0]*ep[0] - ep[1]*ep[1]
                               - ep[2]*ep[2] + ep[3]*ep[3];
                    for (j = 0; j < 3; j++)
                    {
                        initial[j] = 0.0;
                        for (k = 0; k < 3; k++)
                            initial[j] += AI[j][k]*crds_com[k];
                    }
                    for (i = 0; i< 4; i++)
                        ep[i] = euler_params(oldhs)[i];
                    A[0][0] =   ep[0]*ep[0] + ep[1]*ep[1]
                              - ep[2]*ep[2] - ep[3]*ep[3];
                    A[0][1] = 2.0 * (ep[1]*ep[2] - ep[0]*ep[3]);
                    A[0][2] = 2.0 * (ep[1]*ep[3] + ep[0]*ep[2]);
                    A[1][0] = 2.0 * (ep[1]*ep[2] + ep[0]*ep[3]);
                    A[1][1] =   ep[0]*ep[0] - ep[1]*ep[1]
                              + ep[2]*ep[2] - ep[3]*ep[3];
                    A[1][2] = 2.0 * (ep[2]*ep[3] - ep[0]*ep[1]);
                    A[2][0] = 2.0 * (ep[1]*ep[3] - ep[0]*ep[2]);
                    A[2][1] = 2.0 * (ep[2]*ep[3] + ep[0]*ep[1]);
                    A[2][2] =   ep[0]*ep[0] - ep[1]*ep[1]
                              - ep[2]*ep[2] + ep[3]*ep[3];
                    for (j = 0; j < 3; j++)
                    {
                        Coords(newp)[j] = rotation_center(oldhs)[j];
                        for (k = 0; k < 3; k++)
                            Coords(newp)[j] += A[j][k]*initial[k];
                    }
		}
		else
		    for (i = 0; i < dim; ++i)
                        Coords(newp)[i] = Coords(oldp)[i] + 
					dt*(vel[i] + oldst->vel[i])*0.5;
		for (i = 0; i < dim; ++i)
                {
                    newst->vel[i] = vel[i];
                    FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,
                                        Coords(newp),front);
		}
	    }
        }
        else
        {
            fourth_order_point_propagate(front,NULL,oldp,newp,oldhse,
                                    oldhs,dt,vel);
        }
	
        for (i = 0; i < dim; ++i) newst->vel[i] = vel[i];
	
        if(!debugging("collision_off"))
        {
            /* copy newst to the other STATE; used in collision solver */
            if (gas_comp(negative_component(oldhs)))
                std::copy(newst, newst+1, (STATE*)right_state(newp));
            else if (gas_comp(positive_component(oldhs)))
                std::copy(newst, newst+1, (STATE*)left_state(newp));
        }
        return;
}       /* end end rgbody_point_propagate_in_vacuum */

static void set_state_max_speed(
	Front *front,
	STATE *state,
	double *coords)
{
	int dim = front->rect_grid->dim;
	double s = 0.0;
	for (int i = 0; i < dim; ++i)
    {
	    s += sqr(state->momn[i]/state->dens);
    }
    s = sqrt(s);
	double c = EosSoundSpeed(state);
	set_max_front_speed(dim,s+c,NULL,coords,front);
}	/* end set_state_max_speed */

// Flux of Riemann solution of Burgers equation u_t + uu_x = 0

double burger_flux(	
	double ul,
	double um,
	double ur)
{
	double u_Rl,u_Rr;
	if (ul < um)
	{
	    if (ul > 0.0) u_Rl = ul;
	    else if (um < 0.0) u_Rl = um;
	    else u_Rl = 0.0;
	}
	else
	{
	    if (ul + um > 0.0) u_Rl = ul;
	    else u_Rl = um;
	}

	if (um < ur)
	{
	    if (um > 0.0) u_Rr = um;
	    else if (ur < 0.0) u_Rr = ur;
	    else u_Rr = 0.0;
	}
	else
	{
	    if (um + ur > 0.0) u_Rr = um;
	    else u_Rr = ur;
	}
	return 0.5*(u_Rr*u_Rr - u_Rl*u_Rl);
}	/* end flux */

// Flux of Riemann solution of linear equation u_t + au_x = 0

extern double linear_flux(	
	double a,
	double ul,
	double um,
	double ur)
{
	if (a > 0.0)
	    return a*(um - ul);
	else
	    return a*(ur - um);
}	/* end upwind_flux */

extern void reflectVectorThroughPlane(
	double *vec,
	double *nor,
	double *vec_ref,
	int dim)
{
	int i;
	double vec_nor[MAXD];
	for (i = 0; i < dim; ++i)
	{
	    vec_nor[i] = vec[i]*fabs(nor[i]);
	    vec_ref[i] = vec[i] - 2.0*vec_nor[i];
	}	
}	/* end reflectVectorThroughPlane */

extern boolean reflectNeumannState(
	Front *front,
	HYPER_SURF *hs,
	double *coords,
	COMPONENT comp,
	SWEEP *m_vst,
	STATE *state)
{
	int i,dim = front->rect_grid->dim;
	double coordsref[MAXD],nor[MAXD];
	double momn[MAXD];

	if (!FrontReflectPointViaNeumannBdry(coords,coordsref,nor,comp,
				hs,front))
	{
	    printf("ERROR: in appendGhostBuffer(), cannot reflect point!\n");
	    return NO;
	}
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->dens,getStateDens,
					&state->dens,NULL);
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->engy,getStateEngy,
					&state->engy,NULL);
	FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->pres,getStatePres,
					&state->pres,NULL);
	for (i = 0; i < dim; ++i)
	{
	    FT_IntrpStateVarAtCoords(front,comp,coordsref,m_vst->momn[i],
				getStateMom[i],&momn[i],NULL);
	}
        reflectVectorThroughPlane(momn,nor,state->momn,dim);
	return YES;
}	/* end reflectNeumannState */	

extern void findGhostState(
	STATE intfc_st,
	STATE inter_st,
	STATE *ghost_st)
{
	double vel[MAXD];
	vel[0] = inter_st.momn[0]/inter_st.dens;
	vel[1] = inter_st.momn[1]/inter_st.dens;
	vel[2] = inter_st.momn[2]/inter_st.dens;


	ghost_st->dens = intfc_st.dens;
	ghost_st->pres = intfc_st.pres;
	ghost_st->momn[0] = intfc_st.dens*vel[0];
	ghost_st->momn[1] = intfc_st.dens*vel[1];
	ghost_st->momn[2] = intfc_st.dens*vel[2];
	ghost_st->engy = EosEnergy(ghost_st);
}	/* end findGhostState */

/* 	Function computes velocity of center of mass and
 *  	angular velocity of a rigid body, must be a closed curve. 
*/
extern	void cfluid_compute_force_and_torque(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	switch (fr->rect_grid->dim)
	{
	case 2:
	    return cfluid_compute_force_and_torque2d(fr,hs,dt,force,torque);
	case 3:
	    return cfluid_compute_force_and_torque3d(fr,hs,dt,force,torque);
	}
}	/* end cfluid_compute_force_and_torque */

static	void cfluid_compute_force_and_torque2d(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	RECT_GRID *gr = computational_grid(fr->interf);
	double f[MAXD],rr[MAXD];
	double t,pres;
	double bnor[MAXD],posn[MAXD];
    double area;
	BOND *b;
	boolean pos_side;
	int i,dim = gr->dim;
	EQN_PARAMS *cFparams = (EQN_PARAMS*)fr->extra1;
	double *gravity = cFparams->gravity;
	CURVE *curve = Curve_of_hs(hs);

	if (debugging("rigid_body"))
	    (void) printf("Entering cfluid_compute_force_and_torque2d()\n");

	if (gas_comp(negative_component(curve)))
	    pos_side = NO;
	else 
	    pos_side = YES;

	for (i = 0; i < dim; ++i)
	{
	    force[i] = 0.0;
	}
	*torque = 0.0;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    if (force_on_hse(Hyper_surf_element(b),Hyper_surf(curve),gr,
			&pres,bnor,posn,pos_side))
	    {
            area = bond_length(b);
            double mag_bnor = Mag2d(bnor);
	    	for (i = 0; i < dim; ++i)
	    	{
		        f[i] = pres*area*bnor[i]/mag_bnor;
	    	    rr[i] = posn[i] - rotation_center(curve)[i];
	    	    force[i] += f[i];
	    	}
	    	Cross2d(rr,f,t);
	    	*torque += t;
	    }
	}
	 /* Add gravity to the total force */
	if (motion_type(curve) != ROTATION)
	{
	    for (i = 0; i < dim; ++i)
	    	force[i] += gravity[i]*total_mass(curve);
	}
	if (debugging("rigid_body"))
	{
	    (void) printf("Leaving cfluid_compute_force_and_torque2d()\n");
	    (void) printf("total_force = %f %f\n",force[0],force[1]);
	    (void) printf("torque = %f\n",*torque);
	}
}	/* end cfluid_compute_force_and_torque2d */

#define         MAX_TRI_FOR_INTEGRAL            100
static	void cfluid_compute_force_and_torque3d(
	Front *fr,
	HYPER_SURF *hs,
	double dt,
	double *force,
	double *torque)
{
	RECT_GRID *gr = computational_grid(fr->interf);
	double f[MAXD],rr[MAXD];
	double t[MAXD],tdir,pres;
	double tnor[MAXD],posn[MAXD];
    double area;
	TRI *tri;
	boolean pos_side;
	int i,dim = gr->dim;
	EQN_PARAMS *cFparams = (EQN_PARAMS*)fr->extra1;
	double *gravity = cFparams->gravity;
	SURFACE *surface = Surface_of_hs(hs);

	if (gas_comp(negative_component(surface)))
	    pos_side = NO;
	else 
	    pos_side = YES;

	for (i = 0; i < dim; ++i)
	{
	    force[i] = 0.0;
	    torque[i] = 0.0;
	}
	for (tri = first_tri(surface); !at_end_of_tri_list(tri,surface); 
			tri = tri->next)
	{
	    if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
			&pres,tnor,posn,pos_side))
	    {
            area = tri_area(tri);
            double mag_tnor = Mag3d(tnor);
	    	for (i = 0; i < dim; ++i)
	    	{
		        f[i] = pres*area*tnor[i]/mag_tnor;
	    	    force[i] += f[i];
		        rr[i] = posn[i] - rotation_center(surface)[i];
		    }
		
            Cross3d(rr,f,t);
		    tdir = Dot3d(t,(rotation_direction(hs)));
	    	for (i = 0; i < dim; ++i)
		    {
		        t[i] = tdir*rotation_direction(hs)[i];
		        torque[i] += t[i];
		    }
	    }
	}
	 /* Add gravity to the total force */
	if (motion_type(surface) != ROTATION)
	{
	    for (i = 0; i < dim; ++i)
	    	force[i] += gravity[i]*total_mass(surface);
	}
	if (debugging("rigid_body"))
	{
	    printf("In cfluid_compute_force_and_torque3d()\n");
	    printf("total_force = %f %f %f\n",force[0],force[1],force[2]);
	    printf("torque = %f %f %f\n",torque[0],torque[1],torque[2]);
	}
}	/* end cfluid_compute_force_and_torque3d */


static boolean force_on_hse(
	HYPER_SURF_ELEMENT *hse,	/* Bond (2D) or tri (3D) */
	HYPER_SURF *hs,			/* Curve (2D) or surface (3D) */
	RECT_GRID *gr,			/* Rectangular grid */
	double *pres,		/* Average pressure */
	double *nor,		/* normal vector, pointing onto body */
	double *posn,		/* Position of the pressure */
	boolean pos_side)	/* Is the body on the positive side of hs? */
{
	int dim = gr->dim;
	switch (dim)
	{
	case 2: 
	    return force_on_hse2d(hse,hs,gr,pres,nor,posn,pos_side);
	case 3: 
	    return force_on_hse3d(hse,hs,gr,pres,nor,posn,pos_side);
	default: 
	    return NO; 
	}
	
}	/* end force_on_hse */

static boolean force_on_hse2d(
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	RECT_GRID *gr,
	double *pres,
	double *nor,
	double *posn,
	boolean pos_side)
{
	double crds1[MAXD],crds2[MAXD];
	double p1,p2;
	Locstate s1,s2;
	BOND *b = Bond_of_hse(hse);
	CURVE *c = Curve_of_hs(hs);
	double *L = gr->L;
	double *U = gr->U;
	int i;
	
	/* Get pressure at two end points of the bond */
	if (b->start == c->start->posn)
	    s1 = pos_side ? right_start_state(c) : left_start_state(c);
	else
	    s1 = pos_side ? right_state(b->start) : left_state(b->start);
	if (b->end == c->end->posn)
	    s2 = pos_side ? right_end_state(c) : left_end_state(c);
	else
	    s2 = pos_side ? right_state(b->end) : left_state(b->end);

	p1 = getStatePres(s1);	p2 = getStatePres(s2);
	for (i = 0; i < 2; ++i)
	{
	    crds1[i] = Coords(b->start)[i];
	    crds2[i] = Coords(b->end)[i];
	}

	/* Cut and interpolate if one end is outside the domain */
	for (i = 0; i < 2; ++i)
	{
	    if (crds1[i] <= L[i])
	    {
		if (crds2[i] <= L[i]) return NO; // both ends out
		else
		{
		    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p1 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
		    crds1[i] = L[i];
		}
	    }
	    if (crds1[i] >= U[i])
	    {
		if (crds2[i] >= U[i]) return NO; // both ends out
		else
		{
		    crds1[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p1 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
		    crds1[i] = U[i];
		}
	    }
	}
	for (i = 0; i < 2; ++i)
	{
	    if (crds2[i] <= L[i])
	    {
		if (crds1[i] <= L[i]) return NO; // both ends out
		else
		{
		    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],L[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p2 = intrp_between(crds1[i],crds2[i],L[i],p1,p2);
		    crds2[i] = L[i];
		}
	    }
	    if (crds2[i] >= U[i])
	    {
		if (crds1[i] >= U[i]) return NO; // both ends out
		else
		{
		    crds2[(i+1)%2] = intrp_between(crds1[i],crds2[i],U[i],
				crds1[(i+1)%2],crds2[(i+1)%2]);
		    p2 = intrp_between(crds1[i],crds2[i],U[i],p1,p2);
		    crds2[i] = U[i];
		}
	    }
	}
	nor[0] = pos_side ? crds1[1] - crds2[1] : crds2[1] - crds1[1];
	nor[1] = pos_side ? crds2[0] - crds1[0] : crds1[0] - crds2[0];
	*pres = 0.5*(p1 + p2);
	posn[0] = 0.5*(crds1[0] + crds2[0]);
	posn[1] = 0.5*(crds1[1] + crds2[1]);
	return YES;
}	/* end force_on_hse2d */

static boolean force_on_hse3d(
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	RECT_GRID *gr,
	double *pres,
	double *nor,
	double *posn,
	boolean pos_side)
{
        TRI *t = Tri_of_hse(hse);
	POINT *point;
	Locstate sl,sr;
        int i,j,dim = gr->dim;

	*pres = 0.0;
	for (i = 0; i < 3; ++i)
	    posn[i] = 0.0;
	for (i = 0; i < 3; ++i)
	{
	    point = Point_of_tri(t)[i];
	    for (j = 0; j < dim; ++j)
		posn[j] += Coords(point)[j];
	    FT_GetStatesAtPoint(point,hse,hs,&sl,&sr);
	    if (pos_side)
		*pres += getStatePres(sr);
	    else
		*pres += getStatePres(sl);
	}
	*pres /= 3.0;
	for (i = 0; i < dim; ++i)
	{
	    nor[i] = pos_side ? -Tri_normal(t)[i] : Tri_normal(t)[i];
	    posn[i] /= 3.0;
	}
	
    //TODO: Investigate this comment
    //
    /* Need to treat subdomain boundary */
	
    return YES;
}	/* end force_on_hse3d */

static double intrp_between(
	double x1,
	double x2,
	double x,
	double y1,
	double y2)
{
	double y;
	if (x1 == x2) return y1;
	y = y1 + (y2 - y1)/(x2 - x1)*(x - x1);
	return y;
}

static void get_variable_bdry_params(
	int dim,
	FILE *infile,
	POINTER *func_params)
{
	static VAR_BDRY_PARAMS params;
	int nbr_pistons;
	double half_angular_width;
	int i;

	params.dim=dim;
	CursorAfterString(infile,"Enter the center of the circle:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&params.center[i]);
            (void) printf("%f ",params.center[i]);
	}
	(void) printf("\n");

	CursorAfterString(infile,"Enter number of pistons:");
	fscanf(infile,"%d",&nbr_pistons);
	(void) printf("%d\n",nbr_pistons);

	params.number_pistons = nbr_pistons;

	FT_VectorMemoryAlloc((POINTER*)&params.angles_pistons,nbr_pistons+1,
				sizeof(double));
	for (i = 0; i < nbr_pistons+1; ++i)
		params.angles_pistons[i] = 2*PI*i/nbr_pistons;
	params.number_pistons = nbr_pistons + 1;

	CursorAfterString(infile,"Enter angular width of pistons:");
	fscanf(infile,"%lf",&half_angular_width);
	(void) printf("%f\n",half_angular_width);
	params.half_angular_width = half_angular_width;

	CursorAfterString(infile,
		"Enter radial velocity, density and pressure at piston:");
	fscanf(infile,"%lf %lf %lf",&params.bdry_vel,&params.bdry_dens,
					&params.bdry_pres);
	(void) printf("%f %f %f\n",params.bdry_vel,params.bdry_dens,
					params.bdry_pres);
	CursorAfterString(infile,
		"Enter time duration of the piston:");
	fscanf(infile,"%lf",&params.jet_duration_time);
	(void) printf("%f\n",params.jet_duration_time);

	*func_params = (POINTER)&params;
}	/* end get_variable_bdry_params */

extern void restart_set_dirichlet_bdry_function(Front *front)
{
    INTERFACE *intfc = front->interf;
    int i;
    BOUNDARY_STATE  *bstate;
    const char *s;
	CURVE **curves;
	SURFACE **surfs;
	int comp;
	STATE *state;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	int dim = Dimension(intfc);

    for (i = 0; i < num_bstates(intfc); ++i)
    {
        bstate = bstate_list(intfc)[i];
        if (bstate == NULL) continue;
        
        s = bstate->_boundary_state_function_name;
        if (s == NULL) continue;
        
        if (strcmp(s,"cF_flowThroughBoundaryState") == 0)
            bstate->_boundary_state_function = cF_flowThroughBoundaryState;
        else if (strcmp(s,"cF_constantWithWhiteNoise") == 0)
            bstate->_boundary_state_function = cF_constantWithWhiteNoise;
        else if (strcmp(s,"cF_supersonicOutflowState") == 0)
            bstate->_boundary_state_function = cF_supersonicOutflowState;
        else if (strcmp(s,"cF_farfieldBoundaryState") == 0)
            bstate->_boundary_state_function = cF_farfieldBoundaryState;
        else if (strcmp(s,"cF_variableBoundaryState") == 0)
            bstate->_boundary_state_function = cF_variableBoundaryState;
        else
        {
            printf("\n\nUnrecognized Dirichlet boundary type!\n\n");
            LOC(); clean_up(EXIT_FAILURE);
        }
    }

	switch (dim)
	{
	case 2:
	    intfc_curve_loop(intfc,curves)
	    {
	    	if (wave_type(*curves) != DIRICHLET_BOUNDARY) continue;
	    	state = (STATE*)boundary_state(Hyper_surf(*curves));
	    	if (state == NULL) continue;
	    	comp = gas_comp(positive_component(*curves)) ?
                                positive_component(*curves) :
                                negative_component(*curves);
	    	state->eos = &(eqn_params->eos[comp]);
	    }
	    break;
	case 3:
	    intfc_surface_loop(intfc,surfs)
	    {
	    	if (wave_type(*surfs) != DIRICHLET_BOUNDARY) continue;
	    	state = (STATE*)boundary_state(Hyper_surf(*surfs));
	    	if (state == NULL) continue;
	    	comp = gas_comp(positive_component(*surfs)) ?
                                positive_component(*surfs) :
                                negative_component(*surfs);
	    	state->eos = &(eqn_params->eos[comp]);
	    }
	    break;
	}
}       /* end restart_set_dirichlet_bdry_function */

extern void read_open_end_bdry_data(
	char *inname,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	int i,j,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	char msg[256],string[256];
	static OPEN_PIPE_PARAMS pipe_params;

	if (dim != 3) return;
	for (i = 0; i < dim; ++i)
	for (j = 0; j < 2; ++j)
	{
	    if (rect_boundary_type(intfc,i,j) == OPEN_BOUNDARY)
	    {
		printf("Available open boundary types are\n");
		printf("\tOpen pipe (p)\n");
		sprintf(msg,"For open boundary in dir %d side %d",i,j);
		CursorAfterString(infile,msg);
		CursorAfterString(infile,"Enter boundary function type: ");
		fscanf(infile,"%s",string);
		(void) printf("%s\n",string);
		switch (string[0])
		{
		case 'p':
		case 'P':
		    pipe_params.dir = i;
		    pipe_params.side = j;
		    CursorAfterString(infile,"Enter center of pile: ");
		    fscanf(infile,"%lf %lf %lf",&pipe_params.center[0],
						&pipe_params.center[1],
						&pipe_params.center[2]);
		    (void) printf("%f %f %f\n",pipe_params.center[0],
						pipe_params.center[1],
						pipe_params.center[2]);
		    CursorAfterString(infile,"Enter mean radius of pile: ");
		    fscanf(infile,"%lf",&pipe_params.radius);
		    (void) printf("%f\n",pipe_params.radius);
		    CursorAfterString(infile,"Enter inside boundary type: ");
		    fscanf(infile,"%s",string);
		    (void) printf("%s\n",string);
		    switch (string[0])
		    {
		    case 'D':
			pipe_params.in_pipe_bdry = DIRICHLET_BOUNDARY;
		    	CursorAfterString(infile,
				"Constant state (c) or flow through (f)? ");
		    	fscanf(infile,"%s",string);
		    	(void) printf("%s\n",string);
			if (string[0] == 'c' || string[0] == 'C')
			{
			    pipe_params.in_flow_through = NO;
			    CursorAfterString(infile,"Enter density: ");
		    	    fscanf(infile,"%lf",
					&pipe_params.state[0].dens);
			    (void) printf("%f\n",pipe_params.state[0].dens);
			    CursorAfterString(infile,"Enter pressure: ");
		    	    fscanf(infile,"%lf",
					&pipe_params.state[0].pres);
			    (void) printf("%f\n",pipe_params.state[0].pres);
			    CursorAfterString(infile,"Enter velocity: ");
		    	    fscanf(infile,"%lf %lf %lf",
					&pipe_params.state[0].vel[0],
					&pipe_params.state[0].vel[1],
					&pipe_params.state[0].vel[2]);
			    (void) printf("%f %f %f\n",
					pipe_params.state[0].vel[0],
					pipe_params.state[0].vel[1],
					pipe_params.state[0].vel[2]);
			}
			else if (string[0] == 'f' || string[0] == 'F')
			{
			    pipe_params.in_flow_through = YES;
			}
			break;
		    case 'N':
			pipe_params.in_pipe_bdry = NEUMANN_BOUNDARY;
			break;
		    default:
			(void) printf("Unknown boundary type\n");
			clean_up(ERROR);
		    }
		    CursorAfterString(infile,"Enter outside boundary type: ");
		    fscanf(infile,"%s",string);
		    (void) printf("%s\n",string);
		    switch (string[0])
		    {
		    case 'D':
			pipe_params.out_pipe_bdry = DIRICHLET_BOUNDARY;
		    	CursorAfterString(infile,
				"Constant state (c) or flow through (f)? ");
		    	fscanf(infile,"%s",string);
		    	(void) printf("%s\n",string);
			if (string[0] == 'c' || string[0] == 'C')
			{
			    pipe_params.out_flow_through = NO;
			    CursorAfterString(infile,"Enter density: ");
		    	    fscanf(infile,"%lf",
					&pipe_params.state[1].dens);
			    (void) printf("%f\n",pipe_params.state[1].dens);
			    CursorAfterString(infile,"Enter pressure: ");
		    	    fscanf(infile,"%lf",
					&pipe_params.state[1].pres);
			    (void) printf("%f\n",pipe_params.state[1].pres);
			    CursorAfterString(infile,"Enter velocity: ");
		    	    fscanf(infile,"%lf %lf %lf",
					&pipe_params.state[1].vel[0],
					&pipe_params.state[1].vel[1],
					&pipe_params.state[1].vel[2]);
			    (void) printf("%f %f %f\n",
					pipe_params.state[1].vel[0],
					pipe_params.state[1].vel[1],
					pipe_params.state[1].vel[2]);
			}
			else if (string[0] == 'f' || string[0] == 'F')
			{
			    pipe_params.out_flow_through = YES;
			}
			break;
		    case 'N':
			pipe_params.out_pipe_bdry = NEUMANN_BOUNDARY;
			break;
		    default:
			(void) printf("Unknown boundary type\n");
			clean_up(ERROR);
		    }
		    front->open_end_params = (POINTER)&pipe_params;
		    front->open_end_func = pipe_end_func;
		    break;
		default:
		    (void) printf("Unknown open boundary function type\n");
		    clean_up(ERROR);
		}
	    }
	}
}	/* end read_open_end_bdry_data */

static void pipe_end_func(
	Front* front,
	POINTER func_params,
	int *ic,
	COMPONENT comp,
	int idir,
	int side,
	int *bdry_type,
	Locstate state)
{
	OPEN_PIPE_PARAMS* d_params = (OPEN_PIPE_PARAMS*)func_params;
	EQN_PARAMS* eqn_params = (EQN_PARAMS*) front->extra1;
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	double *top_L = top_grid->L;
	double *top_U = top_grid->U;
	double *top_h = top_grid->h;
	int *lbuf = front->rect_grid->lbuf;
	int *ubuf = front->rect_grid->ubuf;
	int *top_gmax = top_grid->gmax;
	double *c = d_params->center;
	double r = d_params->radius;
	double coords[MAXD],dist = 0.0;
	int i, ic_nb[MAXD], index_nb;

	((STATE*)state)->dim = 3;
	((STATE*)state)->eos = &(eqn_params->eos[comp]);

	if (!gas_comp(comp))
	{
	    ((STATE*)state)->dens = 0.0;
	    for (i = 0; i < 3; ++i)
		((STATE*)state)->vel[i] = ((STATE*)state)->momn[i] = 0.0;
	    ((STATE*)state)->pres = 0.0;
	    ((STATE*)state)->engy = 0.0;
	    *bdry_type = PASSIVE_BOUNDARY;
	    return;
	}

	for (i = 0; i < 3; ++i)
	{
	    ic_nb[i] = ic[i];
	    coords[i] = top_L[i] + ic[i] * top_h[i];
	    dist += sqr(coords[i] - c[i]);
	}
	dist = sqrt(dist);
	if (dist < r)
	{
	    *bdry_type = d_params->in_pipe_bdry;
	    switch (*bdry_type)
	    {
	    case DIRICHLET_BOUNDARY:
		if (d_params->in_flow_through)
		{
		    switch (side)
		    {
		    case 0:
			ic_nb[idir] = lbuf[idir];
			break;
		    case 1:
			ic_nb[idir] = top_gmax[idir] - ubuf[idir];
		    default:
			(void) printf("Undefined side \n");
			clean_up(ERROR);
		    }
		    index_nb = d_index3d(ic_nb[0],ic_nb[1],ic_nb[2],top_gmax);
		    ((STATE*)state)->dens = eqn_params->dens[index_nb];
		    for (i = 0; i < 3; ++i)
		    {
			((STATE*)state)->vel[i] = eqn_params->vel[i][index_nb];
			((STATE*)state)->momn[i] = eqn_params->mom[i][index_nb];
		    }
		    ((STATE*)state)->pres = eqn_params->pres[index_nb];
		    ((STATE*)state)->engy = eqn_params->engy[index_nb];
		}
		else
		{
		    ((STATE*)state)->dens = d_params->state[0].dens;
		    ((STATE*)state)->pres = d_params->state[0].pres;
		    for (i = 0; i < 3; ++i)
		    {
			((STATE*)state)->vel[i] = d_params->state[0].vel[i];
			((STATE*)state)->momn[i] = d_params->state[0].dens * 
					d_params->state[0].vel[i];
		    }
		    ((STATE*)state)->engy = EosEnergy((STATE*)state);
		}
		break;
	    case NEUMANN_BOUNDARY:
		 switch (side)
                 {
                 case 0:
                     ic_nb[idir] = lbuf[idir] * 2 - 1 - ic[idir];
                     break;
                 case 1:
                     ic_nb[idir] = (top_gmax[idir] - ubuf[idir]) * 2 
						+ 1 - ic[idir];
                 default:
                     (void) printf("Undefined side \n");
                     clean_up(ERROR);
                 }
		index_nb = d_index3d(ic_nb[0],ic_nb[1],ic_nb[2],top_gmax);
		((STATE*)state)->dens = eqn_params->dens[index_nb];
		for (i = 0; i < 3; ++i)
		{
		    ((STATE*)state)->vel[i] = eqn_params->vel[i][index_nb];
		    ((STATE*)state)->momn[i] = eqn_params->mom[i][index_nb];
		    if (i == idir)
		    {
		        ((STATE*)state)->vel[i] *= -1.0;
			((STATE*)state)->momn[i] *= -1.0;
		    }
		}
		((STATE*)state)->engy = EosEnergy((STATE*)state);
		break;
	    default:
		(void) printf("Unknown inner pipe boundary type\n");
		clean_up(ERROR);
	    }
	}
	else
	{
	    *bdry_type = d_params->out_pipe_bdry;
	    switch (*bdry_type)
	    {
	    case DIRICHLET_BOUNDARY:
		if (d_params->out_flow_through)
		{
		    switch (side)
		    {
		    case 0:
			ic_nb[idir] = lbuf[idir];
			break;
		    case 1:
			ic_nb[idir] = top_gmax[idir] - ubuf[idir];
		    default:
			(void) printf("Undefined side \n");
			clean_up(ERROR);
		    }
		    index_nb = d_index3d(ic_nb[0],ic_nb[1],ic_nb[2],top_gmax);
		    ((STATE*)state)->dens = eqn_params->dens[index_nb];
		    for (i = 0; i < 3; ++i)
		    {
			((STATE*)state)->vel[i] = eqn_params->vel[i][index_nb];
			((STATE*)state)->momn[i] = eqn_params->mom[i][index_nb];
		    }
		    ((STATE*)state)->pres = eqn_params->pres[index_nb];
		    ((STATE*)state)->engy = eqn_params->engy[index_nb];
		}
		else
		{
		    ((STATE*)state)->dens = d_params->state[1].dens;
		    ((STATE*)state)->pres = d_params->state[1].pres;
		    for (i = 0; i < 3; ++i)
		    {
			((STATE*)state)->vel[i] = d_params->state[1].vel[i];
			((STATE*)state)->momn[i] = d_params->state[1].dens * 
					d_params->state[1].vel[i];
		    }
		    ((STATE*)state)->engy = EosEnergy((STATE*)state);
		}
		break;
	    case NEUMANN_BOUNDARY:
		 switch (side)
                 {
                 case 0:
                     ic_nb[idir] = lbuf[idir] * 2 - 1 - ic[idir];
                     break;
                 case 1:
                     ic_nb[idir] = (top_gmax[idir] - ubuf[idir]) * 2 
						+ 1 - ic[idir];
                 default:
                     (void) printf("Undefined side \n");
                     clean_up(ERROR);
                 }
		index_nb = d_index3d(ic_nb[0],ic_nb[1],ic_nb[2],top_gmax);
		((STATE*)state)->dens = eqn_params->dens[index_nb];
		for (i = 0; i < 3; ++i)
		{
		    ((STATE*)state)->vel[i] = eqn_params->vel[i][index_nb];
		    ((STATE*)state)->momn[i] = eqn_params->mom[i][index_nb];
		    if (i == idir)
		    {
		        ((STATE*)state)->vel[i] *= -1.0;
			((STATE*)state)->momn[i] *= -1.0;
		    }
		}
		((STATE*)state)->engy = EosEnergy((STATE*)state);
		break;
	    default:
		(void) printf("Unknown outer pipe boundary type\n");
		clean_up(ERROR);
	    }
	}

	return;
}
