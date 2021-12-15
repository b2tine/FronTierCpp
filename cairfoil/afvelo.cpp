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

#include "airfoil.h"

struct VERTICAL_PARAMS
{
	double cen[MAXD];
	double v0;
	double stop_time;
};

struct RANDOMV_PARAMS
{
	double v0[MAXD];
	double stop_time;
};

struct FIXAREA_PARAMS
{
	int num_pts;
	int *global_ids;
	double vel[MAXD];
    double stop_time;
};

struct SHAPE_PARAMS 
{
	int shape_id;
	double L[2];
	double U[2];
	double cen[2];
	double R[2];
};

static boolean within_shape(SHAPE_PARAMS,double*);

static void initVelocityFunc(FILE*,Front*);
static int zero_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);
static int toroidal_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);
static int parabolic_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);
static int singular_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);
static int vertical_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);
static int random_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);
static int marker_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);

static void init_fixarea_params(Front*,FILE*,FIXAREA_PARAMS*);
static void restart_fixarea_params(Front*,FILE*,FIXAREA_PARAMS*);
static void init_fixpoint_params(Front*,FILE*,FIXAREA_PARAMS*);

static void convert_to_point_mass(Front*, AF_PARAMS*);
static int countSurfPoints(INTERFACE* intfc);
static int countStringPoints(INTERFACE* intfc, boolean is_parachute_system);

static void checkSetGoreNodes(INTERFACE*);
static void set_gore_node(NODE*);

void setMotionParams(Front* front)
{
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	
    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	af_params->no_fluid = NO;
    if (CursorAfterStringOpt(infile,
            "Entering yes to turn off fluid solver: "))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
            af_params->no_fluid = YES;
    }
    fclose(infile);
    
    setFabricPropagators(front);
    setFabricParams(front);
	
    front->_scatter_front_extra = scatterAirfoilExtra;
}

void setFabricPropagators(Front* front)
{
    int dim = FT_Dimension();
	FILE *infile = fopen(InName(front),"r");
	char string[100];

	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    af_params->spring_model = MODEL2;//TODO: Remove -- there is no MODEL#
	
    front->vfunc = nullptr;
    front->vparams = nullptr;
	
    if (af_params->no_fluid == YES)
	{
	    front->curve_propagate = airfoil_curve_propagate;
	    front->node_propagate = airfoil_node_propagate;
	    initVelocityFunc(infile,front);
	}
	else
	{
	    front->_point_propagate = airfoil_point_propagate;
        front->curve_propagate = airfoil_curve_propagate;
        front->node_propagate = airfoil_node_propagate;
	}

	if (af_params->no_fluid == YES || af_params->is_parachute_system == NO)
	{
        CursorAfterString(infile,"Enter interior propagator:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (dim == 2)
	    {
	    	switch (string[0])
	    	{
	    	case 'n':
	    	case 'N':
	    	    front->tan_curve_propagate = nullptr;
	    	    break;
	    	case 'f':
	    	case 'F':
	    	    front->tan_curve_propagate = fixed_length_tan_curve_propagate;
	    	    break;
	    	case 'e':
	    	case 'E':
	    	    if (string[1] == '2')
	    	    	front->tan_curve_propagate = second_order_elastic_curve_propagate;
	    	    else
		        {
	    	    	front->tan_curve_propagate = fourth_order_elastic_curve_propagate;
#if defined(__GPU__)
                    if (CursorAfterStringOpt(infile,"Enter yes to use GPU solver:"))
			        {
                        fscanf(infile,"%s",string);
                        (void) printf("%s\n",string);
                        if (string[0] == 'y' || string[0] == 'Y')
                            af_params->use_gpu = YES;
                    }
#endif
		        }
	    	    break;
	    	default:
                (void) printf("Unknown interior propagator!\n");
                LOC(); clean_up(ERROR);
	    	}
	    }
	    else if (dim == 3)
	    {
	    	switch (string[0])
	    	{
	    	case 'n':
	    	case 'N':
	    	    front->interior_propagate = nullptr;
	    	    break;
	    	case 'e':
	    	case 'E':
	    	    if (string[1] == '2')
	    	    	front->interior_propagate = second_order_elastic_surf_propagate;
	    	    else
                {
                    front->interior_propagate = fourth_order_elastic_surf_propagate;
#if defined(__GPU__)
            		if (CursorAfterStringOpt(infile,"Enter yes to use GPU solver:"))
                    {
                        fscanf(infile,"%s",string);
                        (void) printf("%s\n",string);
                        if (string[0] == 'y' || string[0] == 'Y')
                            af_params->use_gpu = YES;
                    }
#endif
		        }
	    	    break;
	    	case 'p':
	    	case 'P':
	    	    front->interior_propagate = fourth_order_elastic_set_propagate;
                if (CursorAfterStringOpt(infile,"Enter yes for collision substepping:"))
                {
                    fscanf(infile,"%s",string);
                    (void) printf("%s\n",string);
                    if (string[0] == 'y' || string[0] == 'Y')
                    {
                        front->interior_propagate = elastic_set_propagate;
                        
                        CursorAfterString(infile,"Enter nsub per collision step:");
                        fscanf(infile,"%d",&af_params->collsn_step_max_nsub);
                        (void) printf("%d\n",af_params->collsn_step_max_nsub);
                    }
                }
	    	    break;
	    	default:
		    (void) printf("Unknown interior propagator!\n");
		    LOC(); clean_up(ERROR);
	    	}
	    }
	}
	else
    {
        front->interior_propagate = fourth_order_elastic_set_propagate;
        if (CursorAfterStringOpt(infile,"Enter yes for collision substepping:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
            {
                front->interior_propagate = elastic_set_propagate;

                CursorAfterString(infile,"Enter nsub per collision step:");
                fscanf(infile,"%d",&af_params->collsn_step_max_nsub);
                (void) printf("%d\n",af_params->collsn_step_max_nsub);
            }
        }
    }

	af_params->n_sub = 1;
	if (CursorAfterStringOpt(infile,"Enter interior sub step number:"))
    {
        fscanf(infile,"%d",&af_params->n_sub);
        (void) printf("%d\n",af_params->n_sub);
    }

	af_params->ss_dt_relax = 1.0;
	if (CursorAfterStringOpt(infile,"Enter spring solver time step relaxation parameter:"))
    {
        fscanf(infile,"%lf",&af_params->ss_dt_relax);
        (void) printf("%f\n",af_params->ss_dt_relax);
    }
    
    fclose(infile);
}

void setFabricParams(Front* front)
{
	int dim = FT_Dimension();
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	
    EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

    af_params->attach_gores = NO;
    if (dim == 3 && numOfGoreHsbdry(front->interf) != 0)
	{
	    af_params->attach_gores = YES;
	    checkSetGoreNodes(front->interf);
	}

    af_params->use_gpu = NO;
#if defined(__GPU__)
    if (CursorAfterStringOpt(infile,"Enter yes to use GPU solver:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		af_params->use_gpu = YES;
	}
#endif

    if (af_params->no_fluid == NO)
	{
        if (FT_FrontContainWaveType(front,ELASTIC_BOUNDARY))
        {
            eqn_params->with_porosity = NO;
            af_params->with_porosity = NO;
            if(CursorAfterStringOpt(infile,"Enter yes to use porosity:"))
            {
                fscanf(infile,"%s",string);
                (void) printf("%s\n",string);
                if (string[0] == 'y' || string[0] == 'Y')
                {
                    eqn_params->with_porosity = YES;
                    af_params->with_porosity = YES;
                }
            }

            if (eqn_params->with_porosity == YES)
            {
                if (CursorAfterStringOpt(infile,"Enter fabric porosity:"))
                {
                    fscanf(infile,"%lf",&af_params->porosity);
                    (void) printf("%f\n",af_params->porosity);
                }
                eqn_params->porosity = af_params->porosity;

                af_params->poro_scheme = PORO_SCHEME::NORMAL_REFLECTION;
                if (CursorAfterStringOpt(infile,"Enter porosity ghost fluid method:"))
                {
                    fscanf(infile,"%s",string);
                    (void) printf("%s\n",string);
                    if (string[0] == 'n' || string[0] == 'N')
                    {
                        af_params->poro_scheme = PORO_SCHEME::NORMAL_REFLECTION;
                    }
                    else if (string[0] == 'r' || string[0] == 'R')
                    {
                        if (string[1] == 'e' || string[1] == 'E')
                            af_params->poro_scheme = PORO_SCHEME::REFLECTION;
                        else if (string[1] == 'i' || string[1] == 'I')
                            af_params->poro_scheme = PORO_SCHEME::RIEMANN;
                    }
                    else if (string[0] == 'e' || string[0] == 'E')
                    {
                        af_params->poro_scheme = PORO_SCHEME::ERGUN;
                        CursorAfterString(infile,"Enter viscous parameter:");
                        fscanf(infile,"%lf",&af_params->porous_coeff[0]);
                        (void) printf("%f\n",af_params->porous_coeff[0]);
                        CursorAfterString(infile,"Enter inertial parameter:");
                        fscanf(infile,"%lf",&af_params->porous_coeff[1]);
                        (void) printf("%f\n",af_params->porous_coeff[1]);
                        eqn_params->porosity = af_params->porosity;
                        eqn_params->porous_coeff[0] = af_params->porous_coeff[0];
                        eqn_params->porous_coeff[1] = af_params->porous_coeff[1];
                    }
                
                    eqn_params->poro_scheme = af_params->poro_scheme;
                }
            }
            
            CursorAfterString(infile,"Enter area density of canopy:");
            fscanf(infile,"%lf",&af_params->area_dens);
            (void) printf("%f\n",af_params->area_dens);

        }

        af_params->fsi_startstep = 5;
        if (CursorAfterStringOpt(infile,"Enter timestep to activate FSI:"))
        {
            fscanf(infile,"%d",&af_params->fsi_startstep);
            printf("%d\n",af_params->fsi_startstep);
        }
        eqn_params->fsi_startstep = af_params->fsi_startstep;

        if (CursorAfterStringOpt(infile,"Enter yes to disable FSI:"))
        {
            fscanf(infile,"%s",string);
            printf("%s\n",string);

            if (string[0] == 'y' || string[0] == 'Y')
            {
                eqn_params->with_porosity = NO;
                eqn_params->fsi_startstep = 100000000;
                af_params->with_porosity = NO;
                af_params->fsi_startstep = 100000000;
            }
        }

        for (int i = 0; i < dim; ++i)
            af_params->gravity[i] = eqn_params->gravity[i];
	}


    af_params->inflation_assist = false;
    if (CursorAfterStringOpt(infile,"Enter yes to enable inflation assist:"))
    {
        //TODO: should we disable FSI when using this?
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
            af_params->inflation_assist = true;

        if (af_params->inflation_assist)
        {
            CursorAfterString(infile,"Enter uniform normal pressure differential:");
            fscanf(infile,"%lf",&af_params->delta_pres);
            (void) printf("%f\n",af_params->delta_pres);
        }

        //TODO: Need criteria to end inflation assist and switch
        //      back to normal fluid structure interaction routine.
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

            CursorAfterString(infile,"Enter fabric bending stiffness constant:");
            fscanf(infile,"%lf",&af_params->kbs);
            (void) printf("%f\n",af_params->kbs);

            CursorAfterString(infile,"Enter fabric bending damping constant:");
            fscanf(infile,"%lf",&af_params->lambda_bs);
            (void) printf("%f\n",af_params->lambda_bs);

            CursorAfterString(infile,"Enter fabric friction constant:");
            fscanf(infile,"%lf",&af_params->mu_s);
            (void) printf("%f\n",af_params->mu_s);

            if (CursorAfterStringOpt(infile,"Enter fabric thickness:"))
            {
                fscanf(infile,"%lf",&af_params->fabric_thickness);
                (void) printf("%f\n",af_params->fabric_thickness);
            }

            if (CursorAfterStringOpt(infile,"Enter fabric rounding tolerance:"))
            {
                fscanf(infile,"%lf",&af_params->fabric_eps);
                (void) printf("%f\n",af_params->fabric_eps);
            }

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

    af_params->m_l = 0.0;
	if ( (dim == 2 && FT_FrontContainWaveType(front,ELASTIC_STRING)) || 
	     (dim == 3 && FT_FrontContainHsbdryType(front,STRING_HSBDRY)) )
	{
        af_params->strings_present = true;

	    CursorAfterString(infile,"Enter string spring constant:");
            fscanf(infile,"%lf",&af_params->kl);
            (void) printf("%f\n",af_params->kl);
            
            CursorAfterString(infile,"Enter string damping constant:");
            fscanf(infile,"%lf",&af_params->lambda_l);
            (void) printf("%f\n",af_params->lambda_l);
            
            CursorAfterString(infile,"Enter string friction constant:");
            fscanf(infile,"%lf",&af_params->mu_l);
            (void) printf("%f\n",af_params->mu_l);
            
            if (CursorAfterStringOpt(infile,"Enter string thickness:"))
            {
                fscanf(infile,"%lf",&af_params->string_thickness);
                (void) printf("%f\n",af_params->string_thickness);
            }

            if (CursorAfterStringOpt(infile,"Enter string rounding tolerance:"))
            {
                fscanf(infile,"%lf",&af_params->string_eps);
                (void) printf("%f\n",af_params->string_eps);
            }

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


    //For collision solver inelastic impulse control
    if (CursorAfterStringOpt(infile,"Enter collision inelastic impulse relaxation parameter:"))
    {
        fscanf(infile,"%lf",&af_params->inelastic_impulse_coeff);
        (void) printf("%f\n",af_params->inelastic_impulse_coeff);
    }

    //For collision solver elastic impulse control
    if (CursorAfterStringOpt(infile,"Enter overlap coefficient:"))
    {
        fscanf(infile,"%lf",&af_params->overlap_coefficient);
        (void) printf("%f\n",af_params->overlap_coefficient);
    }
    
    //strain limiting parameters
    if (CursorAfterStringOpt(infile,"Enter strain limit:"))
    {
        fscanf(infile,"%lf",&af_params->strain_limit);
        (void) printf("%f\n",af_params->strain_limit);
    }

    if (CursorAfterStringOpt(infile,"Enter compressive strain limit:"))
    {
        fscanf(infile,"%lf",&af_params->compressive_strain_limit);
        (void) printf("%f\n",af_params->compressive_strain_limit);
    }

    if (CursorAfterStringOpt(infile,"Enter strain rate limit:"))
    {
        fscanf(infile,"%lf",&af_params->strainrate_limit);
        (void) printf("%f\n",af_params->strainrate_limit);
    }

    if (CursorAfterStringOpt(infile,"Enter strain velocity constraint tolerance:"))
    {
        fscanf(infile,"%lf",&af_params->strain_vel_tol);
        (void) printf("%f\n",af_params->strain_vel_tol);
    }


	if (dim == 3 && af_params->is_parachute_system == YES)
	{
	    af_params->m_g = af_params->m_s;
        if (af_params->attach_gores == YES)
	    {
            af_params->gores_present = true;

		    CursorAfterString(infile,"Enter gore spring constant:");
        	fscanf(infile,"%lf",&af_params->kg);
        	(void) printf("%f\n",af_params->kg);
        	CursorAfterString(infile,"Enter gore friction constant:");
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
	    
        af_params->unequal_coeff = 1.0;
	    if (CursorAfterStringOpt(infile,"Enter unequal coefficient:"))
	    {
		    fscanf(infile,"%lf",&af_params->unequal_coeff);
		    (void) printf("%f\n",af_params->unequal_coeff);
	    }
	    
        af_params->unequal_strings_num = 0;
	    if (CursorAfterStringOpt(infile,
				"Enter number of unequal strings:"))
	    {
            fscanf(infile,"%d",&af_params->unequal_strings_num);
            (void) printf("%d\n",af_params->unequal_strings_num);
	    }

        if (af_params->unequal_strings_num > 0)
	    {
    		FT_VectorMemoryAlloc((POINTER*)&af_params->unequal_strings_gindex,
                    af_params->unequal_strings_num,sizeof(int));
            
            int *us_gindex = af_params->unequal_strings_gindex;
            
            if (CursorAfterStringOpt(infile,"Enter ID of unequal strings:"))
            {
                for (int i = 0; i < af_params->unequal_strings_num; ++i)
                {
                    fscanf(infile,"%d",&us_gindex[i]);
                    (void) printf("%d ",us_gindex[i]);
                }
                (void) printf("\n");
            }
            else
            {
                for (int i = 0; i < af_params->unequal_strings_num; ++i)
                    us_gindex[i] = i;
            }
	    }
	    
        af_params->break_strings_time = -1.0;
	    if (CursorAfterStringOpt(infile,"Enter time to break strings:"))
	    {
            fscanf(infile,"%lf",&af_params->break_strings_time);
            (void) printf("%f\n",af_params->break_strings_time);
	    }
	    
        af_params->break_strings_num = 0;
	    if (CursorAfterStringOpt(infile,"Enter number of strings to break:"))
	    {
            fscanf(infile,"%d",&af_params->break_strings_num);
            (void) printf("%d\n",af_params->break_strings_num);
	    }
	    
        if (af_params->break_strings_num > 0)
	    {
            FT_VectorMemoryAlloc((POINTER*)&af_params->break_strings_gindex,
                    af_params->break_strings_num,sizeof(int));
		
            int *bs_gindex = af_params->break_strings_gindex;
		
            if (CursorAfterStringOpt(infile,"Enter ID of strings to break:"))
            {
                for (int i = 0; i < af_params->break_strings_num; ++i)
                {
                    fscanf(infile,"%d",&bs_gindex[i]);
                    (void) printf("%d ",bs_gindex[i]);
                }
                (void) printf("\n");
            }
            else
            {
                for (int i = 0; i < af_params->break_strings_num; ++i)
                    bs_gindex[i] = i;
            }
	    }
    }

    /*
    af_params->num_smooth_layers = 1;
	if (CursorAfterStringOpt(infile,"Enter number of smooth layers:"))
	{
        fscanf(infile,"%d",&af_params->num_smooth_layers);
        (void) printf("%d\n",af_params->num_smooth_layers);
	}
    */
	
    if (af_params->is_parachute_system == YES && !af_params->rgb_payload)
    {
        printf("setMotionParams(): pointmass payload detected\n");
        CursorAfterString(infile,"Enter payload:");
        fscanf(infile,"%lf",&af_params->payload);
        (void) printf("%f\n",af_params->payload);
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

	printf("canopy points count (fabric+gore)  = %d, "
		"string points count = %d\n", countSurfPoints(front->interf),
        countStringPoints(front->interf, af_params->is_parachute_system));

	printf("fabric point mass  = %f, string point mass = %f, "
		"gore point mass = %f\n", af_params->m_s, af_params->m_l, af_params->m_g);

	fclose(infile);
}	/* end setFabricParams */


static int countSurfPoints(INTERFACE* intfc)
{
	int surf_count = 0;
	int num_fabric_pts = 0;
	
    for (SURFACE** s = intfc->surfaces; s && *s; ++s)
    {
        if (wave_type(*s) == ELASTIC_BOUNDARY)
        {
            num_fabric_pts += I_NumOfSurfPoints(*s);
            surf_count++;
        }
    }
	
    return num_fabric_pts;
}

//TODO: Need to adapt for general case -- not just for pointmass
//      and single attachement point on forebody for all strings
static int countStringPoints(INTERFACE* intfc, boolean is_parachute_system)
{
	int num_str_pts = 0;
    for (CURVE** c = intfc->curves; c && *c; ++c)
    {
        if (FT_Dimension() == 3 && hsbdry_type(*c) == STRING_HSBDRY)
            num_str_pts += I_NumOfCurvePoints(*c);
        else if (FT_Dimension() == 2 && wave_type(*c) == ELASTIC_STRING)
            num_str_pts += I_NumOfCurvePoints(*c);
        else
            continue;
    
        if (is_parachute_system == YES)
            num_str_pts -= 2; //exclude curve boundaries
    }
	

    //TODO: account for non point load runs -- add number of rg_string_nodes instead
    if (is_parachute_system == YES)
	    num_str_pts += 1; //load node
	
    return num_str_pts;
}

static void initVelocityFunc(
	FILE *infile,
	Front *front)
{
    // velocity function parameters
	static VELO_FUNC_PACK velo_func_pack;
	static VORTEX_PARAMS *vortex_params;
    static BIPOLAR_PARAMS *dv_params;
	static VERTICAL_PARAMS *vert_params;
	static RANDOMV_PARAMS *randv_params;
	static TOROIDAL_PARAMS *toro_params;
	static PARABOLIC_PARAMS *para_params;
	static SINGULAR_PARAMS *sing_params;
	static FIXAREA_PARAMS *fixarea_params;
    
    EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int dim = front->rect_grid->dim;
	char string[100];

	if (af_params->no_fluid == YES)
	{
	    front->curve_propagate = airfoil_curve_propagate;
	    front->node_propagate = airfoil_node_propagate;
	    velo_func_pack.point_propagate = airfoil_point_propagate;
	    (void) printf("Available velocity functions are:\n");
	    (void) printf("\tVortex velocity (R)\n");
	    (void) printf("\tDouble vortex velocity (D)\n");
	    (void) printf("\tVertical velocity (V)\n");
	    (void) printf("\tToroidal velocity (T)\n");
	    (void) printf("\tParabolic velocity (P)\n");
	    (void) printf("\tSingular velocity (S)\n");
	    (void) printf("\tZero velocity (Z)\n");
	    (void) printf("\tFixed area velocity (FA)\n");
	    (void) printf("\tFixed point velocity (FP)\n");
	    (void) printf("\tFree fall velocity (FF)\n");
   
        CursorAfterString(infile,"Enter velocity function:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    
        switch (string[0])
        {
        case 'r':
        case 'R':
            if (string[1] == 'o' || string[1] == 'O')
            {
                FT_ScalarMemoryAlloc((POINTER*)&vortex_params,
                            sizeof(VORTEX_PARAMS));
                front->max_time = 0.4;
                front->movie_frame_interval = 0.02;
                vortex_params->dim = 2;
                vortex_params->type[0] = 'M';
                vortex_params->cos_time = 0;
                vortex_params->cen[0] = 0.5;
                vortex_params->cen[1] = 0.25;
                vortex_params->rad = 0.15;
                vortex_params->time = 0.5*front->max_time;
                velo_func_pack.func_params = (POINTER)vortex_params;
                velo_func_pack.func = vortex_vel;
            }
            else if (string[1] == 'a' || string[1] == 'A')
            {
                FT_ScalarMemoryAlloc((POINTER*)&randv_params,
                        sizeof(RANDOMV_PARAMS));
                CursorAfterString(infile,"Enter random amplitude:");
                for (int i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&randv_params->v0[i]);
                    (void) printf("%fi ",randv_params->v0[i]);
                }
                (void) printf("\n");
                
                CursorAfterString(infile,"Enter stop motion time:");
                fscanf(infile,"%lf",&randv_params->stop_time);
                (void) printf("%f\n",randv_params->stop_time);
                velo_func_pack.func_params = (POINTER)randv_params;
                velo_func_pack.func = random_velo;
            }
            else
            {
                (void) printf("ERROR: need either RO or RA\n");
                LOC(); clean_up(ERROR);
            }
            break;
            
        case 'd':
        case 'D':
	    	FT_ScalarMemoryAlloc((POINTER*)&dv_params, sizeof(BIPOLAR_PARAMS));
            dv_params->cen1[0] = 0.25;
            dv_params->cen1[1] = 0.25;
            dv_params->cen2[0] = 0.75;
            dv_params->cen2[1] = 0.25;
            dv_params->i1 = -0.5;
            dv_params->i2 =  0.5;
            velo_func_pack.func_params = (POINTER)dv_params;
            velo_func_pack.func = double_vortex_vel;
            break;
            
        case 'v':
        case 'V':
	    	FT_ScalarMemoryAlloc((POINTER*)&vert_params,sizeof(VERTICAL_PARAMS));
            CursorAfterString(infile,"Enter center velocity:");
                fscanf(infile,"%lf",&vert_params->v0);
                (void) printf("%f\n",vert_params->v0);
            CursorAfterString(infile,"Enter stop motion time:");
                fscanf(infile,"%lf",&vert_params->stop_time);
                (void) printf("%f\n",vert_params->stop_time);
            CursorAfterString(infile,"Enter center of vertical motion:");
                fscanf(infile,"%lf %lf",&vert_params->cen[0],&vert_params->cen[1]);
                (void) printf("%f %f\n",vert_params->cen[0],vert_params->cen[1]);
            velo_func_pack.func_params = (POINTER)vert_params;
            velo_func_pack.func = vertical_velo;
            break;
            
        case 't':
        case 'T':
	    	FT_ScalarMemoryAlloc((POINTER*)&toro_params,sizeof(TOROIDAL_PARAMS));
            CursorAfterString(infile,"Enter center of toroidal motion:");
                fscanf(infile,"%lf %lf %lf",&toro_params->tcen[0],
                    &toro_params->tcen[1],&toro_params->tcen[2]);
                (void) printf("%f %f %f\n",toro_params->tcen[0],
                    toro_params->tcen[1],toro_params->tcen[2]);
            CursorAfterString(infile,"Enter distance to poloidal center:");
                fscanf(infile,"%lf",&toro_params->R0);
                (void) printf("%f\n",toro_params->R0);
            CursorAfterString(infile,"Enter velocity magnitude:");
                fscanf(infile,"%lf",&toro_params->v0);
                (void) printf("%f\n",toro_params->v0);
            CursorAfterString(infile,"Enter stop motion time:");
                fscanf(infile,"%lf",&toro_params->stop_time);
                (void) printf("%f\n",toro_params->stop_time);
            velo_func_pack.func_params = (POINTER)toro_params;
            velo_func_pack.func = toroidal_velo;
            break;
            
        case 'p':
        case 'P':
	    	FT_ScalarMemoryAlloc((POINTER*)&para_params,sizeof(PARABOLIC_PARAMS));
            CursorAfterString(infile,"Enter center of parabolic velocity:");
                fscanf(infile,"%lf %lf",&para_params->cen[0],
                    &para_params->cen[1]);
                (void) printf("%f %f\n",para_params->cen[0],
                    para_params->cen[1]);
            CursorAfterString(infile,"Enter center velocity:");
                fscanf(infile,"%lf",&para_params->v0);
                (void) printf("%f\n",para_params->v0);
            CursorAfterString(infile,"Enter downward concavity:");
                fscanf(infile,"%lf",&para_params->a);
                (void) printf("%f\n",para_params->a);
            CursorAfterString(infile,"Enter stop motion time:");
                fscanf(infile,"%lf",&para_params->stop_time);
                (void) printf("%f\n",para_params->stop_time);
            velo_func_pack.func_params = (POINTER)para_params;
            velo_func_pack.func = parabolic_velo;
            break;
            
        case 's':
        case 'S':
	    	FT_ScalarMemoryAlloc((POINTER*)&sing_params,sizeof(SINGULAR_PARAMS));
            CursorAfterString(infile,"Enter center of velocity:");
                fscanf(infile,"%lf %lf",&sing_params->cen[0],&sing_params->cen[1]);
                (void) printf("%f %f\n",sing_params->cen[0],sing_params->cen[1]);
            CursorAfterString(infile,"Enter center velocity:");
                fscanf(infile,"%lf",&sing_params->v0);
                (void) printf("%f\n",sing_params->v0);
            CursorAfterString(infile,"Enter radius of center:");
                fscanf(infile,"%lf",&sing_params->R);
                (void) printf("%f\n",sing_params->R);
            CursorAfterString(infile,"Enter stop motion time:");
                fscanf(infile,"%lf",&sing_params->stop_time);
                (void) printf("%f\n",sing_params->stop_time);
            velo_func_pack.func_params = (POINTER)sing_params;
            velo_func_pack.func = singular_velo;
            break;
            
        case 'f':
        case 'F':
	    	FT_ScalarMemoryAlloc((POINTER*)&fixarea_params,sizeof(FIXAREA_PARAMS));
            if (string[1] == 'a' || string[1] == 'A')
            {
                if (!front->f_basic->RestartRun)
                    init_fixarea_params(front,infile,fixarea_params);
                else
                    restart_fixarea_params(front,infile,fixarea_params);
            }
            else if (string[1] == 'p' || string[1] == 'P')
            {
                init_fixpoint_params(front,infile,fixarea_params);
            }
            else if (string[1] == 'f' || string[1] == 'F')
            {
                velo_func_pack.func_params = NULL;
                velo_func_pack.func = NULL;
            }
        
            velo_func_pack.func_params = (POINTER)fixarea_params;
            velo_func_pack.func = marker_velo;
            break;
            
        case 'z':
        case 'Z':
            velo_func_pack.func_params = NULL;
            velo_func_pack.func = zero_velo;
            break;
            
        default:
            printf("ERROR initVelocityFunc(): unknown velocity function!\n");
            LOC(); clean_up(EXIT_FAILURE);
            break;
        }
	}

    if (CursorAfterStringOpt(infile,"Enter gravity:"))
    {
        for (int i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",af_params->gravity+i);
            printf(" %f",af_params->gravity[i]);
            eqn_params->gravity[i] = af_params->gravity[i];
        }
        printf("\n");
    }

	FT_InitFrontVeloFunc(front,&velo_func_pack);
}	/* end initVelocityFunc */

static int zero_velo(
	POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	vel[0] = vel[1] = vel[2] = 0.0;
	return YES;
}	/* end zero_velo */

static int random_velo(
	POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	RANDOMV_PARAMS *randv_params = (RANDOMV_PARAMS*)params;
	double *v0 = randv_params->v0;
	double stop_time = randv_params->stop_time;
	unsigned short int xsubi[3];

	if (front->time >= stop_time)
	{
	    vel[0] = vel[1] = vel[2] = 0.0;
	}
	else
	{
	    xsubi[0] = 7256;
	    xsubi[1] = 764;
	    xsubi[2] = 2163;
	    vel[0] = v0[0]*(2.0*erand48(xsubi) - 1.0);
	    vel[1] = v0[1]*(2.0*erand48(xsubi) - 1.0);
	    vel[2] = v0[2]*(2.0*erand48(xsubi) - 1.0);
	}
	return YES;
}	/* end random_velo */

static int vertical_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	VERTICAL_PARAMS *vert_params = (VERTICAL_PARAMS*)params;
    double* cen = vert_params->cen;
	double v0 = vert_params->v0;

    double *coords = Coords(p);
	double dist = sqrt(sqr(coords[0] - cen[0]) + sqr(coords[1] - cen[1]));
	double v_vert = (0.15 - dist)/0.15;

	vel[0] = vel[1] = 0.0;
	double stop_time = vert_params->stop_time;
	if (front->time < stop_time)
	    vel[2] = v_vert*v0;
	else
	    vel[2] = 0.0;
	return YES;
}       /* end vertical_velo */

static int toroidal_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	TOROIDAL_PARAMS *toro_params = (TOROIDAL_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = toro_params->v0;
	double stop_time = toro_params->stop_time;
	double tcoords[2]; 		/* toroidal coords */
	double pvel[2];		/* projected velocity */
	double *tcen = toro_params->tcen; /* toroidal center of vel */
	double R0 = toro_params->R0; /* radial dist of poloidal center */
	double d1,d2;
	double s1,s2;
	double dx1,dy1;
	double dx2,dy2;
	double dtol = 0.000001*front->rect_grid->h[0];

	if (front->time >= stop_time)
	{
	    for (i = 0; i < dim; ++i)
		vel[i] = 0.0;
	    return YES;
	}
	/* Project 3D to 2D */
	tcoords[0] = sqrt(sqr(coords[0] - tcen[0]) + sqr(coords[1] - tcen[1]));
	tcoords[1] = coords[2] - tcen[2];

	dx1 = tcoords[0] - R0;
	dx2 = tcoords[0] + R0;
	dy1 = dy2 = tcoords[1];

	d1 = sqr(dx1) + sqr(dy1);
	d2 = sqr(dx2) + sqr(dy2);
	s1 = v0;
	s2 = -v0;
	if (d1 < dtol || d2 < dtol)
	    pvel[0] = pvel[1] = 0.0;
	else if (front->time < stop_time)
	{
	    pvel[0] =  s1*dy1/d1 + s2*dy2/d2;
	    pvel[1] = -s1*dx1/d1 - s2*dx2/d2;
	} 
	else
	    pvel[0] = pvel[1] = 0.0;
	    
	/* give it back to 3D */
	vel[0] = pvel[0]*(coords[0]-tcen[0])/tcoords[0];
	vel[1] = pvel[0]*(coords[1]-tcen[1])/tcoords[0];
	vel[2] = pvel[1];
	return YES;
}       /* end toroidal_velo */

static int parabolic_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	PARABOLIC_PARAMS *para_params = (PARABOLIC_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = para_params->v0;
	double a = para_params->a;
	double *cen = para_params->cen;
	double R_sqr = 0.0;

	for (i = 0; i < dim-1; ++i)
	{
	    R_sqr += sqr(coords[i] - cen[i]);
	    vel[i] = 0.0;
	}
	vel[dim-1] = -0.5*a*R_sqr + v0;
	return YES;
}	/* end parabolic_velo */

static int singular_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	SINGULAR_PARAMS *para_params = (SINGULAR_PARAMS*)params;
        double *coords = Coords(p);
	int i,dim = front->rect_grid->dim;
	double v0 = para_params->v0;
	double R = para_params->R;
	double *cen = para_params->cen;
	double r = 0.0;
	for (i = 0; i < dim-1; ++i)
	{
	    r += sqr(coords[i] - cen[i]);
	    vel[i] = 0.0;
	}
	r = sqrt(r);
	if (r < R)
	    vel[dim-1] = v0;
	else
	    vel[dim-1] = 0.0;
	return YES;
}	/* end sigular_velo */

static void init_fixarea_params(
	Front *front,
	FILE *infile,
	FIXAREA_PARAMS *fixarea_params)
{
	char string[100];
	SHAPE_PARAMS sparams;
	int num_pts,dim = front->rect_grid->dim;
	static REGISTERED_PTS *registered_pts;
    SURFACE **s;
    TRI *tri;
	POINT *p;
	INTERFACE *intfc = front->interf;

	(void) printf("Available initial areas are:\n");
	(void) printf("\tRectangle (R)\n");
	(void) printf("\tEllipse (E)\n");
	
    CursorAfterString(infile,"Enter initial shape of fixed area:");
    fscanf(infile,"%s",string);
    (void) printf("%s\n",string);
	
    switch (string[0])
	{
	case 'r':
	case 'R':
	    sparams.shape_id = 0;
	    CursorAfterString(infile,"Enter rectangle lower bounds:");
            fscanf(infile,"%lf %lf",&sparams.L[0],&sparams.L[1]);
            (void) printf("%f %f\n",sparams.L[0],sparams.L[1]);
	    CursorAfterString(infile,"Enter rectangle upper bounds:");
            fscanf(infile,"%lf %lf",&sparams.U[0],&sparams.U[1]);
            (void) printf("%f %f\n",sparams.U[0],sparams.U[1]);
	    break;
	case 'e':
	case 'E':
	    sparams.shape_id = 1;
	    CursorAfterString(infile,"Enter center of ellipse:");
            fscanf(infile,"%lf %lf",&sparams.cen[0],&sparams.cen[1]);
            (void) printf("%f %f\n",sparams.cen[0],sparams.cen[1]);
	    CursorAfterString(infile,"Enter radii of ellipse:");
            fscanf(infile,"%lf %lf",&sparams.R[0],&sparams.R[1]);
            (void) printf("%f %f\n",sparams.R[0],sparams.R[1]);
	    break;
	}
	
    CursorAfterString(infile,"Enter area velocity:");
    for (int i = 0; i < dim; ++i)
    {
        fscanf(infile,"%lf",&fixarea_params->vel[i]);
        (void) printf("%f ",fixarea_params->vel[i]);
    }
    (void) printf("\n");

	num_pts = 0;
    if (!front->f_basic->RestartRun)
    {
        /* Count number of registered points */
        num_pts = 0;
        reset_sort_status(intfc);
        intfc_surface_loop(intfc,s)
        {
            if (wave_type(*s) != ELASTIC_BOUNDARY &&
                wave_type(*s) != ELASTIC_STRING) continue;

            surf_tri_loop(*s,tri)
            {
                for (int i = 0; i < 3; ++i)
                {
                    p = Point_of_tri(tri)[i];
                    if (sorted(p)) continue;
            
                    sorted(p) = YES;
                    if (within_shape(sparams,Coords(p))) num_pts++;	
                }
            }
        }

        FT_ScalarMemoryAlloc((POINTER*)&registered_pts,
                sizeof(REGISTERED_PTS));
        FT_VectorMemoryAlloc((POINTER*)&registered_pts->global_ids,
                num_pts,sizeof(int));
        
        registered_pts->num_pts = num_pts;

        /* Record registered points */
        num_pts = 0;	
        reset_sort_status(intfc);
        intfc_surface_loop(intfc,s)
        {
            if (wave_type(*s) != ELASTIC_BOUNDARY &&
                wave_type(*s) != ELASTIC_STRING) continue;
        
            (*s)->extra = (POINTER)registered_pts;
            surf_tri_loop(*s,tri)
            {
                for (int i = 0; i < 3; ++i)
                {
                    p = Point_of_tri(tri)[i];
                    if (sorted(p)) continue;
                    
                    sorted(p) = YES;
                    if (within_shape(sparams,Coords(p)))
                    {
                        registered_pts->global_ids[num_pts] = Gindex(p);	
                        num_pts++;	
                    }
                }
            }
        }
    }
    /*else
    {
        //Restart runs
        intfc_surface_loop(intfc,s)
        {
            if (wave_type(*s) != ELASTIC_BOUNDARY &&
                wave_type(*s) != ELASTIC_STRING) continue;
            if ((*s)->extra == NULL) continue;
            registered_pts = (REGISTERED_PTS*)(*s)->extra;
            num_pts = registered_pts->num_pts;
            break;
        }
    }*/
    
    if (num_pts == 0) return;

    fixarea_params->num_pts = num_pts;
    FT_VectorMemoryAlloc((POINTER*)&fixarea_params->global_ids,num_pts,sizeof(int));

    for (int i = 0; i < num_pts; ++i)
    {
        fixarea_params->global_ids[i] = registered_pts->global_ids[i];	
    }
}	/* end init_fixarea_params */

static void restart_fixarea_params(
    Front* front,
    FILE *infile,
    FIXAREA_PARAMS *fixarea_params)
{
    REGISTERED_PTS* registered_pts;
    INTERFACE *intfc = front->interf;
    SURFACE **s;
    
    intfc_surface_loop(intfc,s)
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_STRING) continue;
        
        if ((*s)->extra == NULL) continue;
	    
        registered_pts = (REGISTERED_PTS*)(*s)->extra;
        break;
    }

	int num_pts = registered_pts->num_pts;
    fixarea_params->num_pts = num_pts;
    FT_VectorMemoryAlloc((POINTER*)&fixarea_params->global_ids,num_pts,sizeof(int));

    for (int i = 0; i < num_pts; ++i)
    {
        fixarea_params->global_ids[i] = registered_pts->global_ids[i];	
    }
}
    

/*
//TODO: Replacing with version that supports restart runs.
//      Remove when restart runs are stable.
//
static void init_fixarea_params(
	Front *front,
	FILE *infile,
	FIXAREA_PARAMS *fixarea_params)
{
	char string[100];
	SHAPE_PARAMS sparams;
	int i,num_pts,dim = front->rect_grid->dim;
	static REGISTERED_PTS *registered_pts;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	SURFACE *surf;

	(void) printf("Available initial areas are:\n");
	(void) printf("\tRectangle (R)\n");
	(void) printf("\tEllipse (E)\n");
	CursorAfterString(infile,"Enter initial shape of fixed area:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	switch (string[0])
	{
	case 'r':
	case 'R':
	    sparams.shape_id = 0;
	    CursorAfterString(infile,"Enter rectangle lower bounds:");
            fscanf(infile,"%lf %lf",&sparams.L[0],&sparams.L[1]);
            (void) printf("%f %f\n",sparams.L[0],sparams.L[1]);
	    CursorAfterString(infile,"Enter rectangle upper bounds:");
            fscanf(infile,"%lf %lf",&sparams.U[0],&sparams.U[1]);
            (void) printf("%f %f\n",sparams.U[0],sparams.U[1]);
	    break;
	case 'e':
	case 'E':
	    sparams.shape_id = 1;
	    CursorAfterString(infile,"Enter center of ellipse:");
            fscanf(infile,"%lf %lf",&sparams.cen[0],&sparams.cen[1]);
            (void) printf("%f %f\n",sparams.cen[0],sparams.cen[1]);
	    CursorAfterString(infile,"Enter radii of ellipse:");
            fscanf(infile,"%lf %lf",&sparams.R[0],&sparams.R[1]);
            (void) printf("%f %f\n",sparams.R[0],sparams.R[1]);
	    break;
	}
	CursorAfterString(infile,"Enter area velocity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&fixarea_params->vel[i]);
            (void) printf("%f ",fixarea_params->vel[i]);
        }
        (void) printf("\n");
	
    num_pts = 0;
	next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
        if (wave_type(hs) != ELASTIC_BOUNDARY &&
            wave_type(hs) != ELASTIC_STRING) continue;
        if (within_shape(sparams,Coords(p))) num_pts++;	
    }
	FT_VectorMemoryAlloc((POINTER*)&fixarea_params->global_ids,num_pts,
				sizeof(int));
	FT_ScalarMemoryAlloc((POINTER*)&registered_pts,sizeof(REGISTERED_PTS));
	FT_VectorMemoryAlloc((POINTER*)&registered_pts->global_ids,num_pts,
				sizeof(int));
	fixarea_params->num_pts = num_pts;
	registered_pts->num_pts = num_pts;

	num_pts = 0;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            //TODO: Add hsbdry_type() == STRING_HSBDRY, GORE_HSBDRY etc.
            //      to facilitate other kinds initialization procedures.
            if (wave_type(hs) != ELASTIC_BOUNDARY &&
		wave_type(hs) != ELASTIC_STRING) 
		continue;
	    if (within_shape(sparams,Coords(p)))
	    {
		fixarea_params->global_ids[num_pts] = Gindex(p);	
		registered_pts->global_ids[num_pts] = Gindex(p);	
		num_pts++;
	    	surf = Surface_of_hs(hs);
	    	surf->extra = (POINTER)registered_pts;
	    }
	}
}*/	/* end init_fixarea_params */

static boolean within_shape(
	SHAPE_PARAMS sparams,
	double *coords)
{
	double dist;
	switch (sparams.shape_id)
	{
	case 0:
	    if (coords[0] < sparams.L[0]) return NO;
	    if (coords[0] > sparams.U[0]) return NO;
	    if (coords[1] < sparams.L[1]) return NO;
	    if (coords[1] > sparams.U[1]) return NO;
	    return YES;
	case 1:
	    dist = sqr((coords[0] - sparams.cen[0])/sparams.R[0]) +
	           sqr((coords[1] - sparams.cen[1])/sparams.R[1]);
	    if (dist < 1.0) return YES;
	    else return NO;
	}
}	/* within_shape */

static void init_fixpoint_params(
	Front *front,
	FILE *infile,
	FIXAREA_PARAMS *fixarea_params)
{
	int i,dim = front->rect_grid->dim;
	static REGISTERED_PTS *registered_pts;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;
	SURFACE *surf;
	double coords[MAXD];
	double dist,min_dist;

	FT_VectorMemoryAlloc((POINTER*)&fixarea_params->global_ids,1,
				sizeof(int));
	FT_ScalarMemoryAlloc((POINTER*)&registered_pts,sizeof(REGISTERED_PTS));
	FT_VectorMemoryAlloc((POINTER*)&registered_pts->global_ids,1,
				sizeof(int));
	fixarea_params->num_pts = 1;
	registered_pts->num_pts = 1;

	CursorAfterString(infile,"Enter initial fixed point coordinates:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&coords[i]);
            (void) printf("%f ",coords[i]);
        }
        (void) printf("\n");
	CursorAfterString(infile,"Enter point velocity:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&fixarea_params->vel[i]);
            (void) printf("%f ",fixarea_params->vel[i]);
        }
        (void) printf("\n");

	min_dist = HUGE;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY &&
		wave_type(hs) != ELASTIC_STRING) 
		continue;
	    dist = distance_between_positions(coords,Coords(p),dim);
	    if (dist < min_dist)
	    {
		min_dist = dist;
		registered_pts->global_ids[0] = Gindex(p);
		fixarea_params->global_ids[0] = Gindex(p);
	    	surf = Surface_of_hs(hs);
	    	surf->extra = (POINTER)registered_pts;
	    }
        }
}	/* end init_fixpoint_params */

static int marker_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	FIXAREA_PARAMS *fixarea_params = (FIXAREA_PARAMS*)params;
	int dim = front->rect_grid->dim;
	int num_pts = fixarea_params->num_pts;
	int *global_ids = fixarea_params->global_ids;
	
    for (int i = 0; i < num_pts; ++i)
	{
	    if (Gindex(p) == global_ids[i])
	    {
            for (int j = 0; j < dim; ++j)
            {
                vel[j] = fixarea_params->vel[j];
            }
            return YES;
	    }
	}
	
    for (int j = 0; j < dim; ++j) vel[j] = 0.0;
    return YES;
}	/* end marker_velo */

extern void resetFrontVelocity(Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	
    if (CursorAfterStringOpt(infile,"Enter yes to reset front velocity:"))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
    }
    fclose(infile);
        
    if (string[0] == 'y' || string[0] == 'Y')
        zeroFrontVelocity(front);
}

extern void zeroFrontVelocity(Front *front)
{
	INTERFACE *intfc = front->interf;
	POINT *p;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	STATE *sl,*sr;
	int i,dim = front->rect_grid->dim;
	CURVE **c;
	BOND *b;

	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    for (i = 0; i < dim; ++i)
	    {
            p->vel[i] = 0.0;
            p->force[i] = 0.0;
            sl->vel[i] = sr->vel[i] = 0.0;
            sl->impulse[i] = sr->impulse[i] = 0.0;
            sl->fluid_accel[i] = sr->fluid_accel[i] = 0.0;
            sl->other_accel[i] = sr->other_accel[i] = 0.0;
            sl->shear_force[i] = sr->shear_force[i] = 0.0;
            sl->bendforce[i] = sr->bendforce[i] = 0.0;
	    }
	}

	if (dim == 3)
	{
	    for (c = intfc->curves; c && *c; ++c)
	    {
            p = (*c)->start->posn;
            sl = (STATE*)left_state(p);
            sr = (STATE*)right_state(p);
	        for (i = 0; i < dim; ++i)
            {
                p->vel[i] = 0.0;
                p->force[i] = 0.0;
                sl->vel[i] = sr->vel[i] = 0.0;
                sl->impulse[i] = sr->impulse[i] = 0.0;
                sl->fluid_accel[i] = sr->fluid_accel[i] = 0.0;
                sl->other_accel[i] = sr->other_accel[i] = 0.0;
                sl->bendforce[i] = sr->bendforce[i] = 0.0;
            }
    
            for (b = (*c)->first; b != (*c)->last; b = b->next)
            {
                p = b->end;
                sl = (STATE*)left_state(p);
                sr = (STATE*)right_state(p);
                    for (i = 0; i < dim; ++i)
                {
                    p->vel[i] = 0.0;
                    p->force[i] = 0.0;
                    sl->vel[i] = sr->vel[i] = 0.0;
                    sl->impulse[i] = sr->impulse[i] = 0.0;
                    sl->fluid_accel[i] = sr->fluid_accel[i] = 0.0;
                    sl->other_accel[i] = sr->other_accel[i] = 0.0;
                    sl->bendforce[i] = sr->bendforce[i] = 0.0;
                }
            }
    
            p = (*c)->end->posn;
            sl = (STATE*)left_state(p);
            sr = (STATE*)right_state(p);
	        for (i = 0; i < dim; ++i)
            {
                p->vel[i] = 0.0;
                p->force[i] = 0.0;
                sl->vel[i] = sr->vel[i] = 0.0;
                sl->impulse[i] = sr->impulse[i] = 0.0;
                sl->fluid_accel[i] = sr->fluid_accel[i] = 0.0;
                sl->other_accel[i] = sr->other_accel[i] = 0.0;
                sl->bendforce[i] = sr->bendforce[i] = 0.0;
            }
	    }
	}
}	/* end resetFrontVelocity */

//TODO: check if general enough for our purposes or is just a hardcoded special case...
static void convert_to_point_mass(
        Front *front,
        AF_PARAMS *af_params)
{
        INTERFACE *intfc;
        int num_str_pts, num_fabric_pts, num_gore_pts;
        SURFACE **s;
        CURVE **c;
        intfc = front->interf;
        int dim = Dimension(intfc);

        switch (dim)
        {
        case 2:
            num_str_pts = num_fabric_pts = 0;
            for (c = intfc->curves; c && *c; ++c)
            {
                if (wave_type(*c) == ELASTIC_BOUNDARY)
                    num_fabric_pts +=  I_NumOfCurvePoints(*c);
                else if (wave_type(*c) == ELASTIC_STRING)
                {
                    num_str_pts += I_NumOfCurvePoints(*c);
                    if (af_params->is_parachute_system == YES)
                        num_str_pts -= 2; //exclude curve boundary
                }
            }
    
            if (af_params->is_parachute_system == YES)
		        num_str_pts += 1; //load node
	    
            if (num_fabric_pts != 0)
                af_params->m_s = af_params->total_canopy_mass/num_fabric_pts;
            else
                af_params->m_s = 0.001;
            
            if (num_str_pts != 0)
                af_params->m_l = af_params->total_string_mass/num_str_pts;
            else
                af_params->m_l = 0.002;

            break;
        
        case 3:
            num_str_pts = num_fabric_pts = num_gore_pts = 0;
            for (s = intfc->surfaces; s && *s; ++s)
            {
                if (wave_type(*s) == ELASTIC_BOUNDARY)
                    num_fabric_pts += I_NumOfSurfPoints(*s);
            }
            
            for (c = intfc->curves; c && *c; ++c)
            {
                if (hsbdry_type(*c) == STRING_HSBDRY)
                {
                    num_str_pts += I_NumOfCurvePoints(*c); 
                    if (af_params->is_parachute_system == YES)
                        num_str_pts -= 2; //exclude curve boundary
                }
                else if (hsbdry_type(*c) == GORE_HSBDRY)
                {
                    num_gore_pts += I_NumOfCurvePoints(*c);
                }
            }
    
            if (af_params->is_parachute_system == YES)
	        	num_str_pts += 1; //load node
	    
            num_fabric_pts -= num_gore_pts;
	        if (num_fabric_pts != 0)
		        af_params->m_s = af_params->total_canopy_mass/num_fabric_pts;
            else
		        af_params->m_s = 0.001;
            
            if (num_str_pts != 0)
                af_params->m_l = af_params->total_string_mass/num_str_pts;
            else
                af_params->m_l = 0.002;
	    
            if (num_gore_pts != 0)
		        af_params->m_g = af_params->total_gore_mass/num_gore_pts;
	        else
		        af_params->m_g = 0.001;
	    
            break;
        }
}       /* end convert_to_point_mass */

//TODO: Should move to afsetd.cpp
static void checkSetGoreNodes(
	INTERFACE *intfc)
{
	CURVE **c;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY)
	    {
            set_gore_node((*c)->start);
            set_gore_node((*c)->end);
	    }
	}
}	/* end checkSetGoreNodes */


//TODO: Should move to afsetd.cpp
static void set_gore_node(
	NODE *n)
{
	static AF_NODE_EXTRA *extra;
	boolean is_gore_node = NO;
	CURVE **c;

	for (c = n->in_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return;
	    else if (hsbdry_type(*c) == GORE_HSBDRY)
		is_gore_node = YES;
	for (c = n->out_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return;
	    else if (hsbdry_type(*c) == GORE_HSBDRY)
		is_gore_node = YES;
	if (!is_gore_node) return;
	if (n->extra == NULL)
	{
	    if (extra ==  NULL)
	    {
	    	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
		extra->af_node_type = GORE_NODE;
	    }
	    n->extra = (POINTER)extra;
	    n->size_of_extra = sizeof(AF_NODE_EXTRA);
	}
}	/* end set_gore_node */

