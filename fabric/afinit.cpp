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

static void setInitialIntfcAF3d(Front*,LEVEL_FUNC_PACK*,char*);

void setInitialIntfcAF(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	FILE *infile = fopen(inname,"r");
	char string[100];

	level_func_pack->wave_type = ELASTIC_BOUNDARY;
	
    F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
	Fparams->m_comp1 = SOLID_COMP;
    Fparams->m_comp2 = LIQUID_COMP2;

    if (CursorAfterStringOpt(infile,
            "Entering yes to set wave type to FIRST_PHYSICS_WAVE_TYPE: "))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
		level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
	}
	fclose(infile);

    assert(front->rect_grid->dim == 3);
    return setInitialIntfcAF3d(front,level_func_pack,inname);
}	/* end setInitialIntfcAF */

static void setInitialIntfcAF3d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	
    //not even used in this function
        //F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
	
    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int num_canopy;

	level_func_pack->set_3d_bdry = YES;
	
	level_func_pack->neg_component = LIQUID_COMP2;
    level_func_pack->pos_component = LIQUID_COMP2;	
	
    level_func_pack->func_params = NULL;
    level_func_pack->func = NULL;
	level_func_pack->attach_string = NO;		// default

	af_params->is_parachute_system = NO;
	af_params->cut_vent = NO;
	af_params->num_opt_round = 20;
    af_params->spring_model = MODEL1;	// default
	af_params->attach_gores = NO;		// default
	af_params->use_gpu = NO;		    // default
	af_params->gore_len_fac = 1.0;		// default
	
    CursorAfterString(infile,"Enter number of canopy surfaces:");
	fscanf(infile,"%d",&num_canopy);
	(void) printf("%d\n",num_canopy);
	level_func_pack->num_mono_hs = num_canopy;

	(void) printf("Choices of initial surface are:\n");
	(void) printf("\tEllipsoid (E)\n");
	(void) printf("\tPlane (P)\n");
	(void) printf("\tT-10 (T)\n");
	(void) printf("\tNone (N)\n");
	
    CursorAfterString(infile,"Enter initial canopy surface type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'E':
	case 'e':
	    initEllipticSurf(infile,front,level_func_pack);
	    break;
	case 'T':
	case 't':
	    initParabolicSurf(infile,front,level_func_pack);
	    break;
	case 'P':
	case 'p':
	    initPlaneSurf(infile,front,level_func_pack);
	    break;
	case 'A':
    case 'a':
        initAirbag(infile,front,level_func_pack);
	    break;
	case 'N':
	case 'n':
	    break;
	}
	af_params->pert_params.pert_type = NO_PERT;
	(void) printf("Available perturbation types are:\n");
	(void) printf("\tNo perturbation (N)\n");
	(void) printf("\tParallel random perturbation (P)\n");
	(void) printf("\tOrthogonal random perturbation (O)\n");
	(void) printf("\tRadial perturbation (R)\n");
	(void) printf("\tLinear perturbation (L)\n");
	(void) printf("\tSine perturbation (S)\n");
	(void) printf("\tDefault is no perturbation\n");
	if (CursorAfterStringOpt(infile,
	    "Entering perturbation type: "))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
	    {
	    case 'n':
	    case 'N':
		break;
	    case 'p':
	    case 'P':
		af_params->pert_params.pert_type = PARALLEL_RAND_PERT;
		break;
	    case 'o':
	    case 'O':
		af_params->pert_params.pert_type = ORTHOGONAL_RAND_PERT;
		break;
	    case 'r':
	    case 'R':
		af_params->pert_params.pert_type = RADIAL_PERT;
	    	CursorAfterString(infile,"Enter perturbation center:");
	        fscanf(infile,"%lf %lf",&af_params->pert_params.cen[0],
				&af_params->pert_params.cen[1]);
		(void) printf("%f %f\n",af_params->pert_params.cen[0],
				af_params->pert_params.cen[1]);
	    	CursorAfterString(infile,"Enter perturbation radius:");
	        fscanf(infile,"%lf",&af_params->pert_params.pert_radius);
		(void) printf("%f\n",af_params->pert_params.pert_radius);
	    	CursorAfterString(infile,"Enter perturbation amplitude:");
	        fscanf(infile,"%lf",&af_params->pert_params.pert_amp);
		(void) printf("%f\n",af_params->pert_params.pert_amp);
		break;
	    case 'l':
	    case 'L':
		af_params->pert_params.pert_type = LINEAR_PERT;
	    	CursorAfterString(infile,"Enter perturbation direction:");
	        fscanf(infile,"%d",&af_params->pert_params.dir);
		(void) printf("%d\n",af_params->pert_params.dir);
	    	CursorAfterString(infile,"Enter perturbation center:");
	        fscanf(infile,"%lf",&af_params->pert_params.x0);
		(void) printf("%f\n",af_params->pert_params.x0);
	    	CursorAfterString(infile,"Enter perturbation lower end:");
	        fscanf(infile,"%lf",&af_params->pert_params.xl);
		(void) printf("%f\n",af_params->pert_params.xl);
	    	CursorAfterString(infile,"Enter perturbation upper end:");
	        fscanf(infile,"%lf",&af_params->pert_params.xu);
		(void) printf("%f\n",af_params->pert_params.xu);
	    	CursorAfterString(infile,"Enter perturbation amplitude:");
	        fscanf(infile,"%lf",&af_params->pert_params.pert_amp);
		(void) printf("%f\n",af_params->pert_params.pert_amp);
		break;
	    case 's':
	    case 'S':
		af_params->pert_params.pert_type = SINE_PERT;
		break;
	    }
	}

	if (CursorAfterStringOpt(infile,
            "Entering type of spring model: "))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        switch (string[0])
        {
        case '1':
            af_params->spring_model = MODEL1;
            break;
        case '2':
            af_params->spring_model = MODEL2;
            break;
        case '3':
            af_params->spring_model = MODEL3;
            break;
        default:
            break;
        }
    }

	if (CursorAfterStringOpt(infile,
            "Entering number of canopy optimization rounds: "))
    {
        fscanf(infile,"%d",&af_params->num_opt_round);
        (void) printf("%d\n",af_params->num_opt_round);
    }
}	/* end setInitialIntfcAF3d */

EXPORT void read_Fparams(
	char *inname,
	F_PARAMS *Fparams)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int i,dim = Fparams->dim;

	/* defaults numerical schemes */
	Fparams->num_scheme.projc_method = SIMPLE;
	Fparams->num_scheme.advec_method = WENO;
	Fparams->num_scheme.ellip_method = SIMPLE_ELLIP;

	CursorAfterString(infile,"Enter projection type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'S':
	case 's':
	    Fparams->num_scheme.projc_method = SIMPLE;
	    break;
	case 'B':
	case 'b':
	    Fparams->num_scheme.projc_method = BELL_COLELLA;
	    break;
	case 'K':
	case 'k':
	    Fparams->num_scheme.projc_method = KIM_MOIN;
	    break;
	case 'P':
	case 'p':
	    Fparams->num_scheme.projc_method = PEROT_BOTELLA;
	}
	assert(Fparams->num_scheme.projc_method != ERROR_PROJC_SCHEME);
	(void) printf("The default advection order is WENO-Runge-Kutta 4\n");
	Fparams->adv_order = 4;
	if (CursorAfterStringOpt(infile,"Enter advection order:"))
	{
	    fscanf(infile,"%d",&Fparams->adv_order);
	    (void) printf("%d\n",Fparams->adv_order);
	}

	(void) printf("Available elliptic methods are:\n");
	(void) printf("\tSimple elliptic (S)\n");
	(void) printf("\tDouble elliptic (DB)\n");
	if (CursorAfterStringOpt(infile,"Enter elliptic method:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    switch (string[0])
	    {
	    case 'S':
	    case 's':
	    	Fparams->num_scheme.ellip_method = SIMPLE_ELLIP;
	    	break;
	    case 'd':
	    case 'D':
            Fparams->num_scheme.ellip_method = DOUBLE_ELLIP;
	    	break;
        default:
            printf("Elliptic Method Not Implemented\n");
            clean_up(1);
	    }
	}

	for (i = 0; i < dim; ++i) Fparams->U_ambient[i] = 0.0;
        if (CursorAfterStringOpt(infile,"Enter fluid ambient velocity:"))
        {
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf ",&Fparams->U_ambient[i]);
                (void) printf("%f ",Fparams->U_ambient[i]);
            }
            (void) printf("\n");
        }
	Fparams->ub_speed = HUGE;
        if (CursorAfterStringOpt(infile,"Enter upper bound for speed:"))
	{
            fscanf(infile,"%lf ",&Fparams->ub_speed);
            (void) printf("%f\n",Fparams->ub_speed);
	}
	Fparams->total_div_cancellation = NO;
        if (CursorAfterStringOpt(infile,	
		"Enter yes to use total divergence cancellation:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    	Fparams->total_div_cancellation = YES;
	}
        if (CursorAfterStringOpt(infile,
		"Enter density and viscosity of the fluid:"))
        {
            fscanf(infile,"%lf %lf",&Fparams->rho2,&Fparams->mu2);
            (void) printf("%f %f\n",Fparams->rho2,Fparams->mu2);
	}
	Fparams->use_eddy_visc = NO;
        if (CursorAfterStringOpt(infile,
		"Enter yes to use eddy viscosity:"))
        {
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	Fparams->use_eddy_visc = YES;
		(void) printf("Available turbulence models are:\n");
		(void) printf("\tBaldwin-Lomax (B)\n");
		(void) printf("\tMoin (M)\n");
        	CursorAfterString(infile,"Enter turbulence model:");
	    	fscanf(infile,"%s",string);
	    	(void) printf("%s\n",string);
		switch (string[0])
		{
		case 'b':
		case 'B':
		    Fparams->eddy_visc_model = BALDWIN_LOMAX;
        	    CursorAfterString(infile,
			"Enter maximum distance for eddy viscosity:");
            	    fscanf(infile,"%lf",&Fparams->ymax);
            	    (void) printf("%f\n",Fparams->ymax);
		    break;
		case 'm':
		case 'M':
		    Fparams->eddy_visc_model = MOIN;
		    break;
		case 'S':
		case 's':
		    Fparams->eddy_visc_model = SMAGORINSKY;
                    break;
		default:
		    (void) printf("Unknown eddy viscosity model!\n");
		    clean_up(ERROR);
		}
	    }
	}
        if (CursorAfterStringOpt(infile,"Enter gravity:"))
        {
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf ",&Fparams->gravity[i]);
                (void) printf("%f ",Fparams->gravity[i]);
            }
            (void) printf("\n");
        }
	Fparams->scalar_field = NO;
        if (CursorAfterStringOpt(infile,"Enter yes to consider scalar field:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                Fparams->scalar_field = YES;
        }
	Fparams->min_speed = 0.0;
        if (CursorAfterStringOpt(infile,
			"Enter minimum speed to limit time step:"))
        {
            fscanf(infile,"%lf ",&Fparams->min_speed);
            (void) printf("%f ",Fparams->min_speed);
            (void) printf("\n");
        }
	fclose(infile);
}	/* end read_Fparams */

EXPORT void read_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100];
	int i,j,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	INTERFACE *intfc = front->interf;
	HYPER_SURF *hs;
	int i_hs = 0;

	(void) printf("Available type of Dirichlet boundary include: \n");
	(void) printf("\tConstant state (C)\n");
	(void) printf("\tFlow through (F)\n");
	(void) printf("\tTime dependent (T)\n");
	(void) printf("\tSplit state (S)\n");
	for (i = 0; i < dim; ++i)
	for (j = 0; j < 2; ++j)
	{
	    if (rect_boundary_type(intfc,i,j) == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
	        if (rect_boundary_type(front->interf,i,j) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,
					DIRICHLET_BOUNDARY,i,j);
		if (j == 0)
		    sprintf(msg,"For lower boundary in %d-th dimension",i);
		else
		    sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		promptForDirichletBdryState(infile,front,&hs,i_hs);
		i_hs++;
	    }
	    else if (rect_boundary_type(intfc,i,j) == MIXED_TYPE_BOUNDARY)
            {
		HYPER_SURF **hss;
		int k,nhs;
                hss = FT_MixedBoundaryHypSurfs(intfc,i,j,DIRICHLET_BOUNDARY,
                                        &nhs);
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
                        promptForDirichletBdryState(infile,front,hss+k,i_hs);
                        i_hs++;
                    }
                }
            }
	}
	fclose(infile);
}	/* end read_iF_dirichlet_bdry_data */

EXPORT void restart_set_dirichlet_bdry_function(Front *front)
{
	INTERFACE *intfc = front->interf;
	int i;
	BOUNDARY_STATE  *bstate;
	const char *s;
	for (i = 0; i < num_bstates(intfc); ++i)
	{
	    bstate = bstate_list(intfc)[i];
	    if (bstate == NULL) continue;
	    s = bstate->_boundary_state_function_name;
	    if (s == NULL) continue;
	    if (strcmp(s,"flowThroughBoundaryState") == 0)
            	bstate->_boundary_state_function = iF_flowThroughBoundaryState;
	}
}	/* end restart_set_dirichlet_bdry_function */

