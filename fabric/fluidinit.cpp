#include "fabric.h"

static void promptForDirichletBdryState(FILE*,Front*,HYPER_SURF**,int);
static void flowThroughBoundaryState3d( double*,HYPER_SURF*,Front*,POINTER,POINTER); 

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

static void promptForDirichletBdryState(
	FILE *infile,
	Front *front,
	HYPER_SURF **hs,
        int i_hs)
{
	static STATE *state;
	char s[100];
	POINTER func_params;
	int dim = FT_Dimension();
	int k;

	CursorAfterString(infile,"Enter type of Dirichlet boundary:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	switch (s[0])
	{
	case 'c':			// Constant state
	case 'C':
	    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
	    CursorAfterString(infile,"Enter velocity:");
	    for (k = 0; k < dim; ++k)
	    {
		fscanf(infile,"%lf",&state->vel[k]);
		(void) printf("%f ",state->vel[k]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter pressure:");
	    fscanf(infile,"%lf",&state->pres);
	    (void) printf("%f\n",state->pres);

            //state->phi = getPhiFromPres(front,state->pres);
	    
        FT_InsertDirichletBoundary(front,NULL,NULL,
			NULL,(POINTER)state,*hs,i_hs);
	    break;
	default:
	    (void) printf("Unknown Dirichlet boundary!\n");
	    clean_up(ERROR);
	}
}	/* end promptForDirichletBdryState */

EXPORT void restart_set_dirichlet_bdry_function(Front *front)
{
	INTERFACE *intfc = front->interf;
	BOUNDARY_STATE  *bstate;
	const char *s;
	for (int i = 0; i < num_bstates(intfc); ++i)
	{
	    bstate = bstate_list(intfc)[i];
	    if (bstate == NULL) continue;
	    s = bstate->_boundary_state_function_name;
	    if (s == NULL) continue;
	}
}	/* end restart_set_dirichlet_bdry_function */

