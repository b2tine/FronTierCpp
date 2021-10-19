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


static void getFabricState(STATE*,EQN_PARAMS*,double*,COMPONENT);

//TODO: Move to cFinit.cpp and remove cFphys.cpp entirely
void read_cFluid_params(char* inname, EQN_PARAMS* eqn_params)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int i,dim = eqn_params->dim;

	eqn_params->prob_type = ERROR_TYPE;
	eqn_params->tracked = NO;		// Default
	CursorAfterString(infile,"Enter problem type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	
    eqn_params->prob_type = FABRIC; //default
    if (string[0] == 'F' || string[0] == 'f')
    {
        if (string[1] == 'A' || string[1] == 'a')
            eqn_params->prob_type = FABRIC;
    }

    /*
    if (string[0] == 'T' || string[0] == 't')
	{
	    if (string[10] == 'B' || string[10] == 'b')
	    	eqn_params->prob_type = TWO_FLUID_BUBBLE;
	    else if (string[10] == 'R' || string[10] == 'r')
	    {
		if (string[11] == 'T' || string[11] == 't')
	    	    eqn_params->prob_type = TWO_FLUID_RT;
		else if (string[11] == 'M' || string[11] == 'm')
	    	    eqn_params->prob_type = TWO_FLUID_RM;
	    }
	} 
	else if (string[0] == 'F' || string[0] == 'f')
	{
            if (string[12] == 'C' || string[12] =='c')
            {
                if (string[13] == 'I' || string[13] =='i')
                    eqn_params->prob_type = FLUID_SOLID_CIRCLE;
                else if (string[13] == 'Y' || string[13] =='y')
                    eqn_params->prob_type = FLUID_SOLID_CYLINDER;
            }
            else if (string[12] == 'R' || string[12] =='r')
                eqn_params->prob_type = FLUID_SOLID_RECT;
            else if (string[12] == 'T' || string[12] =='t')
                eqn_params->prob_type = FLUID_SOLID_TRIANGLE;
        }
	else if (string[0] == 'B' || string[0] == 'b')
	    eqn_params->prob_type = BUBBLE_SURFACE;
	else if (string[0] == 'I' || string[0] == 'i')
	    eqn_params->prob_type = IMPLOSION;
	else if (string[0] == 'P' || string[0] == 'p')
	    eqn_params->prob_type = PROJECTILE;
	else if (string[0] == 'R' || string[0] == 'r')
	    eqn_params->prob_type = RIEMANN_PROB;
	else if (string[0] == 'M' || string[0] == 'm')
	    eqn_params->prob_type = MT_FUSION;
	else if (string[0] == 'O' || string[0] == 'o')
	{
	    if (string[1] == 'B' || string[1] == 'b')
		eqn_params->prob_type = OBLIQUE_SHOCK_REFLECT;
	    else if (string[5] == 'S' || string[5] == 's')
	    	eqn_params->prob_type = ONED_SSINE;
	    else if (string[5] == 'B' || string[5] == 'b')
	    	eqn_params->prob_type = ONED_BLAST;
	    else if (string[5] == 'A' || string[5] == 'a')
	    	eqn_params->prob_type = ONED_ASINE;
	}
    */
        printf("Available numerical schemes are:\n");
        printf("\tTVD_1st_order\n");
        printf("\tTVD_2nd_order\n");
        printf("\tTVD_4th_order\n");
        printf("\tWENO_1st_order\n");
        printf("\tWENO_2nd_order\n");
        printf("\tWENO_4th_order\n");
	CursorAfterString(infile,"Enter numerical scheme for interior solver:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'T':
	case 't':
	    switch (string[4])
	    {
	    case '1':
		eqn_params->num_scheme = TVD_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = TVD_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = TVD_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    break;
	case 'W':
	case 'w':
	    switch (string[5])
	    {
	    case '1':
		eqn_params->num_scheme = WENO_FIRST_ORDER;
		break;
	    case '2':
		eqn_params->num_scheme = WENO_SECOND_ORDER;
		break;
	    case '4':
		eqn_params->num_scheme = WENO_FOURTH_ORDER;
		break;
	    default:
		printf("Numerical scheme %s not implemented!\n",string);
		clean_up(ERROR);
	    }
	    eqn_params->articomp = NO;
	    if (CursorAfterStringOpt(infile,
		"Enter yes to use artificial compression:"))
	    {
	    	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
	    	if (string[0] == 'y' || string[0] == 'Y')
	            eqn_params->articomp = YES;
	    }
	    break;
	default:
	    printf("Numerical scheme %s not implemented!\n",string);
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Enter order of point propagator:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case '1':
	    eqn_params->point_prop_scheme = FIRST_ORDER;
	    break;
	case '2':
	    eqn_params->point_prop_scheme = SECOND_ORDER;
	    break;
	case '4':
	    eqn_params->point_prop_scheme = FOURTH_ORDER;
	    break;
	default:
	    printf("Point propagator order %s not implemented!\n",string);
	    clean_up(ERROR);
	}

    eqn_params->use_eddy_viscosity = false;
    if (CursorAfterStringOpt(infile,"Enter yes to use eddy viscosity:"))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
            eqn_params->use_eddy_viscosity = true;
    }

	eqn_params->use_base_soln = NO;
	if (CursorAfterStringOpt(infile,
		"Enter yes for comparison with base data:"))
	{
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
            	eqn_params->use_base_soln = YES;
	}

	if (eqn_params->use_base_soln == YES)
    {
	    CursorAfterString(infile,"Enter base directory name:");
            fscanf(infile,"%s",eqn_params->base_dir_name);
            (void) printf("%s\n",eqn_params->base_dir_name);
            CursorAfterString(infile,"Enter number of comparing steps:");
            fscanf(infile,"%d",&eqn_params->num_step);
            (void) printf("%d\n",eqn_params->num_step);
            FT_VectorMemoryAlloc((POINTER*)&eqn_params->steps,
                                eqn_params->num_step,sizeof(int));
            for (i = 0; i < eqn_params->num_step; ++i)
            {
                sprintf(string,"Enter index of step %d:",i+1);
                CursorAfterString(infile,string);
                fscanf(infile,"%d",&eqn_params->steps[i]);
                (void) printf("%d\n",eqn_params->steps[i]);
            }
            FT_ScalarMemoryAlloc((POINTER*)&eqn_params->f_basic,
                                sizeof(F_BASIC_DATA));
            eqn_params->f_basic->dim = dim;
	}

    if (eqn_params->prob_type == ERROR_TYPE)
    {
        printf("eqn_params->prob_type == ERROR_TYPE\n");
        LOC(); clean_up(ERROR);
    }

	fclose(infile);

	if (eqn_params->use_base_soln == YES)
	    FT_ReadComparisonDomain(inname,eqn_params->f_basic);
}	/* end read_cFluid_params */

//TODO: should move to another file -- 
void set_cFluid_params(char* inname, EQN_PARAMS* eqn_params)
{
	switch (eqn_params->prob_type)
	{
    case FABRIC:
        setFabricParams(inname,eqn_params);
        break;
    /*
	case TWO_FLUID_RT:
	    setRayleiTaylorParams(inname,eqn_params);
	    break;
	case TWO_FLUID_RM:
	case TWO_FLUID_RM_RAND:
	    setRichtmyerMeshkovParams(inname,eqn_params);
	    break;
	case TWO_FLUID_BUBBLE:
	    setBubbleParams(inname,eqn_params);
	    break;
	case IMPLOSION:
	    setImplosionParams(inname,eqn_params);
	    break;
	case MT_FUSION:
	    setMTFusionParams(inname,eqn_params);
	    break;
	case PROJECTILE:
	case FLUID_SOLID_CIRCLE:
	case FLUID_SOLID_RECT:
	case FLUID_SOLID_TRIANGLE:
	case FLUID_SOLID_CYLINDER:
	    setProjectileParams(inname,eqn_params);
	    break;
	case RIEMANN_PROB:
	    setRiemProbParams(inname,eqn_params);
	    break;
	case ONED_BLAST:
	case ONED_SSINE:
	case ONED_ASINE:
	    setOnedParams(inname,eqn_params);
	    break;
	case OBLIQUE_SHOCK_REFLECT:
	    setRichtmyerMeshkovParams(inname,eqn_params);
	    break;
        */
	default:
	    printf("\nIn set_cFluid_Params(), only FABRIC problem type supported currently!\n");
	        //printf("In set_cFluid_Params(), unknown problem type!\n");
	    LOC(); clean_up(ERROR);
	}
}	/* end set_cFluid_params */

//TODO: Change to setChannelFlowParams() and match to cFluid dir version
void setFabricParams(char* inname, EQN_PARAMS* eqn_params)
{
    int dim = eqn_params->dim;
	FILE *infile = fopen(inname,"r");
	double pinf,einf,gamma;
	char s[100], str[256];

	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
	             GAS_COMP1);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP1]).gamma = gamma;
	(eqn_params->eos[GAS_COMP1]).pinf = pinf;
	(eqn_params->eos[GAS_COMP1]).einf = einf;
	
	sprintf(str, "Enter gamma, pinf, einf of the fluid with comp %d:", 
		     GAS_COMP2);
	CursorAfterString(infile,str);
	fscanf(infile,"%lf %lf %lf",&gamma,&pinf,&einf);
	(void) printf("%f %f %f\n",gamma,pinf,einf);
	(eqn_params->eos[GAS_COMP2]).gamma = gamma;
	(eqn_params->eos[GAS_COMP2]).pinf = pinf;
	(eqn_params->eos[GAS_COMP2]).einf = einf;

	CursorAfterString(infile,"Enter density and pressure of ambient air1:");
	fscanf(infile,"%lf %lf",&eqn_params->rho1,&eqn_params->p1);
	(void) printf("%f %f\n",eqn_params->rho1,eqn_params->p1);

	CursorAfterString(infile,"Enter density and pressure of ambient air2:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->p2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->p2);

	CursorAfterString(infile,"Enter density and viscosity of the fluid:");
	fscanf(infile,"%lf %lf",&eqn_params->rho2,&eqn_params->mu2);
	(void) printf("%f %f\n",eqn_params->rho2,eqn_params->mu2);

    eqn_params->rho1 = eqn_params->rho2;
    eqn_params->mu1 = eqn_params->mu2;
    eqn_params->p1 = eqn_params->p2;

	CursorAfterString(infile,"Enter gravity:");
	for (int i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&eqn_params->gravity[i]);
	    (void) printf("%f",eqn_params->gravity[i]);
	}
    (void) printf("\n");

	CursorAfterString(infile,"Type yes to track the interface:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    eqn_params->tracked = YES;
	else
	    eqn_params->tracked = NO;
	fclose(infile);
}	/* end setFabricParams */

void G_CARTESIAN::setInitialStates()
{
	switch (eqn_params->prob_type)
	{
    case FABRIC:
        initFabricStates();
        break;
        /*
	case TWO_FLUID_RT:
	    initRayleiTaylorStates();
	    break;
	case TWO_FLUID_RM:
	case TWO_FLUID_RM_RAND:
	    initRichtmyerMeshkovStates();
	    break;
	case TWO_FLUID_BUBBLE:
	    initBubbleStates();
	    break;
	case IMPLOSION:
	    initImplosionStates();
	    break;
	case MT_FUSION:
	    initMTFusionStates();
	    break;
	case PROJECTILE:
        case FLUID_SOLID_CIRCLE:
        case FLUID_SOLID_RECT:
        case FLUID_SOLID_TRIANGLE:
        case FLUID_SOLID_CYLINDER:
	    initProjectileStates();
	    break;
	case RIEMANN_PROB:
	    initRiemProbStates();
	    break;
	case ONED_BLAST:
	    initBlastWaveStates();
	    break;
	case ONED_SSINE:
	    initShockSineWaveStates();
	    break;
	case ONED_ASINE:
	    initAccuracySineWaveStates();
	    break;
	case OBLIQUE_SHOCK_REFLECT:
	    initRichtmyerMeshkovStates();
	    break;
        */
	default:
	    (void) printf("In setInitialStates(), case not implemented!\n");
	    clean_up(ERROR);
	}

	copyMeshStates();
}	/* end setInitialStates */

//TODO: rename to initChannelFlowStates() like in cFluid dir
void G_CARTESIAN::initFabricStates()
{
	int i,j,k,l,index;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double coords[MAXD];
	COMPONENT comp;
	STATE *sl,*sr,state;
    POINT *p;
    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
	INTERFACE *intfc = front->interf;

    double *mu = field.mu;
	double *dens = field.dens;
	double *engy = field.engy;
	double *pres = field.pres;
	double **momn = field.momn;
	double **vel = field.vel;

    m_dens[0] = eqn_params->rho1;
    m_dens[1] = eqn_params->rho2;

    m_mu[0] = eqn_params->mu1;
    m_mu[1] = eqn_params->mu2;

    m_comp[0] = SOLID_COMP;
    m_comp[1] = GAS_COMP2;

    next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
        //TODO: check if wave_type(hs) == ELASTIC_BOUNDARY ??
        FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
        getFabricState(sl,eqn_params,Coords(p),negative_component(hs));
        getFabricState(sr,eqn_params,Coords(p),positive_component(hs));
	}

	FT_MakeGridIntfc(front);
	setDomain();

	switch (dim)
	{
	case 2:
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
            index = d_index2d(i,j,top_gmax);
            comp = top_comp[index];
            getRectangleCenter(index,coords);
            
            getFabricState(&state,eqn_params,coords,comp);
            
            mu[index] = state.mu;
            dens[index] = state.dens;
            pres[index] = state.pres;
            engy[index] = state.engy;
            for (l = 0; l < dim; ++l)
            {
                momn[l][index] = state.momn[l];
                vel[l][index] = state.vel[l];
            }
	    }
        break;
	case 3:
	    for (k = imin[2]; k <= imax[2]; ++k)
	    for (j = imin[1]; j <= imax[1]; ++j)
	    for (i = imin[0]; i <= imax[0]; ++i)
	    {
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            getRectangleCenter(index,coords);
            
            getFabricState(&state,eqn_params,coords,comp);
            
            mu[index] = state.mu;
            dens[index] = state.dens;
            pres[index] = state.pres;
            engy[index] = state.engy;
            for (l = 0; l < dim; ++l)
            {
                momn[l][index] = state.momn[l];
                vel[l][index] = state.vel[l];
            }
        }
        break;
	}
	scatMeshStates();
}	/* end initFabricStates */

//TODO: This is same as getAmbientState() -- reuse from cFinit.cpp
static void getFabricState(
	STATE *state,
	EQN_PARAMS *eqn_params,
	double *coords,
	COMPONENT comp)
{
	EOS_PARAMS	*eos;
    double mu1 = eqn_params->mu1;
    double mu2 = eqn_params->mu2;
	double rho1 = eqn_params->rho1;
	double rho2 = eqn_params->rho2;
	double p1 = eqn_params->p1;
	double p2 = eqn_params->p2;
	double *v1 = eqn_params->v1;
	double *v2 = eqn_params->v2;
	int i,dim;
 
	if (debugging("ambient"))
	    printf("Entering getFabricState(), coords = %f %f\n",
				coords[0],coords[1]);
	dim = eqn_params->dim;

	for (i = 0; i < dim; ++i)
    {
	    state->vel[i] = 0.0;
        state->momn[i] = 0.0;
    }
    state->dim = dim;
	eos = &(eqn_params->eos[comp]);
	state->eos = eos;
	
	switch (comp)
	{
	case GAS_COMP1:
        state->mu = mu1;
	    state->dens = rho1;
	    state->pres = p1;
	    for (i = 0; i < dim; ++i)
	    {
	    	state->vel[i] = v1[i];
	    	state->momn[i] = rho1*v1[i];
	    }
	    state->engy = EosEnergy(state);
	    break;
	case GAS_COMP2:
        state->mu = mu2;
	    state->dens = rho2;
	    state->pres = p2;
	    for (i = 0; i < dim; ++i)
	    {
	    	state->vel[i] = v2[i];
	    	state->momn[i] = rho2*v2[i];
	    }
	    state->engy = EosEnergy(state);
	    break;
	case SOLID_COMP:
        state->mu = 0.0;
	    state->dens = 0.0;
	    state->pres = 0.0;
        for (i = 0; i < dim; ++i)
        {
            state->vel[i] = 0.0;
            state->momn[i] = 0.0;
        }
	    state->engy = 0.0;
	    break;
	default:
	    printf("ERROR: Unknown component %d in getFabricState()!\n",
				comp);
	    clean_up(ERROR);
	}
	if (debugging("airfoil_state"))
	    (void) printf("Leaving getFabricState(): state = %d %f %f %f\n",
			comp,state->dens,state->pres,state->engy);
}	/* end getAmbientState */


