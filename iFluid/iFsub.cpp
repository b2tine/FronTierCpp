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

	/*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void rgbody_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void iF_flowThroughBoundaryState2d(double*,HYPER_SURF*,Front*,POINTER,
        		POINTER);
static void iF_flowThroughBoundaryState3d(double*,HYPER_SURF*,Front*,POINTER,
        		POINTER);
static void iF_splitBoundaryState(double*,HYPER_SURF*,Front*,POINTER,POINTER);
static void iF_parabolicBoundaryState(double*,HYPER_SURF*,Front*,
						POINTER,POINTER);
static void get_time_dependent_params(int,FILE*,POINTER*);
static void get_split_state_params(Front*,FILE*,POINTER*);
static void get_parabolic_state_params(Front*,FILE*,POINTER*);
static void addToEnergyFlux(RECT_GRID*,HYPER_SURF*,double*,double*,int,int,
			boolean);
static void promptForDirichletBdryState(FILE*,Front*,HYPER_SURF**,int);

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

static void ifluid_compute_force_and_torque2d(Front*,HYPER_SURF*,double,
                        double*,double*);
static void ifluid_compute_force_and_torque3d(Front*,HYPER_SURF*,double,
                        double*,double*);
static boolean force_on_hse(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,double*,
                                        double*,double*,boolean);
static boolean force_on_hse2d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static boolean force_on_hse3d(HYPER_SURF_ELEMENT*,HYPER_SURF*,RECT_GRID*,
                                        double*,double*,double*,boolean);
static double intrp_between(double,double,double,double,double);
static void setStateViscosity(IF_PARAMS*,STATE*,int);
//static void prompt_for_velocity_func(int,char*,RG_PARAMS*);
//static void sine_vel_func(Front*,POINTER,double*,double*);
static void pipe_end_func(Front*,POINTER,int*,COMPONENT,
                                int,int,int*,Locstate);
static boolean coords_in_subdomain(double*,RECT_GRID*);
static void xgraphAtOldNode(const char*,NODE*,O_CURVE,O_CURVE,O_CURVE);
static int modify_contact_node(NODE*,NODE*,O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
                              O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
                              O_CURVE*,POINT*,BOND*,BOND*,ANGLE_DIRECTION,
                              double,double,RPROBLEM**,Front*,POINTER,
                              double,double*,NODE_FLAG);


extern int next_index_in_dir(
        int* icoords,
        GRID_DIRECTION dir,
        int dim,
        int* top_gmax)
{
	int icrds[MAXD];
	for (int i = 0; i < dim; ++i)
	    icrds[i] = icoords[i];

    switch (dir)
    {
        case WEST:
            icrds[0] -= 1;
            break;
        case EAST:
            icrds[0] += 1;
            break;
        case SOUTH:
            icrds[1] -= 1;
            break;
        case NORTH:
            icrds[1] += 1;
            break;
        case LOWER:
            icrds[2] -= 1;
            break;
        case UPPER:
            icrds[2] += 1;
    }

    int index = d_index(icrds,top_gmax,dim);
	return index;
}

extern double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

extern double getStatePhi(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->phi;
}	/* end getStatePhi */

extern double getStateGradPhiX(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->grad_phi[0];
}	/* end getStateGradPhiX */

extern double getStateGradPhiY(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->grad_phi[1];
}	/* end getStateGradPhiY */

extern double getStateGradPhiZ(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->grad_phi[2];
}	/* end getStateGradPhiZ */

extern double getStateQ(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->q;
}	/* end getStateQ */

extern double getStateVort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort;
}	/* end getStateVort */

extern double getStateXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[0];
}	/* end getStateXvel */

extern double getStateYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[1];
}	/* end getStateYvel */

extern double getStateZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[2];
}	/* end getStateZvel */

extern double getStateOldXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel_old[0];
}	/* end getStateXvel */

extern double getStateOldYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel_old[1];
}	/* end getStateYvel */

extern double getStateOldZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel_old[2];
}	/* end getStateZvel */

extern double getStateXimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[0];
}	/* end getStateXimp */

extern double getStateYimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[1];
}	/* end getStateYimp */

extern double getStateZimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[2];
}	/* end getStateZimp */

extern double getStateMu(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->mu;
}	/* end getStateMu */

extern double getStateDens(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->dens;
}	/* end getStateDens */

extern double getStateTemp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->temperature;
}	/* end getStateMTemp */


extern void read_iF_dirichlet_bdry_data(
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
                hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,i,j);
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
            hss = FT_MixedBoundaryHypSurfs(intfc,i,j,DIRICHLET_BOUNDARY,&nhs);
            printf("Number of Dirichlet boundaries on dir %d side %d: %d\n",i,j,nhs);
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

extern void restart_set_dirichlet_bdry_function(Front *front)
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

static void iF_splitBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	SPLIT_STATE_PARAMS *sparams = (SPLIT_STATE_PARAMS*)params;
	STATE *iF_state = (STATE*)state;
	int dir = sparams->dir;

	if (p0[dir] < sparams->split_coord)
	    *iF_state =  sparams->left_state;
	else
	    *iF_state =  sparams->right_state;
}	/* end iF_splitBoundaryState */

static void iF_parabolicBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
        PARABOLIC_STATE_PARAMS *pparams = (PARABOLIC_STATE_PARAMS*)params;
        STATE *iF_state = (STATE*)state;
        int i,j,i_nb;
        int dir = pparams->dir;
        int dim = FT_Dimension();
        double *L = front->rect_grid->L;
        double *U = front->rect_grid->U;
	double *h = front->rect_grid->h;

        *iF_state = pparams->state;
        for (i = 0; i < dim; i++)
        {
            iF_state->vel[i] = pparams->v_peak[i];
            for (j = 0; j < dim-1; j++)
            {
                i_nb = (dir+j+1)%dim;
		if (p0[i_nb] > 0.5 * (L[i_nb] + U[i_nb]) 
		    && p0[i_nb] < U[i_nb] - 0.6 * h[i_nb])
                    iF_state->vel[i]  *= -16*p0[i_nb]*p0[i_nb]+24*p0[i_nb]-8;
		else
		    iF_state->vel[i] = 0.0;
            }
        }
}       /* end iF_parabolicBoundaryState */

extern void iF_timeDependBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	TIME_DEPENDENT_PARAMS *td_params = (TIME_DEPENDENT_PARAMS*)params;
	STATE *iF_state = (STATE*)state;
	int i,dim = front->rect_grid->dim;
	double time = front->time;
	double *T = td_params->T;
	double omega = td_params->omega;
	double phase = td_params->phase;
	static int step = 0;

	switch (td_params->td_type)
	{
	case CONSTANT:
	    for (i = 0; i < dim; ++i)
	    	iF_state->vel[i] = td_params->v_base[i];
	    iF_state->pres = td_params->p_base;
	    break;
	case PULSE_FUNC:
	    if (time <= T[0])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_base[i];
	    	iF_state->pres = td_params->p_base;
	    }
	    else if (time > T[0] && time <= T[1])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_base[i] + 
				(time - T[0])/(T[1] - T[0])*
				(td_params->v_peak[i] - td_params->v_base[i]);
	    	iF_state->pres = td_params->p_base + 
				(time - T[0])/(T[1] - T[0])*
				(td_params->p_peak - td_params->p_base);
	    }
	    else if (time > T[1] && time <= T[2])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_peak[i];
	    	iF_state->pres = td_params->p_peak;
	    }
	    else if (time > T[2] && time <= T[3])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_peak[i] + 
				(time - T[2])/(T[3] - T[2])*
				(td_params->v_tail[i] - td_params->v_peak[i]);
	    	iF_state->pres = td_params->p_peak + 
				(time - T[2])/(T[3] - T[2])*
				(td_params->p_tail - td_params->p_peak);
	    }
	    if (time > T[3])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_tail[i];
	    	iF_state->pres = td_params->p_tail;
	    }
	    break;
	case SINE_FUNC:
	    for (i = 0; i < dim; ++i)
	    	iF_state->vel[i] = td_params->v_base[i] +
			td_params->v_amp[i]*sin(omega*time + phase);
	    iF_state->pres = td_params->p_base + td_params->p_amp*
				sin(omega*time + phase);
	    break;
	default:
	    (void) printf("In iF_timeDependBoundaryState(), unknown type!\n");
	    clean_up(ERROR);
	}
	if (debugging("time_depend_bdry"))
	{
	    if (step != front->step)
	    {
	    	printf("time = %f  vel = %f %f  p = %f\n",time,iF_state->vel[0],
				iF_state->vel[1],iF_state->pres);
	    	step = front->step;
	    }
	}
}	/* end iF_timeDependBoundaryState */

extern void iF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    return iF_flowThroughBoundaryState2d(p0,hs,front,params,state);
	case 3:
	    return iF_flowThroughBoundaryState3d(p0,hs,front,params,state);
	}
}	/* end iF_flowThroughBoundaryState */

//NEW VERSION
static void iF_flowThroughBoundaryState2d(
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

	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;

	double dir[MAXD];
	double u[3];		// velocity in the sweeping direction
	double v[3][MAXD];	// velocity in the orthogonal direction
	double vort[3];		// vorticity stencil
	double pres[3];		// pressure stencil
	double phi[3];		// phi stencil
	double f_u;		    // u flux in the sweeping direction
	double f_v[MAXD];	// v flux in the orthogonal direction
	double f_vort;		// vort flux
	double f_pres;		// pressure flux
	double f_phi;		// phi flux
	double dn,dt = front->dt;
	int i,j,dim = front->rect_grid->dim;

	STATE *oldst, *newst = (STATE*)state;
	STATE  **sts;
    POINTER sl,sr;
	
    int nrad = 3;
	
	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState2d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
    
	if (comp == negative_component(hs))  
	    oldst = (STATE*)sl;
	else 
	    oldst = (STATE*)sr;

    for (i = 0; i < dim; ++i)
    {
        newst->vel[i] = oldst->vel[i];
        newst->vel_old[i] = oldst->vel[i];
    }
	newst->vort = oldst->vort;
	newst->pres = oldst->pres;
	newst->phi = oldst->phi;
    

    //TODO: Computing phi at the boundary can not be done with this function!
    //
    //      We require the intermediate velocity, u*, of the current time step
    //      to evaluate phi at the flow-through boundary, but this function
    //      gets called before we compute u*.
    //
    //      Need to apply the boundary condition for phi after computeDiffusion()
    //      and before computeProjection().
    //
    //      It appears that we can compute it inside the function ellip.cpp : solve2d().
    //      If successful, the phi computations in this function should be removed.


    //Normal
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);

    /*
    printf("\niF_flowThroughBoundaryState2d() Normal Stencil:\n");
    fprint_general_vector(stdout,"at point",Coords(oldp),dim,"\n");
    fprint_general_vector(stdout,"normal",nsten->nor,dim,"\n");
    fprint_general_vector(stdout,"nsten->pts[0]",nsten->pts[0],dim,"\n");
    fprint_general_vector(stdout,"nsten->pts[1]",nsten->pts[1],dim,"\n");
    fprint_general_vector(stdout,"nsten->pts[2]",nsten->pts[2],dim,"\n\n");
    */

	for (j = 0; j < 3; ++j) u[j] = 0.0;

    for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(nsten->nor,front);
	
    if (debugging("flow_through"))
	{
	    (void) printf("Normal grid size = %f\n",dn);
	    (void) print_Nor_stencil(front,nsten);
	}

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

	    vort[j] = oldst->vort;
	    pres[j] = oldst->pres;
	    phi[j] = oldst->phi;
	}

    STATE s1;
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
        FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
                field->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
	    s1.vel[i] = vtmp;
	}

	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->vort,
                            getStateVort,&vort[2],&oldst->vort);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2],&oldst->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->phi,
                            getStatePhi,&phi[2],&oldst->phi);

    s1.vort = vort[2];
    s1.pres = pres[2];
    s1.phi = phi[2];

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
	f_phi = linear_flux(u[1],phi[0],phi[1],phi[2]);

	for (i = 0; i < dim; ++i)
    {
        newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
    }
	newst->vort -= dt/dn*f_vort;
	newst->pres -= dt/dn*f_pres;
    newst->phi -= dt/dn*f_phi;

    //Tangential
	tsten = FrontGetTanStencils(front,oldp,nrad);

    /*
    printf("\niF_flowThroughBoundaryState2d() Tangent Stencil:\n");
    for (j = 0; j < 3; ++j)
    {
        std::string printmsg = "tsten[0]->p[" + std::to_string(j-1) + "]";
        fprint_general_vector(stdout,printmsg.c_str(),Coords(tsten[0]->p[j-1]),dim,"\n");
    }
    printf("\n\n");
    */

    if (debugging("flow_through"))
	{
	    (void) printf("Ambient component: %d\n",comp);
	    (void) printf("Tangential grid size = %f\n",dn);
	    (void) print_Tan_stencil(front,tsten[0]);
	}

	for (j = 0; j < 3; ++j) u[j] = 0.0;
	
	for (i = 0; i < dim; ++i)
	    dir[i] = tsten[0]->dir[i];
	dn = FT_GridSizeInDir(dir,front);

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
	    phi[j] = sts[j-1]->phi;
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
    {
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
    }
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_phi = linear_flux(u[1],phi[0],phi[1],phi[2]);

	for (i = 0; i < dim; ++i)
    {
	    newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
    }
	newst->vort -= dt/dn*f_vort;
	newst->pres -= dt/dn*f_pres;
    newst->phi -= dt/dn*f_phi;
    
    //Since pressure is usually updated as
    //      
    //      p^{n+1/2} = q + phi^{n+1} - 0.5*mu*(Div_U);
    //
    //      set 
    //
    //      phi^{n+1} = p^{n+1/2} - q    (Div_U = 0 at the boundary) 

    newst->q = oldst->pres;
        //newst->q = 0.0;
    
    //TODO: correct?
    if (iFparams->num_scheme.projc_method == PMI ||
        iFparams->num_scheme.projc_method == PMII)
    {
        newst->phi = newst->pres - newst->q;
    }
    
    if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) printf("Vorticity: %f\n",newst->vort);
	    (void) printf("Pressure: %f\n",newst->pres);
	    (void) printf("Phi: %f\n",newst->phi);
	    (void) printf("q: %f\n",newst->q);
	}
}       /* end iF_flowThroughBoundaryState2d */

//NEW VERSION
static void iF_flowThroughBoundaryState3d(
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
	
    IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
	
    double dir[MAXD];
	double u[3];		// velocity in the sweeping direction
	double v[3][MAXD];	// velocity in the orthogonal direction
	double pres[3];		// pressure stencil
	double phi[3];		// phi stencil
	double f_u;		    // u flux in the sweeping direction
	double f_v[MAXD];	// v flux in the orthogonal direction
	double f_pres;		// pressure flux
	double f_phi;		// phi flux
	double dn;
    
    int dim = front->rect_grid->dim;
    double dt = front->dt;
	
    STATE *oldst, *newst = (STATE*)state;
    STATE  **sts;
	POINTER sl,sr;
	int i,j,k;
	
    int nrad = 3;

	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState3d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
	
    if (comp == negative_component(hs))  
	    oldst = (STATE*)sl;
	else 
	    oldst = (STATE*)sr;
    
    for (i = 0; i < dim; ++i)
    {
        newst->vel[i] = oldst->vel[i];
        newst->vel_old[i] = oldst->vel[i];
    }
    newst->pres = oldst->pres;
    newst->phi = oldst->phi;


    //Normal
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	
    for (j = 0; j < 3; ++j) u[j] = 0.0;
	
    for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);
    
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
	    phi[j] = oldst->phi;
	}

    STATE s1;
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
                field->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
	    s1.vel[i] = vtmp;
	}

	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2],&oldst->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->phi,
                            getStatePhi,&phi[2],&oldst->phi);

    s1.pres = pres[2];
    s1.phi = phi[2];
    
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
    f_phi = linear_flux(u[1],phi[0],phi[1],phi[2]);

    for (i = 0; i < dim; ++i)
    {
        newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
    }
    newst->pres -= dt/dn*f_pres;
    newst->phi -= dt/dn*f_phi;

    if (debugging("flow_through"))
	{
	    (void) print_Nor_stencil(front,nsten);
	    (void) printf("new velocity after normal prop: %f %f %f\n",
			newst->vel[0],newst->vel[1],newst->vel[2]);
	}
	
    
    //Tangential
	tsten = FrontGetTanStencils(front,oldp,nrad);
	if (tsten != nullptr)
    {
	    for (k = 0; k < dim-1; ++k)
        {
            for (j = 0; j < 3; ++j) u[j] = 0.0;

            if (comp == negative_component(hs))  
                sts = (STATE**)tsten[k]->leftst;
            else 
                sts = (STATE**)tsten[k]->rightst;

            for (i = 0; i < dim; ++i)
                dir[i] = tsten[k]->dir[i];
            dn = FT_GridSizeInDir(dir,front);

            if (debugging("flow_through"))
            {
                printf("Ambient component: %d\n",comp);
                printf("Tangential grid size = %f\n",dn);
                printf("For direction %d\n",k);
                print_Tan_stencil(front,tsten[k]);
            }

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
                phi[j] = sts[j-1]->phi;
            }

            f_u = burger_flux(u[0],u[1],u[2]);
            for (i = 0; i < dim; ++i)
            {
                f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
            }
            f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
            f_phi = linear_flux(u[1],phi[0],phi[1],phi[2]);

            for (i = 0; i < dim; ++i)
            {
                newst->vel[i] -= dt/dn*(f_u*dir[i] + f_v[i]);
            }
            
            newst->pres -= dt/dn*f_pres;
            newst->phi -= dt/dn*f_phi;
        }
    }

    newst->q = oldst->pres;
        //newst->q = 0.0;
     
    //TODO: correct?
    if (iFparams->num_scheme.projc_method == PMI ||
        iFparams->num_scheme.projc_method == PMII)
    {
        newst->phi = newst->pres - newst->q;
    }

    if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) printf("Pressure: %f\n",newst->pres);
	    (void) printf("Phi: %f\n",newst->phi);
	    (void) printf("q: %f\n",newst->q);
	}
}	/* end iF_flowThroughBoundaryState3d */

/*
//OLD VERSION
static void iF_flowThroughBoundaryState2d(
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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
	double dir[MAXD];
	double u[3];		// velocity in the sweeping direction
	double v[3][MAXD];	// velocity in the orthogonal direction
	double vort[3];		// vorticity stencil
	double pres[3];		// pressure stencil
	double f_u;		    // u flux in the sweeping direction
	double f_v[MAXD];	// v flux in the orthogonal direction
	double f_vort;		// vort flux
	double f_pres;		// pressure flux
	double dn,dt = front->dt;
	STATE *oldst,*newst = (STATE*)state;
	STATE  **sts;
	POINTER sl,sr;
	int i,j,dim = front->rect_grid->dim;
	int nrad = 2;
	
	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState2d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	
    for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(nsten->nor,front);
	
    if (debugging("flow_through"))
	{
	    (void) printf("Normal grid size = %f\n",dn);
	    (void) print_Nor_stencil(front,nsten);
	}

	if (comp == negative_component(hs))  
	    oldst = (STATE*)sl;
	else 
	    oldst = (STATE*)sr;

	u[1] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    vort[j] = oldst->vort;
	    pres[j] = oldst->pres;
	}

	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
        FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
                field->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
        
        u[1] += vtmp*dir[i];
	    newst->vel[i] = vtmp;
        newst->vel_old[i] = oldst->vel[i];
	}

	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->vort,
                            getStateVort,&vort[2],&oldst->vort);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2],&oldst->pres);

	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	newst->vort = oldst->vort - dt/dn*f_vort;
	newst->pres = oldst->pres - dt/dn*f_pres;

	tsten = FrontGetTanStencils(front,oldp,nrad);

	if (debugging("flow_through"))
	{
	    (void) printf("Ambient component: %d\n",comp);
	    (void) printf("Tangential grid size = %f\n",dn);
	    (void) print_Tan_stencil(front,tsten[0]);
	}

	if (comp == negative_component(hs))  
	    sts = (STATE**)tsten[0]->leftst;
	else 
	    sts = (STATE**)tsten[0]->rightst;

	for (i = 0; i < dim; ++i)
	    dir[i] = tsten[0]->dir[i];
	dn = FT_GridSizeInDir(dir,front);

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 3; ++j)
	{
        //TODO: Why doesn't this seg fault on j = 0 --> sts[-1]?
	    vort[j] = sts[j-1]->vort;
	    pres[j] = sts[j-1]->pres;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += sts[j-1]->vel[i]*dir[i];
	    }
	    for (i = 0; i < dim; ++i)
	    {
		v[j][i] = sts[j-1]->vel[i] - dir[i]*u[j];
	    }
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]);
	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;

    
    //TODO: Is this reasonably correct? -- see rationale below
    newst->phi = newst->pres;
    newst->q = oldst->pres;
        //newst->q = 0.0;
        //state->q = getQFromPres(front,oldst->pres);
        //newst->phi -= newst->q;
    
    //TODO: add a conditional for incorporating q,
    //      when lagged pressure scheme used.
    
    //Since pressure is usually updated as
    //      
    //      p^{n+1/2} = q + phi^{n+1} - 0.5*mu*(Div_U);
    //
    //      set 
    //
    //      phi^{n+1} = p^{n+1/2} - q    (Div_U = 0 at the boundary) 

    //TODO: check div(u*) = 0 valid for outflow??
    

    if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) printf("Vorticity: %f\n",newst->vort);
	    (void) printf("Pressure: %f\n",newst->pres);
	    (void) printf("Phi: %f\n",newst->phi);
	    (void) printf("q: %f\n",newst->q);
	}
}*/       /* end iF_flowThroughBoundaryState2d */

/*
//OLD VERSION
static void iF_flowThroughBoundaryState3d(
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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
	double dir[MAXD];
	double u[3];		// velocity in the sweeping direction
	double v[3][MAXD];	// velocity in the orthogonal direction
	double pres[3];		// pressure stencil
	double f_u;		    // u flux in the sweeping direction
	double f_v[MAXD];	// v flux in the orthogonal direction
	double f_pres;		// pressure flux
	double dn;

    double dt = front->dt;
	STATE *oldst, *newst = (STATE*)state;
	STATE  **sts;
	POINTER sl,sr;
	int i,j,k,dim = front->rect_grid->dim;
	int nrad = 2;


	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState3d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
	oldst = nullptr;
	if (comp == negative_component(hs))  
	    oldst = (STATE*)sl;
	else 
	    oldst = (STATE*)sr;

    //Normal
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

	u[1] = 0.0;
	
    for (j = 0; j < 2; ++j)
	{
	    pres[j] = oldst->pres;
	}

	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			field->vel[i],getStateVel[i],&vtmp,&oldst->vel[i]);
	    u[1] += vtmp*dir[i];
	    newst->vel[i] = vtmp;
        newst->vel_old[i] = oldst->vel[i];
	}
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2],&oldst->pres);

	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	newst->pres = oldst->pres - dt/dn*f_pres;

	if (debugging("flow_through"))
	{
	    (void) print_Nor_stencil(front,nsten);
	    (void) printf("new velocity after normal prop: %f %f %f\n",
			newst->vel[0],newst->vel[1],newst->vel[2]);
	}
	
    //Tangential
	tsten = FrontGetTanStencils(front,oldp,nrad);
	if (tsten == NULL) return; //TODO: Does this always exit here (see next TODO)?

	for (k = 0; k < dim-1; ++k)
	{
	    for (i = 0; i < dim; ++i)
	    	dir[i] = tsten[k]->dir[i];
	    dn = FT_GridSizeInDir(dir,front);

	    if (comp == negative_component(hs))  
	    	sts = (STATE**)tsten[k]->leftst;
	    else 
	    	sts = (STATE**)tsten[k]->rightst;

	    if (debugging("flow_through"))
	    {
	    	printf("Ambient component: %d\n",comp);
	    	printf("Tangential grid size = %f\n",dn);
		    printf("For direction %d\n",k);
	    	print_Tan_stencil(front,tsten[k]);
	    }

	    for (j = 0; j < 3; ++j) u[j] = 0.0;

	    for (j = 0; j < 3; ++j)
	    {
	    	pres[j] = sts[j-1]->pres;

            for (i = 0; i < dim; ++i)
            {
                u[j] += sts[j-1]->vel[i]*dir[i];
            }

            for (i = 0; i < dim; ++i)
	    	{
                v[j][i] = sts[j-1]->vel[i] - u[j]*dir[i];
	    	}
	    }

	    f_u = burger_flux(u[0],u[1],u[2]);
	    for (i = 0; i < dim; ++i)
	    {
	    	f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	    }
	    f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	    for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]);
	    }
	    newst->pres += - dt/dn*f_pres;
	}

    //TODO: Is this reasonably correct? -- see rationale below
    newst->phi = newst->pres;
    newst->q = oldst->pres;
        //newst->q = 0.0;
        //state->q = getQFromPres(front,oldst->pres);
        //newst->phi -= newst->q;

    
    //TODO: add a conditional for incorporating q,
    //      when lagged pressure scheme used.

    
    //Since pressure is usually updated as
    //      
    //      p^{n+1/2} = q + phi^{n+1} - 0.5*mu*(Div_U);
    //
    //      set 
    //
    //      phi^{n+1} = p^{n+1/2} - q    (Div_U = 0 at the boundary) 

    //TODO: check div(u*) = 0 valid for outflow??
    

    if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) printf("Pressure: %f\n",newst->pres);
	    (void) printf("Phi: %f\n",newst->phi);
	    (void) printf("q: %f\n",newst->q);
	}
}*/	/* end iF_flowThroughBoundaryState3d */

extern void ifluid_point_propagate(
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
        case SUBDOMAIN_BOUNDARY:
            return;
	case MOVABLE_BODY_BOUNDARY:
	case ICE_PARTICLE_BOUNDARY:
	    return rgbody_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case NEUMANN_BOUNDARY:
	case GROWING_BODY_BOUNDARY:
	    return neumann_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case DIRICHLET_BOUNDARY:
	    return dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	default:
	    return contact_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	}
}       /* ifluid_point_propagate */

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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        int i, dim = front->rect_grid->dim;
	double *m_pre = field->pres;
	double *m_phi = field->phi;
	double *m_vor = field->vort;
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	double dn,*h = front->rect_grid->h;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (ifluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}

    //TODO: Interpolate viscosity from nearby like it is
    //      done for the pressure (below).
	setStateViscosity(iFparams,newst,comp);
	
    FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;

	for (i = 0; i < dim; ++i)
	{
        Coords(newp)[i] = Coords(oldp)[i];
	    newst->vel[i] = 0.0;
        FT_RecordMaxFrontSpeed(i,0.0,NULL,Coords(newp),front);
	}
	
    FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres,&oldst->pres);
    
    FT_IntrpStateVarAtCoords(front,comp,p1,m_phi,
            getStatePhi,&newst->phi,&oldst->phi);

    //newst->phi = newst->pres;
    newst->q = oldst->pres;
    //newst->q = 0.0;
	
    if (dim == 2)
    {
	    FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort,&oldst->vort);
	}
	FT_RecordMaxFrontSpeed(dim,0.0,NULL,Coords(newp),front);
        return;
}	/* end neumann_point_propagate */

static  void dirichlet_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double speed;
        int i, dim = front->rect_grid->dim;
	STATE *newst = NULL;
	STATE *bstate;
	COMPONENT comp;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;

	if (debugging("dirichlet_bdry"))
	{
	    (void) printf("Entering dirichlet_point_propagate()\n");
	    (void) print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
	}

	if (ifluid_comp(negative_component(oldhs)))
	{
	    newst = (STATE*)left_state(newp);
	    comp = negative_component(oldhs);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    newst = (STATE*)right_state(newp);
	    comp = positive_component(oldhs);
	}
    
    if (newst == NULL) return;	// node point
    
    //TODO: Interpolate viscosity from nearby like it is
    //      done for the pressure in neumann_point_propagate()???
    setStateViscosity(iFparams,newst,comp);

    //Constant State
	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
        for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] = bstate->vel[i];
		    FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),
                    NULL,Coords(newp),front);
	    }
	    speed = mag_vector(newst->vel,dim);
	    FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
        newst->vort = 0.0;

        newst->pres = 0.0;
        newst->phi = 0.0;
	    newst->q = 0.0;
            //newst->pres = bstate->pres;
            //newst->phi = bstate->pres;
	        //newst->q = bstate->pres;

	    if (debugging("dirichlet_bdry"))
	    {
            printf("Preset boundary state:\n");
            print_general_vector("Velocity: ",newst->vel,dim,"\n");
            //printf("Vorticity: %f\n",newst->vort);
            //printf("Pressure: %f\n",newst->pres);
	    }
	}
	else if (boundary_state_function(oldhs))
	{
	    if (strcmp(boundary_state_function_name(oldhs),
		       "flowThroughBoundaryState") == 0)
	    {
            FLOW_THROUGH_PARAMS ft_params;
            oldp->hse = oldhse;
            oldp->hs = oldhs;
	    	ft_params.oldp = oldp;
            ft_params.comp = comp;
	    	(*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
                    (POINTER)&ft_params,(POINTER)newst);	
	    }
	    else if (strcmp(boundary_state_function_name(oldhs),
		       "iF_timeDependBoundaryState") == 0)
	    {
		TIME_DEPENDENT_PARAMS *td_params = (TIME_DEPENDENT_PARAMS*)
				boundary_state_function_params(oldhs);
	    	(*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
                    (POINTER)td_params,(POINTER)newst);	
	    }
	    else if (strcmp(boundary_state_function_name(oldhs),
		       "iF_splitBoundaryState") == 0)
	    {
		SPLIT_STATE_PARAMS *sparams = (SPLIT_STATE_PARAMS*)
				boundary_state_function_params(oldhs);
	    	(*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
                    (POINTER)sparams,(POINTER)newst);	
	    }
            for (i = 0; i < dim; ++i)
		FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
					front);
	    speed = mag_vector(newst->vel,dim);
	    FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
	}
	if (debugging("dirichlet_bdry"))
	    (void) printf("Leaving dirichlet_point_propagate()\n");
        return;
}	/* end dirichlet_point_propagate */

//TODO: Need to set phi here????
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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double *m_pre = field->pres;
	double *m_vor = field->vort;
	double *p0;
	STATE *oldst,*newst;
	POINTER sl,sr;
	double pres,vort;

        (*front->vfunc)(front->vparams,front,oldp,oldhse,oldhs,vel);
        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	oldst = (STATE*)sl;
	p0 = Coords(newp);
	FT_IntrpStateVarAtCoords(front,-1,p0,m_pre,getStatePres,&pres,
				&oldst->pres);
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,-1,p0,m_vor,getStateVort,&vort,
				&oldst->vort);
	}

	newst = (STATE*)left_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	setStateViscosity(iFparams,newst,negative_component(oldhs));
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	}
	newst = (STATE*)right_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	setStateViscosity(iFparams,newst,positive_component(oldhs));
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	}

	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
    double vel[MAXD];
    int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pre = field->pres;
	double *m_phi = field->phi;
	double *m_vor = field->vort;
	double *m_temp = field->temperature;
	double *m_mu = field->mu;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (ifluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}

    //TODO: Interpolate viscosity from nearby like it is
    //      done for the pressure (below)?
	setStateViscosity(iFparams,newst,comp);
	
    FT_NormalAtPoint(oldp,front,nor,comp);
	dn = grid_size_in_direction(nor,h,dim);

	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;

    if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
    {
        if(!debugging("collision_off"))
        {
            for (i = 0; i < dim; ++i)
            {
                newst->x_old[i] = Coords(oldp)[i];
            }
        }
        
        double omega_dt,crds_com[MAXD];
        omega_dt = angular_velo(oldhs)*dt;

        for (i = 0; i < dim; ++i)
        {
            vel[i] = center_of_mass_velo(oldhs)[i];
            crds_com[i] = Coords(oldp)[i]
                + 0.5*(vel[i] + oldst->vel[i])*dt - rotation_center(oldhs)[i];
        }

        if (dim == 2)
        {
            vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
                angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
            vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
                angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
            
            for (i = 0; i < dim; ++i)
            {
                Coords(newp)[i] =
                    Coords(oldp)[i] + dt*(vel[i] + oldst->vel[i])*0.5;
                newst->vel[i] = vel[i];
                FT_RecordMaxFrontSpeed(i,fabs(vel[i]),
                        NULL,Coords(newp),front);
            }
        }
	    else if (dim == 3)
	    {
            vel[0] += -p_angular_velo(oldhs)[2] * crds_com[1]
                        + p_angular_velo(oldhs)[1] * crds_com[2];
            vel[1] += p_angular_velo(oldhs)[2] * crds_com[0]
                        - p_angular_velo(oldhs)[0] * crds_com[2];
            vel[2] += -p_angular_velo(oldhs)[1] * crds_com[0]
                        + p_angular_velo(oldhs)[0] * crds_com[1];

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
            {
                for (i = 0; i < dim; ++i)
                {
                    Coords(newp)[i] = Coords(oldp)[i]
                        + 0.5*(vel[i] + oldst->vel[i])*dt;
                }
            }
    
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
        fourth_order_point_propagate(front,NULL,oldp,newp,
                oldhse,oldhs,dt,vel);
    }
    
    for (i = 0; i < dim; ++i)
    {
        newst->vel[i] = vel[i];
    }

	FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres,&oldst->pres);
	
	FT_IntrpStateVarAtCoords(front,comp,p1,m_phi,
			getStatePhi,&newst->phi,&oldst->phi);
	
    if (m_temp != NULL)
    {
        FT_IntrpStateVarAtCoords(front,comp,p1,m_temp,getStateTemp,
                &newst->temperature,&oldst->temperature);
    }

    if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort,&oldst->vort);
	}
    
    if(!debugging("collision_off"))
    {
        /* copy newst to the other STATE; used in collision solver */
        if (ifluid_comp(negative_component(oldhs)))
            std::copy(newst, newst+1, (STATE*)right_state(newp));
        else if (ifluid_comp(positive_component(oldhs)))
            std::copy(newst, newst+1, (STATE*)left_state(newp));
    }
}	/* end rgbody_point_propagate */

extern void fluid_print_front_states(
	FILE *outfile,
	Front *front)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int dim = intfc->dim;

	fprintf(outfile,"Interface ifluid states:\n");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",getStatePres(sl),
                                getStatePres(sr));
            fprintf(outfile,"%24.18g %24.18g\n",getStatePhi(sl),
                                getStatePhi(sr));
            if (dim == 2)
            {
                fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                                getStateXvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                                getStateYvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateVort(sl),
                                getStateVort(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateXimp(sl),
                                getStateXimp(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYimp(sl),
                                getStateYimp(sr));
            }
            if (dim == 3)
            {
                fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                                getStateXvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                                getStateYvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateZvel(sl),
                                getStateZvel(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateXimp(sl),
                                getStateXimp(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateYimp(sl),
                                getStateYimp(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateZimp(sl),
                                getStateZimp(sr));
            }
        }
}	/* end fluid_print_front_states */

extern void fluid_read_front_states(
	FILE *infile,
	Front *front)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *lstate,*rstate;
	int dim = intfc->dim;

	next_output_line_containing_string(infile,"Interface ifluid states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            lstate = (STATE*)sl;        rstate = (STATE*)sr;
            fscanf(infile,"%lf %lf",&lstate->pres,&rstate->pres);
            fscanf(infile,"%lf %lf",&lstate->phi,&rstate->phi);
            fscanf(infile,"%lf %lf",&lstate->vel[0],&rstate->vel[0]);
            fscanf(infile,"%lf %lf",&lstate->vel[1],&rstate->vel[1]);
            if (dim == 2)
                fscanf(infile,"%lf %lf",&lstate->vort,&rstate->vort);
            if (dim == 3)
                fscanf(infile,"%lf %lf",&lstate->vel[2],&rstate->vel[2]);
            fscanf(infile,"%lf %lf",&lstate->impulse[0],&rstate->impulse[0]);
            fscanf(infile,"%lf %lf",&lstate->impulse[1],&rstate->impulse[1]);
            if (dim == 3)
            	fscanf(infile,"%lf %lf",&lstate->impulse[2],&rstate->impulse[2]);
        }
}	/* end fluid_read_front_states */

extern void read_iFparams(
	char *inname,
	IF_PARAMS *iFparams)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	int i,dim = iFparams->dim;

	/* defaults numerical schemes */
	iFparams->num_scheme.projc_method = SIMPLE;
	iFparams->num_scheme.advec_method = WENO;
	iFparams->num_scheme.ellip_method = SIMPLE_ELLIP;
    
    if (CursorAfterStringOpt(infile,
        "Entering yes to turn off fluid solver: "))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
        {
            fclose(infile);
            return;
        }
    }

    (void) printf("The default advection order is WENO-Runge-Kutta 4\n");
	iFparams->adv_order = 4;
	if (CursorAfterStringOpt(infile,"Enter advection order:"))
	{
	    fscanf(infile,"%d",&iFparams->adv_order);
	    (void) printf("%d\n",iFparams->adv_order);
	}
	iFparams->extrapolate_advection = false;
	if (CursorAfterStringOpt(infile,"Enter yes for advection term extrapolation:"))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
            iFparams->extrapolate_advection = true;
    }

    CursorAfterString(infile,"Enter projection type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'S':
	case 's':
	    iFparams->num_scheme.projc_method = SIMPLE;
	    break;
	case 'B':
	case 'b':
	case '1':
	    iFparams->num_scheme.projc_method = PMI;
	    //iFparams->num_scheme.projc_method = BELL_COLELLA;
	    break;
	case 'K':
	case 'k':
	case '2':
	    iFparams->num_scheme.projc_method = PMII;
	    //iFparams->num_scheme.projc_method = KIM_MOIN;
	    break;
	case 'P':
	case 'p':
	case '3':
	    iFparams->num_scheme.projc_method = PMIII;
	    //iFparams->num_scheme.projc_method = PEROT_BOTELLA;
	}
	assert(iFparams->num_scheme.projc_method != ERROR_PROJC_SCHEME);
	
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
	    	iFparams->num_scheme.ellip_method = SIMPLE_ELLIP;
	    	break;
	    case 'd':
	    case 'D':
            iFparams->num_scheme.ellip_method = DOUBLE_ELLIP;
	    	break;
        default:
            printf("Elliptic Method Not Implemented\n");
            clean_up(1);
	    }
	}

    //TODO: Should move this into the ambient_state() function
	for (i = 0; i < dim; ++i) iFparams->U_ambient[i] = 0.0;
    if (CursorAfterStringOpt(infile,"Enter fluid ambient velocity:"))
    {
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf ",&iFparams->U_ambient[i]);
            (void) printf("%f ",iFparams->U_ambient[i]);
        }
        (void) printf("\n");
    }
	
    iFparams->ub_speed = HUGE;
    if (CursorAfterStringOpt(infile,"Enter upper bound for speed:"))
	{
            fscanf(infile,"%lf ",&iFparams->ub_speed);
            (void) printf("%f\n",iFparams->ub_speed);
	}
	iFparams->total_div_cancellation = NO;
        if (CursorAfterStringOpt(infile,	
		"Enter yes to use total divergence cancellation:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    	iFparams->total_div_cancellation = YES;
	}
        if (CursorAfterStringOpt(infile,
		"Enter density and viscosity of the fluid:"))
        {
            fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
            (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	}

	iFparams->use_eddy_visc = NO;
    if (CursorAfterStringOpt(infile,
                "Enter yes to use eddy viscosity:"))
    {
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	iFparams->use_eddy_visc = YES;

            printf("Available turbulence models are:\n");
            printf("\tBaldwin-Lomax (B)\n");
            printf("\tVreman (V)\n");
            printf("\tKEPSILON (K)\n");
        	
            CursorAfterString(infile,"Enter turbulence model:");
	    	fscanf(infile,"%s",string);
	    	printf("%s\n",string);
		
            switch (string[0])
            {
                case 'b':
                case 'B':
                    iFparams->eddy_visc_model = BALDWIN_LOMAX;
                    CursorAfterString(infile,"Enter maximum distance for eddy viscosity:");
                    fscanf(infile,"%lf",&iFparams->ymax);
                    printf("%f\n",iFparams->ymax);
                    break;
                case 'm'://keeping for backwards compatibility with old input files
                case 'M'://keeping for backwards compatibility with old input files
                case 'v':
                case 'V':
                    iFparams->eddy_visc_model = VREMAN;
                    CursorAfterString(infile,"Enter model constant:");
                    fscanf(infile,"%lf",&iFparams->C_v);
                    printf("%f\n",iFparams->C_v);
                    break;
                case 'S':
                case 's':
                    iFparams->eddy_visc_model = SMAGORINSKY;
                    CursorAfterString(infile,"Enter model constant:");
                    fscanf(infile,"%lf",&iFparams->C_s);
                    printf("%f\n",iFparams->C_s);
                    break;
                case 'K':
                case 'k':
                    iFparams->eddy_visc_model = KEPSILON;
                    break;
                default:
                    (void) printf("Unknown eddy viscosity model!\n");
                    clean_up(ERROR);
            }

            //      "Enter yes to use slip wall boundary condition:"
            iFparams->use_no_slip = YES;
            if (CursorAfterStringOpt(infile,"Enter yes to use no-slip boundary condition:"))
            {
                fscanf(infile,"%s",string);
                printf("%s\n",string);
                //if (string[0] == 'y' || string[0] == 'Y')
                if (string[0] == 'n' || string[0] == 'N')
                {
                    iFparams->use_no_slip = NO;
                    //iFparams->use_no_slip = YES;
                }
            }
	    }
	}
   
    
    //TODO: Need these here? or gets handled in iFinit.cpp for 2 phase flow problems?
    CursorAfterStringOpt(infile,"Enter surface tension:");
    fscanf(infile,"%lf",&iFparams->surf_tension);
    printf("%f\n",iFparams->surf_tension);
    CursorAfterStringOpt(infile,"Enter factor of smoothing radius:");
    fscanf(infile,"%lf",&iFparams->smoothing_radius);
    printf("%f\n",iFparams->smoothing_radius);

    for (i = 0; i < dim; ++i) iFparams->gravity[i] = 0.0;
    if (CursorAfterStringOpt(infile,"Enter gravity:"))
    {
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf ",&iFparams->gravity[i]);
            (void) printf("%f ",iFparams->gravity[i]);
        }
        (void) printf("\n");
    }

    iFparams->scalar_field = NO;
    if (CursorAfterStringOpt(infile,"Enter yes to consider scalar field:"))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
            iFparams->scalar_field = YES;
    }

    iFparams->min_speed = 0.0;
    if (CursorAfterStringOpt(infile,
        "Enter minimum speed to limit time step:"))
    {
        fscanf(infile,"%lf ",&iFparams->min_speed);
        (void) printf("%f ",iFparams->min_speed);
        (void) printf("\n");
    }
	fclose(infile);
}	/* end read_iFparams */

extern boolean isDirichletPresetBdry(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp)
{
	HYPER_SURF *hs;
	POINTER intfc_state;
	double crx_coords[MAXD];
	INTERFACE *grid_intfc = front->grid_intfc;

	if (!FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
                                comp,&intfc_state,&hs,crx_coords))
	    return NO;
	if (wave_type(hs) != DIRICHLET_BOUNDARY)
	    return NO;
	if (boundary_state(hs) == nullptr)
	    return NO;
	return YES;
}	/* end isDirichletPresetBdry */

static void get_split_state_params(
	Front *front,
	FILE *infile,
	POINTER *params)
{
	static SPLIT_STATE_PARAMS *split_st_params;
	char string[100];
	int k;
	int dim = FT_Dimension();

	FT_ScalarMemoryAlloc((POINTER*)&split_st_params,
                        sizeof(SPLIT_STATE_PARAMS));
	*params = (POINTER)split_st_params;

	CursorAfterString(infile,"Enter direction of split:");
	fscanf(infile,"%d",&split_st_params->dir);
	(void) printf(" %d\n",split_st_params->dir);
	CursorAfterString(infile,"Enter coordinate of split:");
	fscanf(infile,"%lf",&split_st_params->split_coord);
	(void) printf(" %f\n",split_st_params->split_coord);

	CursorAfterString(infile,"For the left state");
	CursorAfterString(infile,"Enter velocity: ");
	for (k = 0; k < dim; ++k)
	{
	    fscanf(infile,"%lf",&split_st_params->left_state.vel[k]);
	    (void) printf("%f ",split_st_params->left_state.vel[k]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Enter pressure:");
	fscanf(infile,"%lf",&split_st_params->left_state.pres);
	(void) printf("%f\n",split_st_params->left_state.pres);
	if (CursorAfterStringOpt(infile,"Enter temperature:"))
	{
	    fscanf(infile,"%lf",&split_st_params->left_state.temperature);
	    (void) printf("%f\n",split_st_params->left_state.temperature);
	}
	split_st_params->left_state.phi = getPhiFromPres(front,
			split_st_params->left_state.pres);

	CursorAfterString(infile,"For the right state");
	CursorAfterString(infile,"Enter velocity: ");
	for (k = 0; k < dim; ++k)
	{
	    fscanf(infile,"%lf",&split_st_params->right_state.vel[k]);
	    (void) printf("%f ",split_st_params->right_state.vel[k]);
	}
	(void) printf("\n");
	CursorAfterString(infile,"Enter pressure:");
	fscanf(infile,"%lf",&split_st_params->right_state.pres);
	(void) printf("%f\n",split_st_params->right_state.pres);
	if (CursorAfterStringOpt(infile,"Enter temperature:"))
	{
	    fscanf(infile,"%lf",&split_st_params->right_state.temperature);
	    (void) printf("%f\n",split_st_params->right_state.temperature);
	}
	split_st_params->right_state.phi = getPhiFromPres(front,
			split_st_params->right_state.pres);
}	/* end get_split_state_params */

static void get_parabolic_state_params(
	Front *front,
	FILE *infile,
	POINTER *params)
{
	static PARABOLIC_STATE_PARAMS *parab_st_params;
        char string[100];
        int k;
        int dim = FT_Dimension();

        FT_ScalarMemoryAlloc((POINTER*)&parab_st_params,
                        sizeof(PARABOLIC_STATE_PARAMS));
        *params = (POINTER)parab_st_params;

        CursorAfterString(infile,"Enter direction:");
        fscanf(infile,"%d",&parab_st_params->dir);
        (void) printf("%d\n",parab_st_params->dir);

        CursorAfterString(infile,"Enter peak velocity: ");
        for (k = 0; k < dim; ++k)
        {
            fscanf(infile,"%lf",&parab_st_params->v_peak[k]);
            (void) printf("%f ",parab_st_params->v_peak[k]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter pressure:");
        fscanf(infile,"%lf",&parab_st_params->state.pres);
        (void) printf("%f\n",parab_st_params->state.pres);
        if (CursorAfterStringOpt(infile,"Enter temperature:"))
        {
            fscanf(infile,"%lf",&parab_st_params->state.temperature);
            (void) printf("%f\n",parab_st_params->state.temperature);
        }

        parab_st_params->state.phi = getPhiFromPres(front,
                        parab_st_params->state.pres);
}	/* end get_parabolic_state_params */

static void get_time_dependent_params(
	int dim,
	FILE *infile,
	POINTER *params)
{
	static TIME_DEPENDENT_PARAMS *td_params;
	char string[100];
	int i;

	FT_ScalarMemoryAlloc((POINTER*)&td_params,
			sizeof(TIME_DEPENDENT_PARAMS));
	CursorAfterString(infile,"Enter type of time-dependent function:");
	fscanf(infile,"%s",string);
	(void) printf(" %s\n",string);
	switch (string[0])
	{
	case 'C':
	case 'c':
	    td_params->td_type = CONSTANT;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    break;
	case 'P':
	case 'p':
	    td_params->td_type = PULSE_FUNC;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    CursorAfterString(infile,"Enter peak velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_peak[i]);
		(void) printf("%f ",td_params->v_peak[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter peak pressure:");
	    fscanf(infile,"%lf ",&td_params->p_peak);
	    (void) printf("%f\n",td_params->p_peak);
	    CursorAfterString(infile,"Enter tail velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_tail[i]);
		(void) printf("%f ",td_params->v_tail[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter tail pressure:");
	    fscanf(infile,"%lf ",&td_params->p_tail);
	    (void) printf("%f\n",td_params->p_tail);
	    CursorAfterString(infile,"Enter time to rise:");
	    fscanf(infile,"%lf ",&td_params->T[0]);
	    (void) printf("%f\n",td_params->T[0]);
	    CursorAfterString(infile,"Enter time to reach peak:");
	    fscanf(infile,"%lf ",&td_params->T[1]);
	    (void) printf("%f\n",td_params->T[1]);
	    CursorAfterString(infile,"Enter time to fall:");
	    fscanf(infile,"%lf ",&td_params->T[2]);
	    (void) printf("%f\n",td_params->T[2]);
	    CursorAfterString(infile,"Enter time to reach tail:");
	    fscanf(infile,"%lf ",&td_params->T[3]);
	    (void) printf("%f\n",td_params->T[3]);
	    break;
	case 'S':
	case 's':
	    td_params->td_type = SINE_FUNC;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    CursorAfterString(infile,"Enter velocity amplitude:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_amp[i]);
		(void) printf("%f ",td_params->v_amp[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter pressure amplitude:");
	    fscanf(infile,"%lf ",&td_params->p_amp);
	    (void) printf("%f\n",td_params->p_amp);
	    CursorAfterString(infile,"Enter oscillation period:");
	    fscanf(infile,"%lf ",&td_params->omega);
	    (void) printf("%f\n",td_params->omega);
	    td_params->omega = 2.0*PI/td_params->omega;
	    CursorAfterString(infile,"Enter initial phase:");
	    fscanf(infile,"%lf ",&td_params->phase);
	    (void) printf("%f\n",td_params->phase);
	    td_params->phase *= PI/180.0;
	    break;
	default:
	    (void) printf("Unknown type of time-dependent function!\n");
	    clean_up(ERROR);
	}

	*params = (POINTER)td_params;	
}	/* end get_time_dependent_params */

//TODO: see addToEnergyFlux() below
extern void recordBdryEnergyFlux(
	Front *front,
	char *out_name)
{
	RECT_GRID *rgr = front->rect_grid;
	INTERFACE *intfc = front->interf;
	HYPER_SURF *hs;
	double Energy_in,Energy_out;
	int dir,side;
	int dim = rgr->dim;
	static FILE *ein_file,*eout_file;
	char file_name[100];

	if (ein_file == NULL && pp_mynode() == 0)
	{
	    sprintf(file_name,"%s/in_energy.xg",out_name);
	    ein_file = fopen(file_name,"w");
	    fprintf(ein_file,"\"Energy influx\" vs. time\n");
	    sprintf(file_name,"%s/out_energy.xg",out_name);
	    eout_file = fopen(file_name,"w");
	    fprintf(eout_file,"\"Energy outflux\" vs. time\n");
	}
	Energy_in = Energy_out = 0.0;
	for (dir = 0; dir < dim; ++dir)
	for (side = 0; side < 2; ++side)
	{
	    hs = FT_RectBoundaryHypSurf(intfc,DIRICHLET_BOUNDARY,dir,side);
	    if (hs == NULL) continue;
	    if (boundary_state(Hyper_surf(hs)))
	    {
		addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,YES);
	    }
	    else if (boundary_state_function(Hyper_surf(hs)))
	    {
		if (strcmp(boundary_state_function_name(hs),
                       "iF_timeDependBoundaryState") == 0)
		{
		    addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,YES);
		}
		else if (strcmp(boundary_state_function_name(hs),
                       "flowThroughBoundaryState") == 0)
		{
		    addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,NO);
		}
	    }
	}
	pp_global_sum(&Energy_in,1);
	pp_global_sum(&Energy_out,1);

	if (pp_mynode() == 0)
	{
	    fprintf(ein_file,"%f %f\n",front->time,Energy_in);
	    fprintf(eout_file,"%f %f\n",front->time,Energy_out);
	}
}	/* end recordBdryEnergyFlux */


static void addToEnergyFlux(
	RECT_GRID *rgr,
	HYPER_SURF *hs,
	double *Energy_in,
	double *Energy_out,
	int dir,
	int side,
	boolean is_influx)
{
	int i,dim = rgr->dim;
	double *L = rgr->L;	
	double *U = rgr->U;	
	CURVE *c;
	SURFACE *s;
	BOND *b = NULL;
	TRI *t;
	double ave_coord,engy_flux,vel;
	boolean is_outside_hse,is_left_side;
	STATE *sl,*sr,*state;
	POINT *pts[MAXD];

	is_left_side = (ifluid_comp(negative_component(hs))) ? YES : NO;
	switch (dim)
	{
	case 2:
	    c = Curve_of_hs(hs);
	    for (b = c->first; b != NULL; b = b->next)
	    {
		is_outside_hse = NO;
		pts[0] = b->start;
		pts[1] = b->end;
		for (i = 1; i < dim; ++i)
		{
		    ave_coord = (Coords(pts[0])[(dir+i)%dim] +
				 Coords(pts[1])[(dir+i)%dim])/2.0;
		    if (ave_coord < L[(dir+i)%dim] ||
			ave_coord > U[(dir+i)%dim])
			is_outside_hse = YES;
		}
		if (is_outside_hse) continue;
		engy_flux = 0.0;
		for (i = 0; i < 2; ++i)
		{
		    FT_GetStatesAtPoint(pts[i],Hyper_surf_element(b),hs,
				(POINTER*)&sl,(POINTER*)&sr);
		    state = (is_left_side) ? sl : sr;
		    if (is_influx)
		    	vel = (side == 0) ? state->vel[dir] : -state->vel[dir];
		    else
		    	vel = (side == 1) ? state->vel[dir] : -state->vel[dir];
		    //engy_flux += 0.5*state->dens*(sqr(state->vel[0]) +
		    engy_flux += 0.5*(sqr(state->vel[0]) +
					sqr(state->vel[1]))*vel;
		}
		engy_flux *= bond_length(b)/2.0;
		if (is_influx)
		    *Energy_in += engy_flux;
		else
		    *Energy_out += engy_flux;
	    }
	    break;
	case 3:
	    s = Surface_of_hs(hs);
	    for (t = first_tri(s); !at_end_of_tri_list(t,s); t = t->next)
	    {
		is_outside_hse = NO;
		for (i = 0; i < 3; ++i)
		    pts[i] = Point_of_tri(t)[i];
		for (i = 1; i < dim; ++i)
		{
		    ave_coord = (Coords(pts[0])[(dir+i)%dim] +
				 Coords(pts[1])[(dir+i)%dim] +
				 Coords(pts[2])[(dir+i)%dim])/3.0;
		    if (ave_coord < L[(dir+i)%dim] ||
			ave_coord > U[(dir+i)%dim])
			is_outside_hse = YES;
		}
		if (is_outside_hse) continue;
		engy_flux = 0.0;
		for (i = 0; i < 3; ++i)
		{
		    FT_GetStatesAtPoint(pts[i],Hyper_surf_element(t),hs,
				(POINTER*)&sl,(POINTER*)&sr);
		    state = (is_left_side) ? sl : sr;
		    if (is_influx)
		    	vel = (side == 0) ? state->vel[dir] : -state->vel[dir];
		    else
		    	vel = (side == 1) ? state->vel[dir] : -state->vel[dir];
		    //engy_flux += 0.5*state->dens*(sqr(state->vel[0]) +
		    engy_flux += 0.5*(sqr(state->vel[0]) +
					sqr(state->vel[1]) +
					sqr(state->vel[2]))*vel;
		}
		engy_flux *= bond_length(b)/3.0;
		if (is_influx)
		    *Energy_in += engy_flux;
		else
		    *Energy_out += engy_flux;
	    }
	    break;
	}
}	/* end addToEnergyFlux */

//TODO: To remove if not used anywhere
extern double p_jump(
	POINTER params,
	int D,
	double *coords)
{
	return 0.0;
}	/* end p_jump */

extern double grad_p_jump_n(
	POINTER params,
	int D,
	double *N,
	double *coords)
{
	return 0.0;
}	/* end grad_p_jump_n */

extern double grad_p_jump_t(
	POINTER params,
	int D,
	int i,
	double *N,
	double *coords)
{
	return 0.0;
}	/* end grad_p_jump_t */

//TODO: This should be used for updating boundary states only.
//      Should not be used for initialization.
extern double getPhiFromPres(
        Front *front,
        double pres)
{
    IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
    /*
    IF_FIELD* field = iFparams->field;
    double* div_U = field->div_U;
    double* mu = field->mu;
    double* q = field->q;
    */

    switch (iFparams->num_scheme.projc_method)
    {
        case PMI:
        case PMII:
            /*if (!isbdry)
                return pres - q[index] + 0.5*mu[index]*div_U[index];
            else
                return 0.0;*/
            return pres;
        case PMIII:
        case SIMPLE:
            /*if (!isbdry)
                return pres + 0.5*mu[index]*div_U[index];
            else
                return pres;*/
            return 0.0;
        default:
            (void) printf("Unknown projection type\n");
            clean_up(0);
    }
}       /* end getPhiFromPres */

extern double getQFromPres(
        Front *front,
        double pres)
{
    IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;

    switch (iFparams->num_scheme.projc_method)
    {
        case PMI:
        case PMII:
            return pres;
        case PMIII:
        case SIMPLE:
            return 0.0;
        default:
            (void) printf("Unknown projection type\n");
            clean_up(EXIT_FAILURE);
    }
}       /* end getPhiFromPres */

extern double getPressure(
        Front *front,
        double *coords,
        double *base_coords)
{
        INTERFACE *intfc = front->interf;
        int i,dim = Dimension(intfc);
        POINT *p0;
        double pres,pres0;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *g = iFparams->gravity;
        double rho = iFparams->rho2;
        boolean hyper_surf_found = NO;

        return 0.0;
        
        //TODO: Does below work???
        //
        //      It appears this early return of 0.0 may have been
        //      hardcoded when the initialization in parachute.cpp
        //      was being worked on and never set back to a general
        //      mode of operation -- they set the l_cartesian->getInitialState
        //      function pointer to zero_state() which zeros the
        //      pressure just like it is here.
        //
        //      Should attempt to restore this functionality and
        //      experiment with some other initial condtions, which
        //      would include the initial boundary conditions at
        //      the inlet/outlet at appears.

        pres0 = 0.0;
        //pres0 = 1.0;
        
        if (dim == 2)
        {
            CURVE **c;
            for (c = intfc->curves; c && *c; ++c)
            {
                if (wave_type(*c) == DIRICHLET_BOUNDARY &&
                    boundary_state(*c) != NULL)
                {
                    p0 = (*c)->first->start;
                    pres0 = getStatePres(boundary_state(*c));
                    hyper_surf_found = YES;
                    break;
                }
            }
        }
        else if (dim == 3)
        {
            SURFACE **s;
            for (s = intfc->surfaces; s && *s; ++s)
            {
                if (wave_type(*s) == DIRICHLET_BOUNDARY &&
                    boundary_state(*s) != NULL)
                {
                    p0 = Point_of_tri(first_tri(*s))[0];
                    pres0 = getStatePres(boundary_state(*s));
                    hyper_surf_found = YES;
                    break;
                }
            }
        }
        
        pres = pres0;
        //return pres;
        
        //TODO: 
        if (hyper_surf_found)
        {
            //NOTE: Assume g = {0,0,-9.8} is standard gravity.
            //      Then if the inlet pressure is prescribed at the
            //      lower z boundary, points in the domain above
            //      will have lower pressure.
            for (i = 0; i < dim; ++i)
                pres += rho*(coords[i] - Coords(p0)[i])*g[i];
                //pres -= rho*(coords[i] - Coords(p0)[i])*g[i];
        }
        else if (base_coords != NULL)
        {
            for (i = 0; i < dim; ++i)
                pres += rho*(coords[i] - Coords(p0)[i])*g[i];
                //pres -= rho*(coords[i] - Coords(p0)[i])*g[i];
        }

        return pres;
}       /* end getPressure */

extern int ifluid_find_state_at_dual_crossing(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	int comp,
	POINTER *state,
	HYPER_SURF **hs,
	double *crx_coords)
{
	boolean status;
	INTERFACE *grid_intfc = front->comp_grid_intfc;
	status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
				comp,state,hs,crx_coords);
	if (status == NO) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == NEUMANN_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == MOVABLE_BODY_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == ICE_PARTICLE_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY) 
	{
	    if (boundary_state(*hs))
	    	return CONST_V_PDE_BOUNDARY;
	    else
	    	return CONST_P_PDE_BOUNDARY;
	}
}	/* ifluid_find_state_at_crossing */

extern  void ifluid_compute_force_and_torque(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        switch (fr->rect_grid->dim)
        {
        case 2:
            return ifluid_compute_force_and_torque2d(fr,hs,dt,force,torque);
        case 3:
            return ifluid_compute_force_and_torque3d(fr,hs,dt,force,torque);
        }
}       /* end ifluid_compute_force_and_torque */

static  void ifluid_compute_force_and_torque2d(
        Front *front,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = computational_grid(front->interf);
        double f[MAXD],rr[MAXD];
        double t,pres;
        double posn[MAXD],bnor[MAXD];
        double area;
        BOND *b;
        boolean pos_side;
        int i,dim = gr->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *gravity = iFparams->gravity;
        CURVE *curve = Curve_of_hs(hs);

        if (debugging("rigid_body"))
            (void) printf("Entering ifluid_compute_force_and_torque2d()\n");

        if (ifluid_comp(negative_component(curve)))
            pos_side = NO;
        else
            pos_side = YES;

        for (i = 0; i < dim; ++i)
        {
            force[i] = 0.0;
        }
        *torque = 0.0;

	//if (front->step > 5)
	if (front->step > iFparams->fsi_startstep)
	{
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
                        //NOTE: posn = 0.5*(Coords(b->start)[i] + Coords(b->end)[i])
                        //      in most cases.
                        force[i] += f[i];
                    }
                    Cross2d(rr,f,t);
                    *torque += t;
                }
            }
	}

    //TODO: Need to add rotational motion
         
        /* Add gravity to the total force */
        if (motion_type(curve) != ROTATION)
        {
            for (i = 0; i < dim; ++i)
                force[i] += gravity[i]*total_mass(curve);
        }

        if (debugging("rigid_body"))
        {
            (void) printf("Leaving ifluid_compute_force_and_torque2d()\n");
            (void) printf("total_force = %f %f\n",force[0],force[1]);
            (void) printf("torque = %f\n",*torque);
        }
}       /* end ifluid_compute_force_and_torque2d */

#define         MAX_TRI_FOR_INTEGRAL            100
static  void ifluid_compute_force_and_torque3d(
        Front *front,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = computational_grid(front->interf);
        double f[MAXD],rr[MAXD];
        double t[MAXD],tdir,pres;
        double posn[MAXD],tnor[MAXD];
        double area;
        TRI *tri;
        boolean pos_side;
        int i,dim = gr->dim;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double *gravity = iFparams->gravity;
        SURFACE *surface = Surface_of_hs(hs);
	CURVE **c;
        NODE *rg_string_nodes[10];
        int j, k, num = 0;
	NODE **n;
	BOND *b;
	double tri_cen[MAXD];

    if (debugging("rigid_body"))
        (void) printf("Entering ifluid_compute_force_and_torque3d()\n"); 
    
    if (ifluid_comp(negative_component(surface)))
        pos_side = NO;
    else
        pos_side = YES;

    for (i = 0; i < dim; ++i)
    {
        force[i] = 0.0;
        torque[i] = 0.0;
    }

	/* count in the force and torque on the RG_STRING_NODE */
	intfc_node_loop(front->interf, n)
	{
	    for (k = 0; k < dim; ++k)
        {
            if (Coords((*n)->posn)[k] <= gr->L[k] ||
                Coords((*n)->posn)[k] > gr->U[k]) break;
        }

        if (k != dim || (*n)->extra == NULL) continue;

        node_out_curve_loop(*n,c)
	    {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) break;
	    }

	    if (c == NULL || (*c) == NULL)
	    {
            node_in_curve_loop(*n,c)
            {
                if (hsbdry_type(*c) == PASSIVE_HSBDRY) break;
            }
	    }

	    if (c == NULL || (*c) == NULL) continue;

	    b = (*c)->first;
	    if (wave_type(b->_btris[0]->surface) == MOVABLE_BODY_BOUNDARY)
        {
            rg_string_nodes[num++] = *n;
        }
	}
        
    for (j = 0; j < num; ++j)
    {
        POINT *p = rg_string_nodes[j]->posn;
        if (!coords_in_subdomain(Coords(p),gr)) continue;

        for (i = 0; i < dim; ++i)
        {
            force[i] += p->force[i];//Computed from setSpecialNodeForce()
            rr[i] = Coords(p)[i] - rotation_center(surface)[i];
        }

        Cross3d(rr, p->force, t);
        for (i = 0; i < dim; ++i)
            torque[i] += t[i];
    
        if (debugging("rigid_body"))
        {
            printf("rg_string_nodes coords = %f %f %f\n", 
                Coords(p)[0], Coords(p)[1], Coords(p)[2]);
            printf("rg_string_nodes force = %f %f %f\n", 
                p->force[0], p->force[1], p->force[2]);
        }
    }
	/* end of counting the force on RG_STRING_NODE */

	//if (front->step > 5)
	if (front->step > iFparams->fsi_startstep)
	{
        for (tri = first_tri(surface); !at_end_of_tri_list(tri,surface); tri = tri->next)
        {
            for (i = 0; i < dim; ++i)
            {
                tri_cen[i] = (Coords(Point_of_tri(tri)[0])[i] +
                      Coords(Point_of_tri(tri)[1])[i] +
                      Coords(Point_of_tri(tri)[2])[i])/3.0;
            }

            if (!coords_in_subdomain(tri_cen,gr)) continue;

            if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
                    &pres,tnor,posn,pos_side))
            {
                area = tri_area(tri);
                    //area = 0.5*Mag3d(tnor);
                double mag_tnor = Mag3d(tnor);
                for (i = 0; i < dim; ++i)
                {
                    f[i] = pres*area*tnor[i]/mag_tnor;
                    force[i] += f[i];
                    rr[i] = posn[i] - rotation_center(surface)[i];
                }
                
                Cross3d(rr,f,t);
                //tdir = Dot3d(t,(rotation_direction(hs)));
                for (i = 0; i < dim; ++i)
                {
                    //t[i] = tdir*rotation_direction(hs)[i];
                    torque[i] += t[i];
                }
            }
        }
	}


        //TODO: force computation should include effects of shear stress from
        //      turbulence model + wall functions (see to computeDiffusionCN() todos).
        
    
        /* Add gravity to the total force */
        if (motion_type(surface) != ROTATION &&
	        motion_type(surface) != PRESET_ROTATION)
        {
            for (i = 0; i < dim; ++i)
                force[i] += gravity[i]*total_mass(surface)/num_clips(surface);
        }

        if (debugging("rigid_body"))
        {
            printf("In ifluid_compute_force_and_torque3d()\n");
            printf("total_force = %f %f %f\n",force[0],force[1],force[2]);
            printf("torque = %f %f %f\n",torque[0],torque[1],torque[2]);
            printf("# of rg_string_node in processor %d = %d\n", pp_mynode(), num);
	        printf("number of clips = %d \n", num_clips(surface));
        }
        
        if (debugging("rigid_body"))
            printf("Leaving ifluid_compute_force_and_torque3d()\n"); 
}       /* end ifluid_compute_force_and_torque3d */

static boolean force_on_hse(
        HYPER_SURF_ELEMENT *hse,        /* Bond (2D) or tri (3D) */
        HYPER_SURF *hs,                 /* Curve (2D) or surface (3D) */
        RECT_GRID *gr,                  /* Rectangular grid */
        double *pres,           /* Average pressure */
        double *nor,           /* normal vector pointing into body */
        double *posn,           /* Position of the pressure */
        boolean pos_side)       /* Is the body on the positive side of hs? */
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

}       /* end force_on_hse */

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

        p1 = getStatePres(s1);  p2 = getStatePres(s2);
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
}       /* end force_on_hse2d */

static boolean force_on_hse3d(
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        RECT_GRID *gr,
        double *pres,
        double *tnor,
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
            //normal points toward the surface
            tnor[i] = pos_side ? -Tri_normal(t)[i] : Tri_normal(t)[i];
            posn[i] /= 3.0;
        }
        /* Need to treat subdomain boundary */
        return YES;
}       /* end force_on_hse3d */

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

static void setStateViscosity(
	IF_PARAMS *iFparams,
	STATE *state,
	int comp)
{
	switch (comp)
	{
	case LIQUID_COMP1:
	    state->mu = iFparams->mu1;
	    break;
	case LIQUID_COMP2:
	    state->mu = iFparams->mu2;
	    break;
	default:
	    state->mu = 0.0;
	}
}

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
	    FT_InsertDirichletBoundary(front,NULL,NULL,
			NULL,(POINTER)state,*hs,i_hs);
	    
        //TODO: Can not prescribe pressure with velocity in current formulation
        state->phi = 0.0;
        state->q = 0.0;
            //state->phi = state->pres;
            //state->q = getQFromPres(front,state->pres);

	    break;
	case 'f':			// Flow through state
	case 'F':
	    FT_InsertDirichletBoundary(front,
			iF_flowThroughBoundaryState,"flowThroughBoundaryState",
			NULL,NULL,*hs,i_hs);
	    break;
	case 't':			// Time dependent state
	case 'T':
	    get_time_dependent_params(dim,infile,&func_params);
	    FT_InsertDirichletBoundary(front,iF_timeDependBoundaryState,
			"iF_timeDependBoundaryState",func_params,NULL,*hs,i_hs);
	    break;
	case 's':			// Split state
	case 'S':
	    get_split_state_params(front,infile,&func_params);
	    FT_InsertDirichletBoundary(front,iF_splitBoundaryState,
			"iF_splitBoundaryState",func_params,NULL,*hs,i_hs);
	    break;
	case 'p':
	case 'P':
	    get_parabolic_state_params(front,infile,&func_params);
	    FT_InsertDirichletBoundary(front,iF_parabolicBoundaryState,
			"iF_splitBoundaryState",func_params,NULL,*hs,i_hs);
	    break;
	default:
	    (void) printf("Unknown Dirichlet boundary!\n");
	    clean_up(ERROR);
	}
}	/* end promptForDirichletBdryState */

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
	IF_PARAMS* ifparams = (IF_PARAMS*)front->extra1;
	IF_FIELD* field = ifparams->field;
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

	if (!ifluid_comp(comp))
	{
	    for (i = 0; i < 3; ++i)
		((STATE*)state)->vel[i] = 0.0;
	    ((STATE*)state)->pres = 0.0;
	    *bdry_type = PASSIVE_BOUNDARY;
	    return;
	}

	for (i = 0; i < 3; ++i)
	{
	    ic_nb[i] = ic[i];
	    coords[i] = top_L[i] + ic[i] * top_h[i];
	    if (i != idir)
	        dist += sqr(coords[i] - c[i]);
	}
	dist = sqrt(dist);
	if (dist < r)
	{
	    *bdry_type = d_params->in_pipe_bdry;
	    switch (*bdry_type)
	    {
	    case DIRICHLET_BOUNDARY:
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

		if (d_params->in_flow_through)
		{
		    ((STATE*)state)->dens = field->rho[index_nb];
		    for (i = 0; i < 3; ++i)
			((STATE*)state)->vel[i] = field->vel[i][index_nb];
		    ((STATE*)state)->pres = 0.0;
		    ((STATE*)state)->phi = 0.0;
		}
		else
		{
		    for (i = 0; i < 3; ++i)
			((STATE*)state)->vel[i] = d_params->state[0].vel[i];
		    ((STATE*)state)->dens = (comp == LIQUID_COMP1) ? 
					    ifparams->rho1 : ifparams->rho2;
		    ((STATE*)state)->pres = field->pres[index_nb]; 
		    ((STATE*)state)->phi = field->phi[index_nb];
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
		((STATE*)state)->dens = field->rho[index_nb];
		for (i = 0; i < 3; ++i)
		{
		    ((STATE*)state)->vel[i] = field->vel[i][index_nb];
		    if (i == idir)
		    {
		        ((STATE*)state)->vel[i] *= -1.0;
		    }
		}
		((STATE*)state)->pres = field->pres[index_nb]; 
		((STATE*)state)->phi = field->phi[index_nb];
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

		if (d_params->out_flow_through)
		{
		    ((STATE*)state)->dens = field->rho[index_nb];
		    for (i = 0; i < 3; ++i)
			((STATE*)state)->vel[i] = field->vel[i][index_nb];
		    ((STATE*)state)->pres = 0.0;
		    ((STATE*)state)->phi = 0.0;
		}
		else
		{
		    for (i = 0; i < 3; ++i)
			((STATE*)state)->vel[i] = d_params->state[1].vel[i];
		    ((STATE*)state)->dens = (comp == LIQUID_COMP1) ? 
					    ifparams->rho1 : ifparams->rho2;
		    ((STATE*)state)->pres = field->pres[index_nb];
		    ((STATE*)state)->phi = field->phi[index_nb];
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
		((STATE*)state)->dens = field->rho[index_nb];
		for (i = 0; i < 3; ++i)
		{
		    ((STATE*)state)->vel[i] = field->vel[i][index_nb];
		    if (i == idir)
		    {
		        ((STATE*)state)->vel[i] *= -1.0;
		    }
		}
		((STATE*)state)->pres = field->pres[index_nb];
		((STATE*)state)->phi = field->phi[index_nb];
		break;
	    default:
		(void) printf("Unknown outer pipe boundary type\n");
		clean_up(ERROR);
	    }
	}

	return;
}

static boolean coords_in_subdomain(
	double *coords,
	RECT_GRID *gr)
{
	int i,dim = gr->dim;
	for (i = 0; i < dim; ++i)
	{
	    if (coords[i] < gr->L[i] || coords[i] >= gr->U[i])
		return NO;
	}
	return YES;
}	/* end coords_in_subdomain */

extern void setContactNodeType(Front *front)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	NODE **n;
	if (iFparams->surf_tension == 0.0) return;
	intfc_node_loop(front->interf,n)
	{
	    if (node_type(*n) == NEUMANN_NODE)
		node_type(*n) = CONTACT_NODE;
	}
}	/* end setContactNodeType */

extern int contact_node_propagate(
	Front           *fr,
        POINTER         wave,
        NODE            *oldn,
        NODE            *newn,
        RPROBLEM        **rp,
        double          dt,
        double          *dt_frac,
        NODE_FLAG       flag)
{
	O_CURVE		Oldcphys, Newcphys;	/* the physical curve */
	O_CURVE		Oldcahead;              /* ahead bdry wrt oldn */
	O_CURVE		Newcahead;              /* crspd of Oldcahead */
	O_CURVE		Oldcbehind;             /* behind bdry wrt oldn */
	O_CURVE		Newcbehind;             /* crspd of Oldcbehind */
	O_CURVE		Newcaprop;              /* ahead bdry wrt newn */
	O_CURVE		Oldcaprop;              /* crspd of Newcaprop */
	O_CURVE		Newcbprop;              /* behind bdry wrt newn */
	O_CURVE		Oldcbprop;              /* crspd of Newcbprop */
	BOND		*crossbphys;		/* intersecting two bonds */
	BOND		*crossbahead;		/* on newcphys, newcahead */
	double		tcr_phys,tcr_ahead;	/* fractional dist to cross */
	ANGLE_DIRECTION	i_to_prop_dir;          /* ang dir - inc to forward
						   wrt dir of motion */
	SIDE		propagation_side;       /* side of cahead defined by
						   newn posn*/
	SIDE		inc_side;               /* side of cahead defined by
						   cphys */
	static POINT	*pc = NULL;		/* crossing point */
	COMPONENT       propagation_comp;	/* comp containing newn */
	COMPONENT	ahead_comp;             /* comp on ahead side of
						   cphys wrt motion */
	static POINT	*newp = NULL;		/* new node propagation point */
	CURVE		*ca;
	int		status;			/* diagnostic return value */
	boolean opposite_dir_tried = NO;
	i_to_prop_dir = ANGLE_DIRECTION_NOT_SET;
	
	if (pc == NULL)
	{
	    pc = Static_point(fr->interf);
	    newp = Static_point(fr->interf);
	}
	zero_scalar(&Oldcphys,sizeof(O_CURVE));
	zero_scalar(&Newcphys,sizeof(O_CURVE));
	zero_scalar(&Oldcahead,sizeof(O_CURVE));
	zero_scalar(&Newcahead,sizeof(O_CURVE));
	zero_scalar(&Oldcbehind,sizeof(O_CURVE));
	zero_scalar(&Newcbehind,sizeof(O_CURVE));

	Oldcphys.curve = find_physical_curve_at_node(oldn,&Oldcphys.orient);
	Newcphys.curve = find_physical_curve_at_node(newn,&Newcphys.orient);
	
                /* Identify curves and components */

	find_propagation_orientation(fr,wave,oldn,newn,newp,&Oldcphys,dt,
				     &i_to_prop_dir,&Oldcahead,&Newcahead,
				     &Oldcbehind,&Newcbehind,&inc_side,
				     &propagation_side,&ahead_comp,
				     &propagation_comp);
	xgraphAtOldNode("oldn-nb.xg",oldn,Oldcphys,Oldcahead,Oldcbehind);
        xgraphAtOldNode("newn-nb-0.xg",newn,Newcphys,Newcahead,Newcbehind);

	copy_o_curve(&Newcaprop,&Newcahead);
	copy_o_curve(&Oldcaprop,&Oldcahead);
	copy_o_curve(&Newcbprop,&Newcbehind);
	copy_o_curve(&Oldcbprop,&Oldcbehind);

	     /* Identify new position of node */

	if (inc_side == propagation_side)
	{
	    printf("Incident side is positive\n");
	    ca = Newcaprop.curve;
	    status = D_extend_crossing_of_two_propagated_curves(
			    &Oldcaprop,&Newcaprop,&Oldcbehind,&Newcbehind,
			    &Oldcphys,&Newcphys,ahead_comp,propagation_comp,
			    pc,&crossbahead,&crossbphys,&tcr_ahead,
			    &tcr_phys,fr,wave,rp,dt,dt_frac,flag);
	    if (ca != Newcaprop.curve)
	    {
		/* Propagation direction reversed*/
		i_to_prop_dir = Opposite_ang_dir(i_to_prop_dir);
		inc_side = (curve_ang_oriented_l_to_r(i_to_prop_dir,Oldcaprop.orient)) ?
			    		NEGATIVE_SIDE : POSITIVE_SIDE;
	    }
	}
	else
	{
	    printf("Incident side is negative\n");
	    status = crossing_of_two_propagated_curves(
			&Oldcphys,&Newcphys,&Oldcaprop,&Newcaprop,
			pc,&crossbphys,&crossbahead,&tcr_phys,
			&tcr_ahead,fr,wave,rp,dt,dt_frac,flag);
	}
	xgraphAtOldNode("newn-nb-1.xg",newn,Newcphys,Newcahead,Newcbehind);
	if (status != GOOD_NODE)
        {
            printf("status != GOOD_NODE\n");
            clean_up(0);
        }
	
	    /* Modify the interface and assign the new states */

	status = modify_contact_node(oldn,newn,&Oldcphys,&Newcphys,&Oldcahead,
			       &Newcahead,&Oldcaprop,&Newcaprop,&Oldcbehind,
			       &Newcbehind,&Oldcbprop,&Newcbprop,pc,
		               crossbphys,crossbahead,i_to_prop_dir,tcr_phys,
		               tcr_ahead,rp,fr,wave,dt,dt_frac,flag);
	xgraphAtOldNode("newn-nb-2.xg",newn,Newcphys,Newcahead,Newcbehind);
	if (status == GOOD_NODE) 
	    printf("status == GOOD_NODE\n");

	propagation_status(newn) = PROPAGATED_NODE;
	return status;
}	/* end contact_node_propagate */

static void xgraphAtOldNode(
	const char *fname,
	NODE *n,
	O_CURVE cphys,
	O_CURVE cahead,
	O_CURVE cbehind)
{
	int count;
	BOND *b;
	FILE *file = fopen(fname,"w");

	fprintf(file,"Next\n");
        fprintf(file,"color=%s\n","red");
        fprintf(file,"shape 1\n");
        fprintf(file,"thickness=1.5\n");
	count = 0;
	if (cphys.orient == POSITIVE_ORIENTATION)
	{
	    for (b = cphys.curve->first; b != NULL && count < 4; b = b->next)
	    {
		fprintf(file,"%20.14f %20.14f\n",Coords(b->start)[0],Coords(b->start)[1]);
		count++;
	    }
	}
	else
	{
	    for (b = cphys.curve->last; b != NULL && count < 4; b = b->prev)
	    {
		fprintf(file,"%20.14f %20.14f\n",Coords(b->end)[0],Coords(b->end)[1]);
		count++;
	    }
	}
	fprintf(file,"Next\n");
        fprintf(file,"color=%s\n","blue");
        fprintf(file,"shape 1\n");
        fprintf(file,"thickness=1.5\n");
	count = 0;
	if (cahead.orient == POSITIVE_ORIENTATION)
	{
	    for (b = cahead.curve->first; b != NULL && count < 4; b = b->next)
	    {
		fprintf(file,"%20.14f %20.14f\n",Coords(b->start)[0],Coords(b->start)[1]);
		count++;
	    }
	}
	else
	{
	    for (b = cahead.curve->last; b != NULL && count < 4; b = b->prev)
	    {
		fprintf(file,"%20.14f %20.14f\n",Coords(b->end)[0],Coords(b->end)[1]);
		count++;
	    }
	}
	fprintf(file,"Next\n");
        fprintf(file,"color=%s\n","green");
        fprintf(file,"shape 1\n");
        fprintf(file,"thickness=1.5\n");
	count = 0;
	if (cbehind.orient == POSITIVE_ORIENTATION)
	{
	    for (b = cbehind.curve->first; b != NULL && count < 4; b = b->next)
	    {
		fprintf(file,"%20.14f %20.14f\n",Coords(b->start)[0],Coords(b->start)[1]);
		count++;
	    }
	}
	else
	{
	    for (b = cbehind.curve->last; b != NULL && count < 4; b = b->prev)
	    {
		fprintf(file,"%20.14f %20.14f\n",Coords(b->end)[0],Coords(b->end)[1]);
		count++;
	    }
	}
	fclose(file);
}	/* end xgraphAtOldNode */


/*
*			modify_contact_node():
*
*	Uses shift_node() and cut_curve() to modify the interface in the
*	neighborhood of a B_node().
*
*	Note: in the special case of crossing a fixed node
*	newcahead and newcbehind are not the correct curves.
*	(oldcahead is also wrong - simply crspd of newcahead).  Technically
*	speaking, the correct curves cannot be computed until the node has
*	been shifted.  In the normal case, all these curves are correct.
*/

static	int modify_contact_node(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		*oldcphys,
	O_CURVE		*newcphys,
	O_CURVE		*oldcahead,
	O_CURVE		*newcahead,
	O_CURVE		*oldcaprop,
	O_CURVE		*newcaprop,
	O_CURVE		*oldcbehind,
	O_CURVE		*newcbehind,
	O_CURVE		*oldcbprop,
	O_CURVE		*newcbprop,
	POINT		*pc,
	BOND		*crossbphys,
	BOND		*crossbahead,
	ANGLE_DIRECTION	i_to_prop_dir,
	double		tcr_phys,
	double		tcr_ahead,
	RPROBLEM	**rp,
	Front		*fr,
	POINTER		wave,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	int		status;
	SIDE		p_side_of_b, b_side_of_p;
	double		V[MAXD];
	static Locstate	left_st_a = NULL, right_st_a = NULL;
	static Locstate	left_st_p = NULL, right_st_p = NULL;
	static POINT	*p2 = NULL;

	if (debugging("contact_node"))
	    printf("Entering contact_node_propagate()\n");

	if (p2 == NULL) 
	{
	    p2 = Static_point(fr->interf);
	    if (fr->sizest)
	    {
	    	alloc_state(fr->interf, &left_st_a,fr->sizest);
	    	alloc_state(fr->interf,&right_st_a,fr->sizest);
	    	alloc_state(fr->interf, &left_st_p,fr->sizest);
	    	alloc_state(fr->interf,&right_st_p,fr->sizest);
	    }
	}

	/* Keep cross point on outer boundary  */
	/* IF outer boundary is appropriate.   */
	/* Modified to handle internal Neumann */
	/* boundaries. A more general approach */
	/* to deal with internal boundaries is */
	/* required.			       */

	if ((is_bdry(oldn) && is_bdry(oldcaprop->curve) &&
	     (to_next_node_only(flag) == YES))
				||
	    ((continue_past_fixed_node(flag) == YES) &&
	     is_bdry(oldcaprop->curve)))
	{
	    nearest_boundary_point(Coords(pc),Coords(pc),fr->rect_grid);
	}

			/* Interpolate physical states */

	left_state_along_bond(tcr_phys,crossbphys,newcphys->curve,left_st_p);
	right_state_along_bond(tcr_phys,crossbphys,newcphys->curve,right_st_p);

			/* Interpolate ahead states */

	left_state_along_bond(tcr_ahead,crossbahead,newcaprop->curve,left_st_a);
	right_state_along_bond(tcr_ahead,crossbahead,newcaprop->curve,
								right_st_a);

	/*  POINT_PROPAGATE CALL ALONE IS INAPPROPRIATE FOR SETTING   */
	/*    STATES FOR SOME FLOWS, CAUSING BOUNDARY PROPAGATION     */
	/* PROBLEMS HENCE, CODE FOLLOWING PT_PROP CALL HAS BEEN ADDED */

		/* Obtain behind states by propagating behind curve */

	point_propagate(fr,wave,oldn->posn,p2,
			Bond_at_node_of_o_curve(oldcbehind),
			oldcbehind->curve,dt,V);
	if (oldcbehind->orient != newcbehind->orient)
	    reverse_states_at_point(p2,fr);

		/* Obtain behind states on physical side */
		/*    from incident (physical) curve     */

	if (curve_ang_oriented_l_to_r(i_to_prop_dir,oldcbehind->orient))
	    p_side_of_b = POSITIVE_SIDE;
	else
	    p_side_of_b = NEGATIVE_SIDE;

	if (curve_ang_oriented_l_to_r(i_to_prop_dir,oldcphys->orient))
	    b_side_of_p = NEGATIVE_SIDE;
	else
	    b_side_of_p = POSITIVE_SIDE;

	if (fr->sizest)
	{
	    if ((wave_type(oldcbehind->curve) == SUBDOMAIN_BOUNDARY) &&
	        (wave_type(newcaprop->curve) == DIRICHLET_BOUNDARY))
	    {
	    /* The node is crossing the corner from a subdomain boundary
	     * to a Dirichlet boundary, so states on both sides of the
	     * new behind boundary need to be copied from behind the
	     * physical curve.
	     */

	        if (p_side_of_b == NEGATIVE_SIDE)
	        {
	    	    ft_assign(left_state(p2),left_st_p,fr->sizest);
	    	    obstacle_state(fr->interf,right_state(p2),fr->sizest);
	        }
	        else
	        {
	    	    obstacle_state(fr->interf,left_state(p2),fr->sizest);
	    	    ft_assign(right_state(p2),right_st_p,fr->sizest);
	        }
	    }
	    else if (p_side_of_b == b_side_of_p)
	    {
	        if (p_side_of_b == NEGATIVE_SIDE)
	    	    ft_assign(left_state(p2),left_st_p,fr->sizest);
	        else
	    	    ft_assign(right_state(p2),right_st_p,fr->sizest);
	    }
	    else
	    {
	        if (p_side_of_b == NEGATIVE_SIDE)
	    	    ft_assign(left_state(p2),right_st_p,fr->sizest);
	        else
	    	    ft_assign(right_state(p2),left_st_p,fr->sizest);
	    }
	}
	if (debugging("contact_node"))
	{
	    print_side("p_side_of_b ",p_side_of_b," ");
	    print_side("b_side_of_p ",b_side_of_p,"\n");
	}

		/* Modify boundary curves in vicinity of node */

	shift_node_past(pc,crossbahead,newcaprop->curve,newcaprop->orient,
			newcbehind->curve,newcbehind->orient,i_to_prop_dir,
			newn,fr,flag,left_st_a,right_st_a,left_state(p2),
			right_state(p2));

	cut_curve(pc,crossbphys,newcphys->curve,newcphys->orient,fr,
				left_st_p,right_st_p);

	if (continue_past_fixed_node(flag) == YES)
	{
	    newcbprop->curve = adjacent_curve(newcphys->curve,newcphys->orient,
					      Opposite_ang_dir(i_to_prop_dir),
					      &newcbprop->orient);
	    copy_o_curve(oldcbprop,oldcaprop);
	}

	if (fr->B_node_bifurcation)
	    status = (*fr->B_node_bifurcation)(fr,wave,oldcphys,newcphys,
			                       oldcahead,newcahead,oldcaprop,
					       newcaprop,oldcbehind,newcbehind,
					       oldcbprop,newcbprop,oldn->posn,
			                       left_st_p,right_st_p,
					       i_to_prop_dir,rp,dt,
					       dt_frac,flag);
	else
	    status = GOOD_NODE;

	if (debugging("contact_node"))
	    printf("Leaving contact_node_propagate()\n");
	return	status;
}		/*end modify_contact_node*/
