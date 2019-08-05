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


#include <iFluid.h>
#include <crystal.h>
#include "subsurf.h"

	/*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);

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
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pre = iFparams->field->pres;
	double *m_vor = iFparams->field->vort;
	double nor[MAXD],tan[MAXD],p1[MAXD];
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
	FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;
	tan[0] = -nor[1]; 	tan[1] = nor[0];

	if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
	{
            double omega_dt,crds_com[MAXD];
            omega_dt = angular_velo(oldhs)*dt;
            for (i = 0; i < dim; ++i)
            {
                vel[i] = center_of_mass_velo(oldhs)[i];
                crds_com[i] = Coords(oldp)[i] +dt*vel[i] - 
			center_of_mass(oldhs)[i];
            }
            vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
            vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
	}
	else
	{
            for (i = 0; i < dim; ++i)
	    	vel[i] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
            Coords(newp)[i] = Coords(newp)[i] + dt*vel[i];
	    newst->vel[i] = vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
	}
	FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres,&oldst->pres);
	if (dim == 2)
        {
	    FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort,&oldst->vort);
	}
	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
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
	STATE *newst;
	STATE *bstate;
	FLOW_THROUGH_PARAMS ft_params;
	COMPONENT comp;

	if (debugging("dirichlet_bdry"))
	{
	    printf("Entering dirichlet_point_propagate()\n");
	    print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
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

	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
            for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] = bstate->vel[i];
		FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
                                        front);
	    }
	    speed = mag_vector(newst->vel,dim);
            FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
            newst->pres = bstate->pres;
            newst->vort = 0.0;

	    if (debugging("dirichlet_bdry"))
	    {
		printf("Preset boundary state:\n");
		print_general_vector("Velocity: ",newst->vel,dim,"\n");
		printf("Pressure: %f\n",newst->pres);
		printf("Vorticity: %f\n",newst->vort);
	    }
	}
	else if (boundary_state_function(oldhs))
	{
	    ft_params.oldp = oldp;
	    ft_params.comp = comp;
	    (*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)&ft_params,(POINTER)newst);	
	    /*
	    for (i = 0; i < dim; ++i)
                FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
                                        front);
            speed = mag_vector(newst->vel,dim);
            FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
	    */
	}
	if (debugging("dirichlet_bdry"))
	    printf("Leaving dirichlet_point_propagate()\n");
        return;
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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double *m_pre = iFparams->field->pres;
	double *m_vor = iFparams->field->vort;
	double *p0;
	STATE *oldst,*newst;
	POINTER sl,sr;
	double pres,vort;

        (*front->vfunc)(front->vparams,front,oldp,oldhse,oldhs,vel);
        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(newp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	oldst = (STATE*)sl;
	p0 = Coords(newp);
	FT_IntrpStateVarAtCoords(front,-1,p0,m_pre,getStatePres,&pres,
				&oldst->pres);
	FT_IntrpStateVarAtCoords(front,-1,p0,m_vor,getStateVort,&vort,
				&oldst->vort);

	newst = (STATE*)left_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	    newst->vel[i] = vel[i];
	newst = (STATE*)right_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	    newst->vel[i] = vel[i];

	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
}	/* end contact_point_propagate */

void read_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	static STATE state;
	HYPER_SURF *hs;
	int i_surf;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		i_surf = 2*i;
		if (rect_boundary_type(front->interf,i,0) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
						i,0);
		sprintf(msg,"For lower boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
		    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs,i_surf);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	    if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		i_surf = 2*i + 1;
		if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
						i,1);
		sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
		    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs,i_surf);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_dirichlet_bdry_data */

void read_ss_dirichlet_bdry_data(
	char *inname,
        Front *front,
        F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,k,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	STATE state;
	HYPER_SURF *hs;
	int i_surf;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		i_surf = 2*i;
	        if (rect_boundary_type(front->interf,i,0) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
					i,0);
		sprintf(msg,"For lower boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf",&state.vel[k]);
			(void) printf("%f ",state.vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state.pres);
		    (void) printf("%f\n",state.pres);
		    CursorAfterString(infile,"Enter solute concentration:");
                    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
				(POINTER)&state,hs,i_surf);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_InsertDirichletBoundary(front,ss_flowThroughBoundaryState,
				"flowThroughBoundaryState",NULL,NULL,hs,i_surf);
		    break;
		}
	    }
            if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		i_surf = 2*i + 1;
                if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)                    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
                                                i,1);
		sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf ",&state.vel[k]);
			(void) printf("%f ",state.vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state.pres);
		    (void) printf("%f\n",state.pres);
		    CursorAfterString(infile,"Enter solute concentration:");
                    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_InsertDirichletBoundary(front,NULL,NULL,NULL,
				(POINTER)&state,hs,i_surf);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_InsertDirichletBoundary(front,ss_flowThroughBoundaryState,
				"flowThroughBoundaryState",NULL,NULL,hs,i_surf);
		    break;
		}
	    }
	}
	fclose(infile);
}	/* end read_ss_dirichlet_bdry_data */

void init_fluid_state_func(
	Front *front,
        Incompress_Solver_Smooth_Basis *l_cartesian)
{
	l_cartesian->getInitialState = zero_state;
}	/* end init_fluid_state_func */

static void zero_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
        int i;
        for (i = 0; i < dim; ++i)
            field->vel[i][index] = 0.0;
        field->pres[index] = 0.0;
}       /* end zero_state */

static double (*getStateVel[MAXD])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

extern void ss_flowThroughBoundaryState(
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
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
	IF_FIELD *iF_field = iFparams->field;
	CRT_FIELD *cR_field = cRparams->field;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double vort[3];		/* vorticity stencil */
	double pres[3];		/* pressure stencil */
	double solu[3];		/* solute stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_vort;		/* vort flux */
	double f_pres;		/* pressure flux */
	double f_solu;		/* solute flux */
	double dn,dt = front->dt;
	STATE *newst = (STATE*)state;
	STATE  **sts;
	int i,j,dim = front->rect_grid->dim;
	int nrad = 2;
	
	if (debugging("flow_through"))
	    printf("Entering ss_flowThroughBoundaryState()\n");

	tsten = FrontGetTanStencils(front,oldp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = tsten[0]->dir[i];
	dn = FT_GridSizeInDir(dir,front);

	if (comp == negative_component(hs))  
	    sts = (STATE**)tsten[0]->leftst;
	else 
	    sts = (STATE**)tsten[0]->rightst;

	if (debugging("flow_through"))
	{
	    printf("Ambient component: %d\n",comp);
	    printf("hs = %p  oldp->hs = %p\n",(POINTER)hs,(POINTER)oldp->hs);
	    printf("Time step = %f  Tangential grid size = %f\n",dt,dn);
	    printf("Tangential direction: ");
	    for (j = 0; j < dim; ++j)
		printf("%f ",tsten[0]->dir[j]);
	    printf("\n");
	    printf("Tan_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    printf("Left points:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",Coords(tsten[0]->p[-i])[j]);
		printf("\n");
	    }
	    printf("Right points:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",Coords(tsten[0]->p[i])[j]);
		printf("\n");
	    }
	}

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 3; ++j)
	{
	    vort[j] = sts[j-1]->vort;
	    pres[j] = sts[j-1]->pres;
	    solu[j] = sts[j-1]->solute;
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
	f_solu = linear_flux(u[1],solu[0],solu[1],solu[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] = sts[0]->vel[i] - dt/dn*(
		f_u*dir[i] + f_v[i]) ;
	newst->vort = sts[0]->vort - dt/dn*f_vort;
	newst->pres = sts[0]->pres - dt/dn*f_pres;
	newst->solute = sts[0]->solute - dt/dn*f_solu;
	
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
	    vort[j] = sts[0]->vort;
	    pres[j] = sts[0]->pres;
	    solu[j] = sts[0]->solute;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += sts[0]->vel[i]*dir[i];
		v[j][i] = sts[0]->vel[i]*(1.0 - dir[i]);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			iF_field->vel[i],getStateVel[i],&vtmp,&sts[0]->vel[i]);
	    u[2] += vtmp*dir[i];
	    v[2][i] = vtmp*(1.0 - dir[i]);
	}
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			iF_field->vort,getStateVort,&vort[2],&sts[1]->vort);
	}
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],iF_field->pres,
                            getStatePres,&pres[2],&sts[1]->pres);
	FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],cR_field->solute,
                            getStateSolute,&solu[2],&sts[1]->solute);

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_solu = linear_flux(u[1],solu[0],solu[1],solu[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;
	newst->solute += - dt/dn*f_solu;
	if (debugging("dirichlet_bdry"))
	{
	    printf("flow through boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
	    printf("Vorticity: %f\n",newst->vort);
	    printf("Solute: %f\n",newst->solute);
	}
}       /* end ss_flowThroughBoundaryState */

extern void read_fluid_params(
	char *inname,
	IF_PARAMS *iFparams)
{
	char string[100];
	FILE *infile = fopen(inname,"r");

	/* defaults numerical schemes */
        iFparams->num_scheme.projc_method = SIMPLE;
        iFparams->num_scheme.advec_method = WENO;
        iFparams->num_scheme.ellip_method = SIMPLE_ELLIP;

	CursorAfterString(infile,"Enter density and viscosity of the fluid:");
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	iFparams->m_comp1 = CRYSTAL_COMP;
	iFparams->m_comp2 = LIQUID_COMP2;
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);
	iFparams->num_scheme.projc_method = ERROR_PROJC_SCHEME;
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
            iFparams->num_scheme.projc_method = BELL_COLELLA;
            break;
        case 'K':
        case 'k':
            iFparams->num_scheme.projc_method = KIM_MOIN;
            break;
        case 'P':
        case 'p':
            iFparams->num_scheme.projc_method = PEROT_BOTELLA;
        }
	assert(iFparams->num_scheme.projc_method != ERROR_PROJC_SCHEME);
	(void) printf("The default advection order is 4\n");
        iFparams->adv_order = 4;
        if (CursorAfterStringOpt(infile,"Enter advection order:"))
        {
            fscanf(infile,"%d",&iFparams->adv_order);
            (void) printf("%d\n",iFparams->adv_order);
        }
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
            default:
                printf("Elliptic Method Not Implemented\n");
                clean_up(1);
            }
        }
	iFparams->ub_speed = HUGE;
	iFparams->skip_neumann_solver = YES;

	fclose(infile);
}	/* end read_fluid_params */
