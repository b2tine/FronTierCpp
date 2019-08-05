#include "airfoil.h"


static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void rgbody_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);



EXPORT double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

EXPORT double getStatePhi(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->phi;
}	/* end getStatePhi */

EXPORT double getStateVort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort;
}	/* end getStateVort */

EXPORT double getStateXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[0];
}	/* end getStateXvel */

EXPORT double getStateYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[1];
}	/* end getStateYvel */

EXPORT double getStateZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[2];
}	/* end getStateZvel */

EXPORT double getStateXimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[0];
}	/* end getStateXimp */

EXPORT double getStateYimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[1];
}	/* end getStateYimp */

EXPORT double getStateZimp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->impulse[2];
}	/* end getStateZimp */

EXPORT double getStateMu(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->mu;
}	/* end getStateMu */

EXPORT double getStateTemp(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->temperature;
}	/* end getStateMTemp */


EXPORT void fluid_point_propagate(
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
}       /* fluid_point_propagate */

static void neumann_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
	F_FIELD *field = Fparams->field;
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
	if (fluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (fluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}

	setStateViscosity(Fparams,newst,comp);
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
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort,&oldst->vort);
	}
	FT_RecordMaxFrontSpeed(dim,0.0,NULL,Coords(newp),front);
        return;
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
        double speed;
        int i, dim = front->rect_grid->dim;
	STATE *newst = NULL;
	STATE *bstate;
	COMPONENT comp;
	F_PARAMS *Fparams = (F_PARAMS*)front->extra1;

	if (debugging("dirichlet_bdry"))
	{
	    (void) printf("Entering dirichlet_point_propagate()\n");
	    (void) print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
	}

	if (fluid_comp(negative_component(oldhs)))
	{
	    newst = (STATE*)left_state(newp);
	    comp = negative_component(oldhs);
	}
	else if (fluid_comp(positive_component(oldhs)))
	{
	    newst = (STATE*)right_state(newp);
	    comp = positive_component(oldhs);
	}
	setStateViscosity(Fparams,newst,comp);
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
	    newst->phi = bstate->phi;
            newst->vort = 0.0;

	    if (debugging("dirichlet_bdry"))
	    {
		(void) printf("Preset boundary state:\n");
		(void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
		(void) printf("Pressure: %f\n",newst->pres);
		(void) printf("Vorticity: %f\n",newst->vort);
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

static void contact_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
	F_FIELD *field = Fparams->field;
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
	setStateViscosity(Fparams,newst,negative_component(oldhs));
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	}
	newst = (STATE*)right_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	setStateViscosity(Fparams,newst,positive_component(oldhs));
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	}

	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
}	/* end contact_point_propagate */

static  void rgbody_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
	F_FIELD *field = Fparams->field;
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pre = field->pres;
	double *m_vor = field->vort;
	double *m_temp = field->temperature;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (fluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (fluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}
	setStateViscosity(Fparams,newst,comp);
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
	for (i = 0; i < dim; ++i) newst->vel[i] = vel[i];
	FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres,&oldst->pres);
	if (m_temp != NULL)
            FT_IntrpStateVarAtCoords(front,comp,p1,m_temp,
                        getStateTemp,&newst->temperature,&oldst->temperature);
	if (dim == 2)
	{
	    FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort,&oldst->vort);
	}
    
    if(!debugging("collision_off"))
    {
        /* copy newst to the other STATE; used in collision solver */
        if (fluid_comp(negative_component(oldhs)))
            std::copy(newst, newst+1, (STATE*)right_state(newp));
        else if (fluid_comp(positive_component(oldhs)))
            std::copy(newst, newst+1, (STATE*)left_state(newp));
    }
    return;
}	/* end rgbody_point_propagate */


EXPORT  void fluid_compute_force_and_torque(
        Front *fr,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
    return fluid_compute_force_and_torque3d(fr,hs,dt,force,torque);
}       /* end fluid_compute_force_and_torque */

static void fluid_compute_force_and_torque3d(
        Front *front,
        HYPER_SURF *hs,
        double dt,
        double *force,
        double *torque)
{
        RECT_GRID *gr = computational_grid(front->interf);
        double f[MAXD],rr[MAXD];
        double t[MAXD],tdir,pres;
        double area[MAXD],posn[MAXD];
        TRI *tri;
        boolean pos_side;
        int i,dim = gr->dim;
        F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
        double *gravity = Fparams->gravity;
        SURFACE *surface = Surface_of_hs(hs);
	CURVE **c;
        NODE *rg_string_nodes[10];
        int j, k, num = 0;
	NODE **n;
	BOND *b;
	double tri_cen[MAXD];

        if (debugging("rigid_body"))
	    (void) printf("Entering fluid_compute_force_and_torque3d()\n"); 
        if (fluid_comp(negative_component(surface)))
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
                    Coords((*n)->posn)[k] > gr->U[k])
                    break;
            }
            if (k != dim || (*n)->extra == NULL) continue;
	    node_out_curve_loop(*n,c)
	    {
                if (hsbdry_type(*c) == PASSIVE_HSBDRY)
		    break;
	    }
	    if (c == NULL || (*c) == NULL)
	    {
		node_in_curve_loop(*n,c)
		{
		    if (hsbdry_type(*c) == PASSIVE_HSBDRY)
			break;
		}
	    }
	    if (c == NULL || (*c) == NULL) continue;
	    b = (*c)->first;
	    if (wave_type(b->_btris[0]->surface) == MOVABLE_BODY_BOUNDARY)
		rg_string_nodes[num++] = *n;
	}
        for (j = 0; j < num; ++j)
        {
            POINT *p = rg_string_nodes[j]->posn;
	    if (!coords_in_subdomain(Coords(p),gr)) 
	    {
		continue;
	    }
            for (i = 0; i < dim; ++i)
            {
                force[i] += p->force[i];
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

	if (front->step > 5)
	{
            for (tri = first_tri(surface); !at_end_of_tri_list(tri,surface);
                        tri = tri->next)
            {
		for (i = 0; i < dim; ++i)
		{
		    tri_cen[i] = (Coords(Point_of_tri(tri)[0])[i] +
				  Coords(Point_of_tri(tri)[1])[i] +
				  Coords(Point_of_tri(tri)[2])[i])/3.0;
		}
	    	if (!coords_in_subdomain(tri_cen,gr)) 
		{
		    continue;
		}
                if (force_on_hse(Hyper_surf_element(tri),Hyper_surf(surface),gr,
                        &pres,area,posn,pos_side))
                {
                    for (i = 0; i < dim; ++i)
                    {
                        f[i] = pres*area[i];
                        force[i] += f[i];
                        rr[i] = posn[i] - rotation_center(surface)[i];
                    }
                    Cross3d(rr,f,t);
//		    tdir = Dot3d(t,(rotation_direction(hs)));
                    for (i = 0; i < dim; ++i)
                    {
//		        t[i] = tdir*rotation_direction(hs)[i];
                        torque[i] += t[i];
                    }
                }
            }
	}
         /* Add gravity to the total force */
        if (motion_type(surface) != ROTATION &&
	    motion_type(surface) != PRESET_ROTATION)
        {
            for (i = 0; i < dim; ++i)
                force[i] += gravity[i]*total_mass(surface)/num_clips(surface);
        }
        if (debugging("rigid_body"))
        {
            printf("In fluid_compute_force_and_torque3d()\n");
            printf("total_force = %f %f %f\n",force[0],force[1],force[2]);
            printf("torque = %f %f %f\n",torque[0],torque[1],torque[2]);
	    printf("# of rg_string_node in processor %d = %d\n", 
			pp_mynode(), num);
	    printf("number of clips = %d \n", num_clips(surface));
	
        }
        if (debugging("rigid_body"))
	    (void) printf("Leaving fluid_compute_force_and_torque3d()\n"); 
}       /* end fluid_compute_force_and_torque3d */
