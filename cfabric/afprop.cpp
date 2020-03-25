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

#include "fabric.h"
#include "collid.h"

#include "cFluid.h"


static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,getStateZvel};

static void string_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void mono_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void gore_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void passive_curve_propagation(Front*,POINTER,CURVE*,CURVE*,double);
static void gore_point_propagate(Front*,POINTER,POINT*,POINT*,BOND*,double);
static void load_node_propagate(Front*,NODE*,NODE*,double);
static void rg_string_node_propagate(Front*,NODE*,NODE*,double);

static void fourth_order_elastic_set_propagate2d(Front*,double);
static void fourth_order_elastic_set_propagate3d(Front*,double);
static void setCollisionFreePoints3d(INTERFACE*);

static void coating_mono_hyper_surf3d(Front*);


extern void fourth_order_elastic_set_propagate(Front* fr, double fr_dt)
{
        switch (fr->rect_grid->dim)
        {
        case 2:
            return fourth_order_elastic_set_propagate2d(fr,fr_dt);
        case 3:
            return fourth_order_elastic_set_propagate3d(fr,fr_dt);
        }
}       /* end fourth_order_elastic_set_propagate */

static void fourth_order_elastic_set_propagate2d(Front* fr, double fr_dt)
{
        CURVE **c,*elastic_curve;
        intfc_curve_loop(fr->interf,c)
        {
            if (wave_type(*c) == ELASTIC_BOUNDARY ||
                wave_type(*c) == ELASTIC_STRING)
            {
                elastic_curve = *c;
                break;
            }
        }
        fourth_order_elastic_curve_propagate(fr,fr,fr->interf,elastic_curve,
                        elastic_curve,fr_dt);
}       /* end fourth_order_elastic_set_propagate2d */

static void fourth_order_elastic_set_propagate3d(Front* fr, double fr_dt)
{
	static ELASTIC_SET geom_set;
	static int size = 0,owner_size,client_size;
	static int *client_size_old, *client_size_new;
        AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
        int i,j,k,n_sub;
        double dt;
        static SPRING_VERTEX *sv;
        static boolean first = YES;
        static GLOBAL_POINT **point_set;
        static GLOBAL_POINT *point_set_store;
	static GLOBAL_POINT **client_point_set_store;
        int dim = FT_Dimension();
        long max_point_gindex = fr->interf->max_point_gindex;
	int owner[MAXD];
	int owner_id = af_params->node_id[0];
        int myid = pp_mynode();
	int gindex;
        INTERFACE *elastic_intfc = NULL;
	double *L = fr->rect_grid->L;
	double *U = fr->rect_grid->U;
	double *h = fr->rect_grid->h;
	double client_L[MAXD],client_U[MAXD];
	static boolean first_break_strings = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;
        
	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate3d()\n");

    geom_set.front = fr;

	if (first_break_strings && break_strings_num > 0 &&
	    break_strings_time >= 0.0 && 
	    fr->time + fr->dt >= break_strings_time)
	{
	    printf("Some strings break! Count and set spring vertex again.\n");
	    first_break_strings = NO;
	    first = YES;
	}

	if (first)
    {
        set_elastic_params(&geom_set,fr_dt);
        if (debugging("step_size"))
            print_elastic_params(geom_set);
    }

    if (fr_dt > geom_set.dt_tol)
    {
        n_sub = (int)(fr_dt/geom_set.dt_tol);
        dt = fr_dt/n_sub;
    }
    else
    {
        n_sub = af_params->n_sub;
        dt = fr_dt/n_sub;
    }
    printf("fr_dt = %f  dt = %f  n_sub = %d\n",fr_dt,dt,n_sub);

	if (first)
	{
            owner[0] = 0;
            owner[1] = 0;
            owner[2] = 0;
	    if (point_set != NULL)
		FT_FreeThese(1, point_set);
	    FT_VectorMemoryAlloc((POINTER*)&point_set,max_point_gindex,
					sizeof(GLOBAL_POINT*));
	    for (i = 0; i < max_point_gindex; ++i)
		point_set[i] = NULL;

	    if (pp_numnodes() > 1)
	    {
                elastic_intfc =
                    FT_CollectHypersurfFromSubdomains(fr,owner,ELASTIC_BOUNDARY);
                collectNodeExtra(fr,elastic_intfc,owner_id);
	    }
	    else
            elastic_intfc = fr->interf;
	    
            start_clock("set_data");
	    if (myid == owner_id)
            {
		if (client_size_old != NULL)
		    FT_FreeThese(3, client_size_old, client_size_new, 
					client_point_set_store);
		FT_VectorMemoryAlloc((POINTER*)&client_size_old,pp_numnodes(),
                                        sizeof(int));
                FT_VectorMemoryAlloc((POINTER*)&client_size_new,pp_numnodes(),
                                        sizeof(int));
                FT_VectorMemoryAlloc((POINTER*)&client_point_set_store,
                                        pp_numnodes(),sizeof(GLOBAL_POINT*));
                for (i = 0; i < pp_numnodes(); i++)
                    client_size_old[i] = client_size_new[i] = 0;

		assembleParachuteSet(elastic_intfc,&geom_set);
		owner_size = geom_set.num_verts;
		if (point_set_store != NULL) 
		    FT_FreeThese(2,point_set_store, sv);
		FT_VectorMemoryAlloc((POINTER*)&point_set_store,owner_size,
                                        sizeof(GLOBAL_POINT));
                FT_VectorMemoryAlloc((POINTER*)&sv,owner_size,
                                        sizeof(SPRING_VERTEX));
		link_point_set(&geom_set,point_set,point_set_store);
	    	count_vertex_neighbors(&geom_set,sv);
	    	set_spring_vertex_memory(sv,owner_size);
	    	set_vertex_neighbors(&geom_set,sv,point_set);
		
                if (elastic_intfc != fr->interf)
                    delete_interface(elastic_intfc);
	    }
	    stop_clock("set_data");
	    first = NO;
	}

	elastic_intfc = fr->interf;
	assembleParachuteSet(elastic_intfc,&geom_set);
	
    if (myid != owner_id)
	{
	    client_size = geom_set.num_verts;
	    if (size < client_size)
	    {
	    	size = client_size;
	    	if (point_set_store != NULL)
		{
		    FT_FreeThese(2,point_set_store,sv);
		}
	    	FT_VectorMemoryAlloc((POINTER*)&point_set_store,size,
                                        sizeof(GLOBAL_POINT));
                FT_VectorMemoryAlloc((POINTER*)&sv,size,sizeof(SPRING_VERTEX));
	    }
	    for (i = 0; i < max_point_gindex; ++i)
                point_set[i] = NULL;
	    link_point_set(&geom_set,point_set,point_set_store);
	    count_vertex_neighbors(&geom_set,sv);
	    set_spring_vertex_memory(sv,client_size);
	    set_vertex_neighbors(&geom_set,sv,point_set);
	    get_point_set_from(&geom_set,point_set);
	    pp_send(5,L,MAXD*sizeof(double),owner_id);
	    pp_send(6,U,MAXD*sizeof(double),owner_id);
	    pp_send(1,&(client_size),sizeof(int),owner_id);
            pp_send(2,point_set_store,client_size*sizeof(GLOBAL_POINT),
					owner_id);
	}
	else
    {
	    size = owner_size;
    }

    CollisionSolver3d* collision_solver;
	if (!debugging("collision_off"))
    {
        collision_solver = new CollisionSolver3d();
        printf("COLLISION DETECTION ON\n");
    }
    else
        printf("COLLISION DETECTION OFF\n");

    if (myid == owner_id)
	{
        if (!debugging("collision_off"))
        {
            setCollisionFreePoints3d(fr->interf);

            collision_solver->assembleFromInterface(fr->interf,fr->dt);
            collision_solver->recordOriginalPosition();
        
            collision_solver->setRestitutionCoef(1.0);
            collision_solver->setVolumeDiff(af_params->vol_diff);
            
            collision_solver->setFabricRoundingTolerance(af_params->fabric_eps);
            collision_solver->setFabricThickness(af_params->fabric_thickness);
            collision_solver->setFabricFrictionConstant(af_params->mu_s);
            collision_solver->setFabricSpringConstant(af_params->ks); 
            collision_solver->setFabricPointMass(af_params->m_s);

            collision_solver->setStringRoundingTolerance(af_params->string_eps);
            collision_solver->setStringThickness(af_params->string_thickness);
            collision_solver->setStringFrictionConstant(af_params->mu_l);
            collision_solver->setStringSpringConstant(af_params->kl); 
            collision_solver->setStringPointMass(af_params->m_l);

            collision_solver->setStrainLimit(af_params->strain_limit);
            collision_solver->setStrainRateLimit(af_params->strainrate_limit);

            collision_solver->gpoints = fr->gpoints;
            collision_solver->gtris = fr->gtris;
        }

        //write to GLOBAL_POINT** point_set
        get_point_set_from(&geom_set,point_set);

	    for (i = 0; i < pp_numnodes(); i++)
	    {
		if (i == myid) continue;
		pp_recv(5,i,client_L,MAXD*sizeof(double));
		pp_recv(6,i,client_U,MAXD*sizeof(double));
		pp_recv(1,i,client_size_new+i,sizeof(int));
		if (client_size_new[i] > client_size_old[i])
		{
		    client_size_old[i] = client_size_new[i];
		    if (client_point_set_store[i] != NULL)
		    	FT_FreeThese(1,client_point_set_store[i]);
	    	    FT_VectorMemoryAlloc((POINTER*)&client_point_set_store[i],
				client_size_new[i], sizeof(GLOBAL_POINT));
		}
		pp_recv(2,i,client_point_set_store[i],
		    client_size_new[i]*sizeof(GLOBAL_POINT));
		copy_from_client_point_set(point_set,client_point_set_store[i],
				client_size_new[i],client_L,client_U);
	    } 

	    start_clock("spring_model");
#if defined(__GPU__)
            if (af_params->use_gpu)
            {
            	if (debugging("trace"))
                    (void) printf("Enter gpu_spring_solver()\n");
                gpu_spring_solver(sv,dim,size,n_sub,dt);
                if (debugging("trace"))
                    (void) printf("Left gpu_spring_solver()\n");
            }
            else
#endif
                generic_spring_solver(sv,dim,size,n_sub,dt);
	    stop_clock("spring_model");

	    for (i = 0; i < pp_numnodes(); i++)
            {
                if (i == myid) continue;
                copy_to_client_point_set(point_set,
                            client_point_set_store[i], client_size_new[i]);
                pp_send(3,client_point_set_store[i],
                            client_size_new[i]*sizeof(GLOBAL_POINT),i);
            }
	}

    if (myid != owner_id)
    {
        pp_recv(3,owner_id,point_set_store,
        client_size*sizeof(GLOBAL_POINT));
    }

	/* Owner send and patch point_set_store from other processors */	
	put_point_set_to(&geom_set,point_set);
	
    /* Calculate the real force on load_node and rg_string_node */
	setSpecialNodeForce(fr, geom_set.kl);

	set_vertex_impulse(&geom_set,point_set);
	set_geomset_velocity(&geom_set,point_set);
	compute_center_of_mass_velo(&geom_set);

	if (!debugging("collision_off"))
    {
        if (myid == owner_id)
        {
            if (FT_Dimension() == 3)
                collision_solver->resolveCollision();
        }
        setSpecialNodeForce(fr,geom_set.kl);

        delete collision_solver;
    }
    
	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate3d()\n");
}	/* end fourth_order_elastic_set_propagate3d() */

static void setCollisionFreePoints3d(INTERFACE* intfc)
{
    POINT *p;
    HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
    SURFACE* surf;
    
    if (intfc->dim == 2) {
        printf("ERROR dim = %d\n",intfc->dim);
        clean_up(ERROR);
    }

    next_point(intfc,NULL,NULL,NULL);
    while(next_point(intfc,&p,&hse,&hs))
    {
        STATE* sl = (STATE*)left_state(p);
        sl->is_fixed = false;
        sl->is_movableRG = false;
        
        if ((surf = Surface_of_hs(hs)) &&
                (is_registered_point(surf,p) ||
                 wave_type(hs) == NEUMANN_BOUNDARY))
        {
            sl->is_fixed = true;
        }
    
        if ((surf = Surface_of_hs(hs)) &&
                (wave_type(hs) == MOVABLE_BODY_BOUNDARY))
        {
            sl->is_movableRG = true;
        }
    }

    CURVE **c;
    BOND* b;
    intfc_curve_loop(intfc,c)
    {
        if (hsbdry_type(*c) != FIXED_HSBDRY)
            continue;

        for (b = (*c)->first; b != (*c)->last; b = b->next)
        {
            STATE* sl = (STATE*)left_state(b->end);
            sl->is_fixed = true;
        }
    }

    NODE** n;
    intfc_node_loop(intfc,n)
    {
        STATE* sl = (STATE*)left_state((*n)->posn);
        sl->is_fixed = false;
        
        AF_NODE_EXTRA* extra = (AF_NODE_EXTRA*)(*n)->extra;
        if (extra && extra->af_node_type == PRESET_NODE)
        {
            sl->is_fixed = true;
        }
        else if ((*n)->hsb && is_fixed_node(*n))
        {
            sl->is_fixed = true;
        }
    }
}       /* setCollisionFreePoints3d() */

extern void elastic_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
        
	if (af_params->no_fluid)
	{
	    fourth_order_point_propagate(front,wave,oldp,newp,oldhse,
				oldhs,dt,V);
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	    return;
	}
	
    EQN_PARAMS *cFparams = (EQN_PARAMS*)front->extra1;
	double *vort = cFparams->vort;
	double **vel = cFparams->vel;
	double *pres = cFparams->pres;

    /*
	F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
    F_FIELD *field = Fparams->field;
	double *vort = field->vort;
	double **vel = field->vel;
	double *pres = field->pres;
    */
	
    COMPONENT base_comp = positive_component(oldhs);
	double pp[MAXD],pm[MAXD],nor[MAXD],h;
	double area_dens = af_params->area_dens;
	double left_nor_speed,right_nor_speed;
	double dv[MAXD];

    int dim = front->rect_grid->dim;
    STATE *newsl,*newsr;
	STATE *sl,*sr;

	sl = (STATE*)left_state(oldp);
	sr = (STATE*)right_state(oldp);
	//FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	newsl = (STATE*)left_state(newp);
	newsr = (STATE*)right_state(newp);

	FT_NormalAtPoint(oldp,front,nor,NO_COMP);
	h = FT_GridSizeInDir(nor,front);
	for (int i = 0; i < dim; ++i)
	{
	    pm[i] = Coords(oldp)[i] - h*nor[i];
	    pp[i] = Coords(oldp)[i] + h*nor[i];
	}

    //Interpolate grid pressure to fabric points
    FT_IntrpStateVarAtCoords(front,base_comp-1,pm,pres,
    getStatePres,&newsl->pres,&sl->pres);
    FT_IntrpStateVarAtCoords(front,base_comp+1,pp,pres,
    getStatePres,&newsr->pres,&sr->pres);

	/* Impulse is incremented by the fluid pressure force */
	for (int i = 0; i < dim; ++i)
	{
	    dv[i] = 0.0;

	    if (debugging("rigid_canopy"))
            dv[i] = 0.0;
        else if (front->step > 5)
            dv[i] = (sl->pres - sr->pres)*nor[i]/area_dens;
            
        newsr->fluid_accel[i] = newsl->fluid_accel[i] = dv[i];
        newsr->other_accel[i] = newsl->other_accel[i] = 0.0;
        newsr->impulse[i] = newsl->impulse[i] = sl->impulse[i];
        newsr->vel[i] = newsl->vel[i] = sl->vel[i];
	}
}       /* elastic_point_propagate */


extern void airfoil_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        if (wave_type(oldhs) == ELASTIC_BOUNDARY ||
            wave_type(oldhs) == ELASTIC_STRING)
        {
            return elastic_point_propagate(
                    front,wave,oldp,newp,oldhse,oldhs,dt,V);
        }
        else
        {
            return cFluid_point_propagate(
                    front,wave,oldp,newp,oldhse,oldhs,dt,V);
            //return ifluid_point_propagate(
              //      front,wave,oldp,newp,oldhse,oldhs,dt,V);
        }
}       /* airfoil_point_propagate */

extern void airfoil_curve_propagate(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	int dim = front->rect_grid->dim;
        if (dim == 2) return;
        switch (hsbdry_type(oldc))
        {
        case STRING_HSBDRY:
            return string_curve_propagation(front,wave,oldc,newc,dt);
        case MONO_COMP_HSBDRY:
            return mono_curve_propagation(front,wave,oldc,newc,dt);
        case GORE_HSBDRY:
            return gore_curve_propagation(front,wave,oldc,newc,dt);
        case PASSIVE_HSBDRY:
            return passive_curve_propagation(front,wave,oldc,newc,dt);
        default:
            return;
        }
}	/* end airfoil_curve_propagate */

static void string_curve_propagation(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;

	if (!is_load_node(oldc->start))
	{
	    oldp = oldc->start->posn;
	    newp = newc->start->posn;
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	}

	if (!is_load_node(oldc->end))
	{
	    oldp = oldc->end->posn;
	    newp = newc->end->posn;
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	}

	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	}
}	/* end string_curve_propagation */

static void gore_curve_propagation(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;

	if (debugging("interact_curve"))
	{
	    (void) printf("Entering gore_curve_propagation()\n");
	}
	oldp = oldc->start->posn;
	newp = newc->start->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	oldp = oldc->end->posn;
	newp = newc->end->posn;
	ft_assign(left_state(newp),left_state(oldp),front->sizest);
	ft_assign(right_state(newp),right_state(oldp),front->sizest);

	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    gore_point_propagate(front,wave,oldp,newp,oldb,dt);
	}
	if (debugging("interact_curve"))
	{
	    (void) printf("Leaving gore_curve_propagation()\n");
	}
}	/* end gore_curve_propagation */

static void gore_point_propagate(
	Front *front,
        POINTER wave,
	POINT *oldp,
	POINT *newp,
	BOND *oldb,
	double dt)
{
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *oldhse;
	HYPER_SURF         *oldhs;
	STATE *sl,*sr,*newsl,*newsr;
	double mag_nor,branch_nor[MAXD],nor[MAXD];
	double pm[MAXD],pp[MAXD],h;
	 
    /*
    F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
	F_FIELD *field = Fparams->field;
	
    double **vel = field->vel;
    double *pres = field->pres;
    */

    EQN_PARAMS *cFparams = (EQN_PARAMS*)front->extra1;
    double **vel = cFparams->vel;
    double *pres = cFparams->pres;
	
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    double area_dens = af_params->area_dens;

	double left_nor_speed,right_nor_speed,dv;
	COMPONENT base_comp;
	double V[MAXD];
	int i;

	if (af_params->no_fluid)
	{
	    for (btris = Btris(oldb); btris && *btris; ++btris)
	    {
	    	oldhse = Hyper_surf_element((*btris)->tri);
	    	oldhs = Hyper_surf((*btris)->surface);
	    }
	    fourth_order_point_propagate(front,wave,oldp,newp,oldhse,
				oldhs,dt,V);
	    ft_assign(left_state(newp),left_state(oldp),front->sizest);
	    ft_assign(right_state(newp),right_state(oldp),front->sizest);
	    return;
	}
	sl = (STATE*)left_state(oldp);		
	sr = (STATE*)right_state(oldp);
	newsl = (STATE*)left_state(newp);	
	newsr = (STATE*)right_state(newp);

	for (i = 0; i < 3; ++i) nor[i] = 0.0;
	for (btris = Btris(oldb); btris && *btris; ++btris)
	{
	    oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    oldp->hs = oldhs = Hyper_surf((*btris)->surface);
	    FT_NormalAtPoint(oldp,front,branch_nor,NO_COMP);
	    base_comp = positive_component(oldhs);
	    for (i = 0; i < 3; ++i) nor[i] += branch_nor[i];
	}
	mag_nor = Mag3d(nor);
	for (i = 0; i < 3; ++i) nor[i] /= mag_nor;
	h = FT_GridSizeInDir(nor,front);
	for (i = 0; i < 3; ++i)
	{
	    pm[i] = Coords(oldp)[i] - h*nor[i];
            pp[i] = Coords(oldp)[i] + h*nor[i];
	}
	FT_IntrpStateVarAtCoords(front,base_comp-1,pm,pres,
                        getStatePres,&newsl->pres,&sl->pres);
        FT_IntrpStateVarAtCoords(front,base_comp+1,pp,pres,
                        getStatePres,&newsr->pres,&sr->pres);
	for (i = 0; i < 3; ++i)
        {
            FT_IntrpStateVarAtCoords(front,base_comp-1,pm,vel[i],
                        getStateVel[i],&newsl->vel[i],&sl->vel[i]);
            FT_IntrpStateVarAtCoords(front,base_comp+1,pp,vel[i],
                        getStateVel[i],&newsr->vel[i],&sr->vel[i]);
        }
        left_nor_speed = Dot3d(newsl->vel,nor);
        right_nor_speed = Dot3d(newsr->vel,nor);
	for (i = 0; i < 3; ++i)
        {
            newsl->vel[i] -= left_nor_speed*nor[i];
            newsr->vel[i] -= right_nor_speed*nor[i];
        }
	/* Impulse is incremented by the fluid pressure force */
        for (i = 0; i < 3; ++i)
        {
	    dv = 0.0;

	    if (front->step > 5)
		dv = (sl->pres - sr->pres)*nor[i]/area_dens;
	    if (debugging("rigid_canopy"))
	    	dv = 0.0;
	    newsr->fluid_accel[i] = newsl->fluid_accel[i] = dv;
	    newsr->other_accel[i] = newsl->other_accel[i] = 0.0;
	    newsr->impulse[i] = newsl->impulse[i] = sl->impulse[i];
	}
}	/* end gore_point_propagate */

static void mono_curve_propagation(
        Front *front,
        POINTER wave,
	CURVE *oldc,
	CURVE *newc,
        double dt)
{
	BOND *oldb,*newb;
	POINT *oldp,*newp;
	double V[MAXD];
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *oldhse;
	HYPER_SURF         *oldhs;

	if (debugging("interact_curve"))
	{
	    (void) printf("Entering mono_curve_propagation()\n");
	}

	/*
	oldb = oldc->first;
	newb = newc->first;
	oldp = oldb->start;
	newp = newb->start;
	for (btris = Btris(oldb); btris && *btris; ++btris)
	{
	    oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    oldp->hs = oldhs = Hyper_surf((*btris)->surface);
	    elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	}
	*/
	for (oldb = oldc->first, newb = newc->first; oldb != NULL;
		oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    for (btris = Btris(oldb); btris && *btris; ++btris)
	    {
	    	oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    	oldp->hs = oldhs = Hyper_surf((*btris)->surface);
		elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	    }
	}
	if (debugging("interact_curve"))
	{
	    (void) printf("Leaving mono_curve_propagation()\n");
	}
}	/* end mono_curve_propagation */

static void load_node_propagate(
	Front *front,
	NODE *oldn,
	NODE *newn,
	double dt)
{
	F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	POINT *oldp,*newp;
	double *g = af_params->gravity;
	double f[MAXD],accel[MAXD];
	double kl = af_params->kl;
	double mass = af_params->payload;
	CURVE **c;
	STATE *sl,*sr,*newsl,*newsr;
	double vec[MAXD],vec_mag;
	BOND *b;
	int i,dim = FT_Dimension();

	if (!is_load_node(oldn)) return;
	oldp = oldn->posn;
	newp = newn->posn;
	sl = (STATE*)left_state(oldp);
	sr = (STATE*)right_state(oldp);
	newsl = (STATE*)left_state(newp);
	newsr = (STATE*)right_state(newp);

	if (debugging("trace"))
	{
	    (void)printf("\nEntering load_node_propagate()\n");
	}
	for (i = 0; i < dim; ++i)
	    f[i] = 0.0;
	node_out_curve_loop(oldn,c)
	{
	    b = (*c)->first;
	    for (i = 0; i < dim; ++i)
	    {
		vec[i] = Coords(b->end)[i] - Coords(b->start)[i];
		vec[i] /= bond_length(b);
		f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
	    }
	}
	node_in_curve_loop(oldn,c)
	{
	    b = (*c)->last;
	    for (i = 0; i < dim; ++i)
	    {
		vec[i] = Coords(b->start)[i] - Coords(b->end)[i];
		vec[i] /= bond_length(b);
		f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    accel[i] = f[i]/mass;
	    newsl->fluid_accel[i] = newsr->fluid_accel[i] = 0.0;
	    newsr->other_accel[i] = newsl->other_accel[i] = accel[i];
	    newsl->impulse[i] = newsr->impulse[i] = sl->impulse[i];
	    newsl->vel[i] = newsr->vel[i] = sl->vel[i] + 
				(accel[i] + g[i]) * dt;
	}
	node_out_curve_loop(newn,c)
	{
	    b = (*c)->first;
	    set_bond_length(b,dim);
	}
	node_in_curve_loop(newn,c)
	{
	    b = (*c)->last;
	    set_bond_length(b,dim);
	}
	
	if (debugging("trace"))
	{
	    (void)printf("accel = %f %f %f\n",accel[0],accel[1],
				accel[2]);
	    (void)printf("Leaving load_node_propagate()\n\n");
	}
}	/* end load_node_propagate */

static void rg_string_node_propagate(
        Front *front,
        NODE *oldn,
        NODE *newn,
        double dt)
{
	F_PARAMS *Fparams = (F_PARAMS*)front->extra1;
        AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	RECT_GRID *gr = computational_grid(front->interf);
        POINT *oldp,*newp;
        double *g = af_params->gravity;
        double f[MAXD],accel[MAXD];
        double kl = af_params->kl;
        double mass = af_params->payload;
        CURVE **c;
        STATE *sl,*sr,*newsl,*newsr;
        double vec[MAXD],vec_mag;
        BOND *b;
        int i,dim = FT_Dimension();
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	POINTER wave;
	double V[MAXD];

        if (!is_rg_string_node(oldn)) return;
	for (i = 0; i < dim; ++i)
	{
	    if (Coords(oldn->posn)[i] <= gr->L[i] || 
		Coords(oldn->posn)[i] > gr->U[i])
		break;
	}
	if (i != dim || oldn->extra == NULL) return;

        if (debugging("trace"))
	{
            (void)printf("\nEntering rg_string_node_propagate()\n");
	}

        oldp = oldn->posn;
        newp = newn->posn;
        sl = (STATE*)left_state(oldp);
        sr = (STATE*)right_state(oldp);
        newsl = (STATE*)left_state(newp);
        newsr = (STATE*)right_state(newp);

        for (i = 0; i < dim; ++i)
            f[i] = 0.0;
	/* calculate the force from the string chords */
	node_out_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
            b = (*c)->first;
            for (i = 0; i < dim; ++i)
            {
                vec[i] = Coords(b->end)[i] - Coords(b->start)[i];
                vec[i] /= bond_length(b);
                f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
            }
        }
        node_in_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
            b = (*c)->last;
            for (i = 0; i < dim; ++i)
            {
                vec[i] = Coords(b->start)[i] - Coords(b->end)[i];
                vec[i] /= bond_length(b);
                f[i] += kl*(bond_length(b) - bond_length0(b))*vec[i];
            }
        }
	/* propagate the nodes along with the rigid body */
	node_out_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY)
	    {
		b = (*c)->first;
		hs = Hyper_surf(b->_btris[0]->surface);
		hse = Hyper_surf_element(b->_btris[0]->tri);
		break;
	    }
        }
        node_in_curve_loop(oldn,c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY)
	    {
		b = (*c)->last;
		hs = Hyper_surf(b->_btris[0]->surface);
		hse = Hyper_surf_element(b->_btris[0]->tri);
		break;
	    }
        }
	if (hs == NULL || hse == NULL)
	{
	    printf("ERROR in rg_string_node_propagate \n");
	    printf("No related hs or hse found");
	    clean_up(ERROR);
	}
	
    cFluid_point_propagate(front,wave,oldp,newp,hse,hs,dt,V);
    //ifluid_point_propagate(front,wave,oldp,newp,hse,hs,dt,V);

	if (dt > 0.0)
	{
	    for (i = 0; i < dim; ++i)
		accel[i] = (Coords(newp)[i] - Coords(oldp)[i] - 
				oldp->vel[i] * dt) * 2.0 / dt / dt;
	}
	else
	{
	    for (i = 0; i < dim; ++i)
		accel[i] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	    accel[i] -= g[i];

        if (debugging("rigid_body"))
        {
	    (void)printf("accel = %f %f %f\n", accel[0], accel[1], accel[2]);
	    (void)printf("old coords = %f %f %f\n", Coords(oldp)[0], 
				Coords(oldp)[1], Coords(oldp)[2]);
	    (void)printf("oldsl velo = %f %f %f\n", sl->vel[0], 
				sl->vel[1], sl->vel[2]);
	    (void)printf("oldsr velo = %f %f %f\n", sr->vel[0], 
				sr->vel[1], sr->vel[2]);
	    (void)printf("new coords = %f %f %f\n", Coords(newp)[0], 
				Coords(newp)[1], Coords(newp)[2]);
	    (void)printf("newsl velo = %f %f %f\n", newsl->vel[0], 
				newsl->vel[1], newsl->vel[2]);
	    (void)printf("newsr velo = %f %f %f\n", newsr->vel[0], 
				newsr->vel[1], newsr->vel[2]);
	}

	/* Do not change coords, but record the force */
	for (i = 0; i < dim; ++i)
	{
	    Coords(newp)[i] = Coords(oldp)[i];
	    newp->force[i] = f[i];
	    newsl->fluid_accel[i] = newsr->fluid_accel[i] = 
					accel[i] - f[i]/mass;
	    newsr->other_accel[i] = newsl->other_accel[i] = f[i]/mass;
	    newsl->impulse[i] = newsr->impulse[i] = sl->impulse[i];
	}

        if (debugging("trace"))
        {
            (void)printf("Leaving rg_string_node_propagate()\n\n");
        }
}	/* end rg_string_node_propagate */

extern int airfoil_node_propagate(
	Front *front,
	POINTER wave,
	NODE *oldn,
	NODE *newn,
	RPROBLEM        **rp,
        double          dt,
        double          *dt_frac,
        NODE_FLAG       flag,
        POINTER         user)
{
	if (is_load_node(oldn))
	    load_node_propagate(front,oldn,newn,dt);
	else if (is_rg_string_node(oldn))
	    rg_string_node_propagate(front,oldn,newn,dt);
	else
	    return GOOD_NODE;
	return GOOD_NODE;
}	/* end airfoil_node_propagate */

static void passive_curve_propagation(
        Front *front,
        POINTER wave,
        CURVE *oldc,
        CURVE *newc,
        double dt)
{
        BOND *oldb,*newb;
        POINT *oldp,*newp;
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	STATE *sl,*sr;
	STATE *oldst,*btrist;
	double V[MAXD];
	int i,dim = FT_Dimension();

        if (debugging("trace"))
        {
            (void) printf("Entering passive_curve_propagation()\n");
        }
        for (oldb = oldc->first, newb = newc->first; oldb != NULL;
                oldb = oldb->next, newb = newb->next)
        {
            oldp = oldb->end;
            newp = newb->end;
	    for (btris = Btris(oldb); btris && *btris; ++btris)
	    {
		hse = Hyper_surf_element((*btris)->tri);
		hs = Hyper_surf((*btris)->surface);
		FT_GetStatesAtPoint(oldp,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		if (fluid_comp(negative_component(hs)))
		{
		    oldst = (STATE*)left_state(oldp);
		    btrist = (STATE*)sl;
		}
		else if (fluid_comp(positive_component(hs)))
		{
		    oldst = (STATE*)right_state(oldp);
		    btrist = (STATE*)sr;
		}
		for (i = 0; i < dim; ++i)
		    btrist->vel[i] = oldst->vel[i];
	    }
            cFluid_point_propagate(front,wave,oldp,newp,
                        Hyper_surf_element(oldb->_btris[0]->tri),
                        Hyper_surf(oldb->_btris[0]->surface),dt,V);
            /*
            ifluid_point_propagate(front,wave,oldp,newp,
                        Hyper_surf_element(oldb->_btris[0]->tri),
                        Hyper_surf(oldb->_btris[0]->surface),dt,V);
            */
        }
        if (debugging("trace"))
        {
            (void) printf("Leaving passive_curve_propagation()\n\n");
        }
}       /* end passive_curve_propagation */

extern void coating_mono_hyper_surf(
	Front *front)
{
    coating_mono_hyper_surf3d(front);
}	/* end coating_mono_hyper_surf */

static void coating_mono_hyper_surf3d(
	Front *front)
{
	INTERFACE *grid_intfc = front->grid_intfc;
	RECT_GRID *top_grid = &topological_grid(grid_intfc);
	struct Table *T = table_of_interface(grid_intfc);
	COMPONENT *top_comp = T->components;
        COMPONENT          comp;
        INTERFACE          *intfc = front->interf;
	double 		   *L = top_grid->L;
	double 		   *h = top_grid->h;
	double             coords[MAXD];
        double             t[MAXD],p[MAXD],vec[MAXD];
	const double 	   *nor;
	SURFACE **s,*immersed_surf;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	COMPONENT base_comp;
	int i,index,nb,index_nb,*top_gmax = top_grid->gmax;
	int dim = top_grid->dim;
	int icoords[MAXD],icn[MAXD],smin[MAXD],smax[MAXD];
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};

	if (debugging("trace"))
	    (void) printf("Entering coating_mono_hyper_surf3d()\n");
	immersed_surf = NULL;
	for (s = grid_intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		immersed_surf = *s;
		comp = base_comp = negative_component(*s);
		break;
	    }
	}

	if (immersed_surf == NULL)
	    return;

	for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
	for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
	for (icoords[2] = 1; icoords[2] < top_gmax[2]; ++icoords[2])
	{
	    index = d_index(icoords,top_gmax,dim);
	    for (i = 0; i < dim; ++i)
            coords[i] = L[i] + icoords[i]*h[i];
	    
        if (nearest_interface_point_within_range(coords,comp,grid_intfc,
			NO_BOUNDARIES,NULL,p,t,&hse,&hs,3))
	    {
		
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
	    	
            nor = Tri_normal(Tri_of_hse(hse));
	    	for (i = 0; i < dim; ++i)
		        vec[i] = coords[i] - p[i];

	    	if (scalar_product(vec,nor,dim) > 0.0)
                top_comp[index] = base_comp + 1;
	    	else
                top_comp[index] = base_comp - 1;
	    }
	}

	for (s = grid_intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		if (base_comp == negative_component(*s))
		{
		    negative_component(*s) = base_comp - 1;
		    positive_component(*s) = base_comp + 1;
		}
	    }
	}
	if (debugging("coat_comp"))
	{
	    icoords[0] = top_gmax[0]/2;
	    for (icoords[2] = 0; icoords[2] <= top_gmax[2]; ++icoords[2])
	    {
	    	for (icoords[1] = 0; icoords[1] <= top_gmax[1]; ++icoords[1])
	    	{
		    index = d_index(icoords,top_gmax,dim);
		    printf("%d",top_comp[index]);
	    	}
	    	printf("\n");
	    }
	}
	if (debugging("immersed_surf") && front->step%1 == 0)
	{
	    int icrd_nb[MAXD],index_nb,n;
	    POINTER l_state,u_state;
	    double crx_coords[MAXD],crx_nb[MAXD];
	    static double *pl,*pu,*vz,*x;
	    FILE *pfile;
	    char pname[200];

	    n = 0;
	    if (pu == NULL)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&pu,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&pl,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&vz,top_gmax[1],sizeof(double));
	    	FT_VectorMemoryAlloc((POINTER*)&x,top_gmax[1],sizeof(double));
	    }
	    icrd_nb[0] = icoords[0] = top_gmax[0]/2;
	    for (icoords[1] = 0; icoords[1] <= top_gmax[1]; ++icoords[1])
	    {
	    	icrd_nb[1] = icoords[1];
	        for (icoords[2] = 2; icoords[2] < top_gmax[2]-1; ++icoords[2])
		{
		    index = d_index(icoords,top_gmax,dim);
	    	    icrd_nb[2] = icoords[2] + 1;
		    index_nb = d_index(icrd_nb,top_gmax,dim);
		    if (top_comp[index] != top_comp[index_nb] &&
			FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
                                UPPER,top_comp[index],&l_state,&hs,crx_coords)
                        &&
                        FT_StateStructAtGridCrossing(front,grid_intfc,icrd_nb,
                                LOWER,top_comp[index_nb],&u_state,&hs,
                                crx_coords))
		    {
			pl[n] = getStatePres(l_state);
			pu[n] = getStatePres(u_state);
			vz[n] = getStateZvel(l_state);
			x[n] = crx_coords[1];
			n++;
		    }
		}
	    }
	    sprintf(pname,"cpres-%d.xg",front->step);
	    pfile = fopen(pname,"w");
	    fprintf(pfile,"\"Lower pressure\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pl[i]);
	    fprintf(pfile,"\n\n\"Upper pressure\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pu[i]);
	    fprintf(pfile,"\n\n\"Pressure difference\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],pl[i]-pu[i]);
	    fclose(pfile);
	    sprintf(pname,"cvelz-%d.xg",front->step);
	    pfile = fopen(pname,"w");
	    fprintf(pfile,"\"Z-velocity\"\n");
	    for (i = 0; i < n; ++i)
		fprintf(pfile,"%f %f\n",x[i],vz[i]);
	    fclose(pfile);
	}
	if (debugging("trace"))
	    (void) printf("Leaving coating_mono_hyper_surf3d()\n");
}	/* end coating_mono_hyper_surf3d */

