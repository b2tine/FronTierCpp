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

#if defined(__GPU__)
#include "airfoil_gpu.cuh"
#endif

#define		MAX_SURF_CURVES		10
#define		MAX_SURF_NODES		20

static void set_equilibrium_mesh2d(Front*);
static void set_equilibrium_mesh3d(Front*);
static void adjust_for_node_type(NODE*,int,STRING_NODE_TYPE,double**,double**,
				double,double,double*);
static void adjust_for_curve_type(CURVE*,int,double**,double**,double,double*);
static void adjust_for_cnode_type(NODE*,int,double**,double**,double,double*);
static void propagate_curve(ELASTIC_SET*,CURVE*,double**);

static void set_special_node_type(NODE*,int,STRING_NODE_TYPE,SPRING_VERTEX*,
				double,double,double*);

static boolean find_crx_between_box_line(double*,double*,double*,double*,double,int,double*);


extern void second_order_elastic_curve_propagate(
	Front           *fr,
        Front           *newfr,
        INTERFACE       *intfc,
        CURVE           *oldc,
        CURVE           *newc,
        double           fr_dt)
{
	static int size = 0;
	static double **x_old,**x_new,**v_old,**v_new,**f_old,**f_new;
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	int i,j,num_pts,count;
	int n,n_sub = af_params->n_sub;
	double dt_tol,dt = fr_dt/(double)n_sub;
	NODE *ns,*ne;
	int is,ie;
        double *g = af_params->gravity;
        double mass, payload = af_params->payload;
        ELASTIC_SET geom_set;
        STRING_NODE_TYPE start_type = af_params->start_type;
        STRING_NODE_TYPE end_type = af_params->end_type;
	void (*compute_node_accel)(ELASTIC_SET*,NODE*,double**,
				double**,double **,int*);
	void (*compute_curve_accel)(ELASTIC_SET*,CURVE*,double**,
				double**,double **,int*);

	switch (af_params->spring_model)
	{
	case MODEL1:
	    compute_curve_accel = compute_curve_accel1;
	    compute_node_accel = compute_node_accel1;
	    break;
	case MODEL2:
	    compute_curve_accel = compute_curve_accel2;
	    compute_node_accel = compute_node_accel2;
	    break;
	case MODEL3:
	    compute_curve_accel = compute_curve_accel3;
	    compute_node_accel = compute_node_accel3;
	    break;
	default:
	    (void) printf("Model function not implemented yet!\n");
	    clean_up(ERROR);
	}


	if (wave_type(newc) != ELASTIC_BOUNDARY &&
	    wave_type(newc) != ELASTIC_STRING)
	    return;
	if (debugging("trace"))
	    (void) printf("Entering "
			"second_order_elastic_curve_propagate()\n");
	dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
	if (dt > dt_tol)
        {
            n_sub = (int)(fr_dt/dt_tol);
            dt = fr_dt/(double)n_sub;
        }

	num_pts = I_NumOfCurvePoints(oldc);
	if (size < num_pts)
	{
	    size = num_pts;
	    if (x_old != NULL)
		free_these(6,x_old,x_new,v_old,v_new,f_old,f_new);
            FT_MatrixMemoryAlloc((POINTER*)&x_old,size,2,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,2,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,2,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,2,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&f_old,size,2,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&f_new,size,2,sizeof(double));
	}

	geom_set.front = fr;
        geom_set.kl = af_params->kl;
        geom_set.lambda_l = af_params->lambda_l;
        geom_set.m_l = mass = af_params->m_l;
        geom_set.dt = dt;

	ns = newc->start;	ne = newc->end;
	is = 0;                 ie = size - 1;

	count = 0;
        compute_node_accel(&geom_set,ns,f_old,x_old,v_old,&count);
        compute_curve_accel(&geom_set,newc,f_old,x_old,v_old,&count);
        compute_node_accel(&geom_set,ne,f_old,x_old,v_old,&count);

	for (n = 0; n < n_sub; ++n)
	{
	    adjust_for_node_type(ns,is,start_type,f_old,v_old,mass,payload,g);
            adjust_for_node_type(ne,ie,end_type,f_old,v_old,mass,payload,g);
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < 2; ++j)
	    {
	    	x_new[i][j] = x_old[i][j] + v_old[i][j]*dt;
	    	v_new[i][j] = v_old[i][j] + f_old[i][j]*dt;
	    }

	    count = 0;
	    assign_node_field(ns,x_new,v_new,&count);
	    assign_curve_field(newc,x_new,v_new,&count);
	    assign_node_field(ne,x_new,v_new,&count);
	    count = 0;
            compute_node_accel(&geom_set,ns,f_new,x_new,v_new,&count);
            compute_curve_accel(&geom_set,newc,f_new,x_new,v_new,&count);
            compute_node_accel(&geom_set,ne,f_new,x_new,v_new,&count);
            adjust_for_node_type(ns,is,start_type,f_new,v_new,mass,payload,g);
            adjust_for_node_type(ne,ie,end_type,f_new,v_new,mass,payload,g);

	    for (i = 0; i < size; ++i)
	    for (j = 0; j < 2; ++j)
	    {
	    	x_new[i][j] = x_old[i][j] + 0.5*dt*(v_old[i][j] + v_new[i][j]);
	    	v_new[i][j] = v_old[i][j] + 0.5*dt*(f_old[i][j] + f_new[i][j]);
	    }
	    count = 0;
	    propagate_curve(&geom_set,newc,x_new);
	    assign_node_field(ns,x_new,v_new,&count);
	    assign_curve_field(newc,x_new,v_new,&count);
	    assign_node_field(ne,x_new,v_new,&count);
	    if (n != n_sub-1)
	    {
		count = 0;
                compute_node_accel(&geom_set,ns,f_old,x_old,v_old,&count);
                compute_curve_accel(&geom_set,newc,f_old,x_old,v_old,&count);
                compute_node_accel(&geom_set,ne,f_old,x_old,v_old,&count);
	    }
	}
	
	if (debugging("trace"))
	    (void) printf("Leaving "
			"second_order_elastic_curve_propagate()\n");
}	/* end second_order_elastic_curve_propagate */

extern void second_order_elastic_surf_propagate(
        Front           *newfr,
        double           fr_dt)
{
	static int size = 0;
	static double **x_old,**x_new,**v_old,**v_new,**f_old,**f_new;
	AF_PARAMS *af_params = (AF_PARAMS*)newfr->extra2;
        double *g = af_params->gravity;
	double mass;
	int i,j,num_pts,count;
	int n,n0,n_sub = af_params->n_sub;
	double dt_tol,dt = fr_dt/(double)n_sub;
	ELASTIC_SET geom_set;
	CURVE **nc,*newc[MAX_SURF_CURVES];
	SURFACE *news,**s;
	NODE *newn[MAX_SURF_NODES];
	int num_nodes,num_curves;	/* Numbers of nodes and curves */
	boolean in_list;
	void (*compute_node_accel)(ELASTIC_SET*,NODE*,double**,
				double**,double **,int*);
	void (*compute_curve_accel)(ELASTIC_SET*,CURVE*,double**,
				double**,double **,int*);
	void (*compute_surf_accel)(ELASTIC_SET*,SURFACE*,double**,
				double**,double **,int*);

	for (s = newfr->interf->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		news = *s;
		break;
	    }
	}

	switch (af_params->spring_model)
	{
	case MODEL1:
	    compute_curve_accel = compute_curve_accel1;
	    compute_node_accel = compute_node_accel1;
	    compute_surf_accel = compute_surf_accel1;
	    break;
	case MODEL2:
	    compute_curve_accel = compute_curve_accel2;
	    compute_node_accel = compute_node_accel2;
	    compute_surf_accel = compute_surf_accel2;
	    break;
	case MODEL3:
	default:
	    (void) printf("Model function not implemented yet!\n");
	    clean_up(ERROR);
	}

	if (debugging("trace"))
	    (void) printf("Entering "
			"second_order_elastic_surf_propagate()\n");

	dt_tol = sqrt((af_params->m_s)/(af_params->ks))/10.0;
        if (af_params->m_l != 0.0 &&
	    dt_tol > sqrt((af_params->m_l)/(af_params->kl))/10.0)
            dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
        if (af_params->m_g != 0.0 &&
	    dt_tol > sqrt((af_params->m_g)/(af_params->kg))/10.0)
            dt_tol = sqrt((af_params->m_g)/(af_params->kg))/10.0;
	if (dt > dt_tol)
        {
            n_sub = (int)(fr_dt/dt_tol);
            dt = fr_dt/(double)n_sub;
        }
	geom_set.ks = af_params->ks;
	geom_set.lambda_s = af_params->lambda_s;
	geom_set.m_s = mass = af_params->m_s;
	geom_set.kg = af_params->kg;
	geom_set.lambda_g = af_params->lambda_g;
	geom_set.m_g = af_params->m_g;
	geom_set.kl = af_params->kl;
	geom_set.lambda_l = af_params->lambda_l;
	geom_set.m_l = af_params->m_l;
	geom_set.front = newfr;
	geom_set.dt = dt;
	if (debugging("step_size"))
	{
	    (void) printf("n_sub = %d\n",n_sub);
	    (void) printf("ks = %f  kl = %f\n",geom_set.ks,geom_set.kl);
	    (void) printf("m_s = %f  m_l = %f\n",geom_set.m_s,geom_set.m_l);
	    (void) printf("lambda_s = %f  lambda_l = %f\n",
				geom_set.lambda_s,geom_set.lambda_l);
	    (void) printf("fr_dt = %f  dt_tol = %20.14f  dt = %20.14f\n",
                                fr_dt,dt_tol,dt);
            (void) printf("Number of interior sub-steps = %d\n",n_sub);
	}

	/* Assume there is only one closed boundary curve */
	num_nodes = num_curves = 0;
	for (nc = news->pos_curves; nc && *nc; ++nc)
	{
	    if (hsbdry_type(*nc) == FIXED_HSBDRY) continue;
	    if (!pointer_in_list((POINTER)(*nc),num_curves,
				(POINTER*)newc))
	    {
	    	newc[num_curves] = *nc;
	    	num_curves++;
	    }
	}
	for (nc = news->neg_curves; nc && *nc; ++nc)
	{
	    if (hsbdry_type(*nc) == FIXED_HSBDRY) continue;
	    if (!pointer_in_list((POINTER)(*nc),num_curves,
				(POINTER*)newc))
	    {
	    	newc[num_curves] = *nc;
	    	num_curves++;
	    }
	}
	for (i = 0; i < num_curves; ++i)
	{
	    if (!pointer_in_list((POINTER)newc[i]->start,num_nodes,
				(POINTER*)newn))
	    {
		newn[num_nodes] = newc[i]->start;
		num_nodes++;
	    }
	    if (is_closed_curve(newc[i])) continue;
	    if (!pointer_in_list((POINTER)newc[i]->end,num_nodes,
				(POINTER*)newn))
	    {
		newn[num_nodes] = newc[i]->end;
		num_nodes++;
	    }
	}

	num_pts = I_NumOfSurfPoints(news);
	if (size < num_pts)
	{
	    size = num_pts;
	    if (v_old != NULL)
	    {
		FT_FreeThese(6,v_old,v_new,x_old,x_new,f_old,f_new);
	    }
	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&f_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&f_new,size,3,sizeof(double));
	}

	count = 0;
	compute_surf_accel(&geom_set,news,f_old,x_old,v_old,&count);
	for (i = 0; i < num_curves; ++i)
	{
	    n0 = count;
	    compute_curve_accel(&geom_set,newc[i],f_old,x_old,v_old,&count);
	    adjust_for_curve_type(newc[i],n0,f_old,v_old,mass,g);
	}
	for (i = 0; i < num_nodes; ++i)
	{
	    n0 = count;
	    compute_node_accel(&geom_set,newn[i],f_old,x_old,v_old,&count);
	    adjust_for_cnode_type(newn[i],n0,f_old,v_old,mass,g);
	}

	for (n = 0; n < n_sub; ++n)
	{
	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
		x_new[i][j] = x_old[i][j] + dt*v_old[i][j];
                v_new[i][j] = v_old[i][j] + dt*f_old[i][j];
            }
	    count = 0;
            assign_surf_field(news,x_new,v_new,&count);
	    for (i = 0; i < num_curves; ++i)
            	assign_curve_field(newc[i],x_new,v_new,&count);
	    for (i = 0; i < num_nodes; ++i)
            	assign_node_field(newn[i],x_new,v_new,&count);
	    count = 0;
	    compute_surf_accel(&geom_set,news,f_new,x_new,v_new,&count);
	    for (i = 0; i < num_curves; ++i)
	    {
		n0 = count;
	    	compute_curve_accel(&geom_set,newc[i],f_new,x_new,v_new,&count);
	    	adjust_for_curve_type(newc[i],n0,f_new,v_new,mass,g);
	    }
	    for (i = 0; i < num_nodes; ++i)
	    {
		n0 = count;
	    	compute_node_accel(&geom_set,newn[i],f_new,x_new,v_new,&count);
	    	adjust_for_cnode_type(newn[i],n0,f_old,v_old,mass,g);
	    }

	    for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
		x_new[i][j] = x_old[i][j] + 0.5*dt*(v_old[i][j] + v_new[i][j]);
                v_new[i][j] = v_old[i][j] + 0.5*dt*(f_old[i][j] + f_new[i][j]);
            }
	    count = 0;
            assign_surf_field(news,x_new,v_new,&count);
	    for (i = 0; i < num_curves; ++i)
	    {
            	assign_curve_field(newc[i],x_new,v_new,&count);
	    }
	    for (i = 0; i < num_nodes; ++i)
            	assign_node_field(newn[i],x_new,v_new,&count);
	    if (n != n_sub-1)
	    {
	    	count = 0;
	    	compute_surf_accel(&geom_set,news,f_old,x_old,v_old,&count);
	    	for (i = 0; i < num_curves; ++i)
		{
		    n0 = count;
	    	    compute_curve_accel(&geom_set,newc[i],f_old,x_old,
						v_old,&count);
	    	    adjust_for_curve_type(newc[i],n0,f_old,v_old,mass,g);
		}
	    	for (i = 0; i < num_nodes; ++i)
		{
		    n0 = count;
	    	    compute_node_accel(&geom_set,newn[i],f_old,x_old,
						v_old,&count);
	    	    adjust_for_cnode_type(newn[i],n0,f_old,v_old,mass,g);
	    	}
	    }
	}

	if (debugging("trace"))
	    (void) printf("Leaving "
			"second_order_elastic_surf_propagate()\n");
}	/* end second_order_elastic_surf_propagate */

extern void set_equilibrium_mesh(
	Front *front)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    set_equilibrium_mesh2d(front);
	    return;
	case 3:
	    set_equilibrium_mesh3d(front);
	    return;
	}
}	/* end set_equilibrium_mesh */

static void set_equilibrium_mesh2d(
	Front *front)
{
	CURVE **c,*curve;
	BOND *b;
	short unsigned int seed[3] = {2,72,7172};
	double len0,total_length = 0.0;
	int i,n = 0;
	INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int dim = front->rect_grid->dim;

	for (c = intfc->curves; c && *c; ++c)
	{
	    if (wave_type(*c) != ELASTIC_BOUNDARY &&
		wave_type(*c) != ELASTIC_STRING) 
		continue;
	    curve = *c;
	    for (b = curve->first; b != NULL; b = b->next)
	    {
		total_length +=  bond_length(b);
		n++;
	    }
	    len0 = total_length/(double)n;
	    for (b = curve->first; b != NULL; b = b->next)
	    {
		b->length0 = len0;
		for (i = 0; i < dim; ++i)
		    b->dir0[i] = (Coords(b->end)[i] - Coords(b->start)[i])
				/b->length0;	
	    }
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		if (af_params->pert_params.pert_type == PARALLEL_RAND_PERT) 
		{
		    double dx_rand;
		    double tan[MAXD];
		    tangent(b->end,b,curve,tan,front); 
		    dx_rand = (erand48(seed) - 0.5)*len0;
		    for (i = 0; i < dim; ++i)
		    	Coords(b->end)[i] += 0.2*dx_rand*tan[i];
		}
		else if (af_params->pert_params.pert_type == 
				ORTHOGONAL_RAND_PERT)
		{
		    double dx_rand;
		    double nor[MAXD];
		    double amp = af_params->pert_params.pert_amp;
		    FT_NormalAtPoint(b->end,front,nor,NO_COMP); 
		    if (amp > 1.0) amp = 1.0;
		    dx_rand = (erand48(seed) - 0.5)*amp*len0;
		    for (i = 0; i < dim; ++i)
		    	Coords(b->end)[i] += dx_rand*nor[i];
		}
		else if (af_params->pert_params.pert_type == SINE_PERT)
		{
		    /* This assumes the curve is horizontal */
		    double amp = af_params->pert_params.pert_amp;
		    double L = Coords(curve->start->posn)[0]; 
		    double U = Coords(curve->end->posn)[0];
		    double x = Coords(b->end)[0];
		    for (i = 1; i < dim; ++i)
		    	Coords(b->end)[i] += amp*sin(PI*(x-L)/(U-L));
		}
	    }
	    for (b = curve->first; b != NULL; b = b->next)
		set_bond_length(b,2);
	    never_redistribute(Hyper_surf(curve)) = YES;
	}
}	/* end set_equilibrium_mesh2d */

//TODO: Compare to set{Curve,Surf}ZeroMesh() functions in cgal.cpp
//      Consolidate functionality if possible.
static void set_equilibrium_mesh3d(
	Front *front)
{
	SURFACE **s,*surf;
	CURVE **c;
	TRI *t;
	int i,j,dir;
	short unsigned int seed[3] = {2,72,7172};
	double max_len,min_len,ave_len,len;
	double vec[3],*cen,radius,amp;
	double x0,xl,xu;
	double dx_rand;
	double count;
	INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	BOND *b;
	double gore_len_fac = af_params->gore_len_fac;

	for (c = intfc->curves; c && *c; ++c)
	{
        if (hsbdry_type(*c) != STRING_HSBDRY ||
            hsbdry_type(*c) != GORE_HSBDRY) continue;

	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
            set_bond_length(b,3);
            b->length0 = bond_length(b);
            
            if (hsbdry_type(*c) == GORE_HSBDRY)
                b->length0 *= gore_len_fac;
            
            for (i = 0; i < 3; ++i)
            {
                b->dir0[i] =
                    (Coords(b->end)[i] - Coords(b->start)[i])/b->length0;	
            }
	    }
    
        never_redistribute(Hyper_surf(*c)) = YES;
	}

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
	    surf = *s;
	    ave_len = 0.0;
	    max_len = 0.0;
	    min_len = HUGE;
	    count = 0.0;
	    for (t = first_tri(surf); !at_end_of_tri_list(t,surf); t = t->next)
	    {
		for (i = 0; i < 3; ++i)
		{
		    t->side_length0[i] = separation(Point_of_tri(t)[i],
			Point_of_tri(t)[(i+1)%3],3);
		    for (j = 0; j < 3; ++j)
		    {
		    	t->side_dir0[i][j] = 
				(Coords(Point_of_tri(t)[(i+1)%3])[j] -
				 Coords(Point_of_tri(t)[i])[j])/
				 t->side_length0[i];
		    }
		    if (max_len < t->side_length0[i]) 
			max_len = t->side_length0[i];
		    if (min_len > t->side_length0[i])
			min_len = t->side_length0[i];
		    ave_len += t->side_length0[i];
		    count += 1.0;
		}
	    }
	    never_redistribute(Hyper_surf(surf)) = YES;
	}
	printf("Original length:\n");
	printf("min_len = %16.12f\n",min_len);
	printf("max_len = %16.12f\n",max_len);
	printf("ave_len = %16.12f\n",ave_len/count);

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
	    surf = *s;
	    switch (af_params->pert_params.pert_type)
	    {
	    case ORTHOGONAL_RAND_PERT:
	        for (t = first_tri(surf); !at_end_of_tri_list(t,surf); 
				t = t->next)
		for (i = 0; i < 3; ++i)
		{
		    for (j = 0; j < 3; ++j)
		    {
			vec[j] = Coords(Point_of_tri(t)[i])[j] -
				 Coords(Point_of_tri(t)[(i+1)%3])[j];
		    }
		    for (j = 0; j < 3; ++j)
		    {
		    	dx_rand = (2.0 + erand48(seed))/3.0;
			vec[j] *= dx_rand;
			Coords(Point_of_tri(t)[j])[2] += vec[j]; 
		    }
		}
		break;
	    case PARALLEL_RAND_PERT:
	        for (t = first_tri(surf); !at_end_of_tri_list(t,surf); 
				t = t->next)
		for (i = 0; i < 3; ++i)
		{
		    for (j = 0; j < 3; ++j)
		    {
			vec[j] = Coords(Point_of_tri(t)[i])[j] -
				 Coords(Point_of_tri(t)[(i+1)%3])[j];
		    }
		    dx_rand = (2.0 + erand48(seed))/3.0;
		    for (j = 0; j < 3; ++j)
		    {
			vec[j] *= dx_rand;
			Coords(Point_of_tri(t)[i])[j] = vec[j] +
				Coords(Point_of_tri(t)[(i+1)%3])[j];
		    }
		}
		break;
	    case RADIAL_PERT:
		cen = af_params->pert_params.cen;
		radius = af_params->pert_params.pert_radius;
		amp = af_params->pert_params.pert_amp;
	        for (t = first_tri(surf); !at_end_of_tri_list(t,surf); 
				t = t->next)
		for (i = 0; i < 3; ++i)
		    sorted(Point_of_tri(t)[i]) = NO;
	        for (t = first_tri(surf); !at_end_of_tri_list(t,surf); 
				t = t->next)
		for (i = 0; i < 3; ++i)
		{
		    POINT *p = Point_of_tri(t)[i];
		    double r;
		    if (Boundary_point(p)) continue;
		    if (sorted(p)) continue;
		    r = sqr(Coords(p)[0] - cen[0]) + 
			sqr(Coords(p)[1] - cen[1]);
		    r = sqrt(r);
		    if (r < radius)
		    	Coords(p)[2] += amp*(1.0 - r/radius);
		    sorted(p) = YES;
		}
		break;
	    case LINEAR_PERT:
		x0 = af_params->pert_params.x0;
		xl = af_params->pert_params.xl;
		xu = af_params->pert_params.xu;
		amp = af_params->pert_params.pert_amp;
		dir = af_params->pert_params.dir;
	        for (t = first_tri(surf); !at_end_of_tri_list(t,surf); 
				t = t->next)
		for (i = 0; i < 3; ++i)
		    sorted(Point_of_tri(t)[i]) = NO;
	        for (t = first_tri(surf); !at_end_of_tri_list(t,surf); 
				t = t->next)
		for (i = 0; i < 3; ++i)
		{
		    POINT *p = Point_of_tri(t)[i];
		    if (sorted(p)) continue;
		    if (Coords(p)[dir] < x0)
			Coords(p)[2] += amp*(Coords(p)[dir] - xl)/(x0 - xl);
		    else
			Coords(p)[2] += amp*(xu - Coords(p)[dir])/(xu - x0);
		    sorted(p) = YES;
		}
		break;
	    case NO_PERT:
	    default:
		break;
	    }
	}
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
	    surf = *s;
	    ave_len = 0.0;
	    max_len = 0.0;
	    min_len = HUGE;
	    count = 0.0;
	    for (t = first_tri(surf); !at_end_of_tri_list(t,surf); t = t->next)
	    {
		for (i = 0; i < 3; ++i)
		{
		    len = separation(Point_of_tri(t)[i],
			Point_of_tri(t)[(i+1)%3],3);
		    if (max_len < len) 
			max_len = len;
		    if (min_len > len)
			min_len = len;
		    ave_len += len;
		    count += 1.0;
		}
	    }
	}
	printf("Perturbed length:\n");
	printf("min_len = %16.12f\n",min_len);
	printf("max_len = %16.12f\n",max_len);
	printf("ave_len = %16.12f\n",ave_len/count);
}	/* end set_equilibrium_mesh3d */

#define 	MAX_NUM_RING1		30

extern void fourth_order_elastic_curve_propagate(
	Front           *fr,
        Front           *newfr,
        INTERFACE       *intfc,
        CURVE           *oldc,
        CURVE           *newc,
        double           fr_dt)
{
	static int size = 0;
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	double mass;
	int i,j,num_pts,count;
	int n,n_sub = af_params->n_sub;
	double dt_tol,dt = fr_dt/(double)n_sub;
	int dim = fr->rect_grid->dim;
	double *g = af_params->gravity;
	double payload = af_params->payload;
	ELASTIC_SET geom_set;
	STRING_NODE_TYPE start_type = af_params->start_type;
	STRING_NODE_TYPE end_type = af_params->end_type;
	void (*compute_node_accel)(ELASTIC_SET*,NODE*,double**,
				double**,double **,int*);
	void (*compute_curve_accel)(ELASTIC_SET*,CURVE*,double**,
				double**,double **,int*);
	static SPRING_VERTEX *sv;
	static boolean first = YES;
	static GLOBAL_POINT **point_set;
        static GLOBAL_POINT *point_set_store;
        static GLOBAL_POINT *local_point_store;
	long max_point_gindex = fr->interf->max_point_gindex;

	if (wave_type(newc) != ELASTIC_BOUNDARY &&
	    wave_type(newc) != ELASTIC_STRING)
	    return;
	if (debugging("trace"))
	    (void) printf("Entering "
			"fourth_order_elastic_curve_propagate()\n");

	dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
	dt_tol = std::min(dt_tol, sqrt((af_params->m_s)/(af_params->ks))/10.0);
	if (dt > dt_tol)
        {
            n_sub = (int)(fr_dt/dt_tol);
            dt = fr_dt/(double)n_sub;
        }

	if (point_set == NULL)
        {
            FT_VectorMemoryAlloc((POINTER*)&point_set,max_point_gindex,
                                        sizeof(GLOBAL_POINT*));
            for (i = 0; i < max_point_gindex; ++i)
                point_set[i] = NULL;
        }
	num_pts = I_NumOfCurvePoints(oldc);
	geom_set.num_surfs = 0;
	geom_set.num_curves = 1;
	geom_set.curves[0] = newc;
	geom_set.num_nodes = 2;
	geom_set.nodes[0] = newc->start;
	geom_set.nodes[1] = newc->end;
	if (size < num_pts)
	{
	    size = num_pts;
	    if (sv != NULL)
	    {
		FT_FreeThese(1,sv);
	    }
            FT_VectorMemoryAlloc((POINTER*)&sv,size,sizeof(SPRING_VERTEX));
	    FT_VectorMemoryAlloc((POINTER*)&point_set_store,size,
                                        sizeof(GLOBAL_POINT));
            link_point_set(&geom_set,point_set,point_set_store);
	    FT_VectorMemoryAlloc((POINTER*)&local_point_store,size,
                                        sizeof(GLOBAL_POINT));
	}

	geom_set.front = fr;
	geom_set.ks = af_params->ks;
        geom_set.lambda_s = af_params->lambda_s;
        geom_set.m_s = af_params->m_s;
	geom_set.kl = af_params->kl;
	geom_set.lambda_l = af_params->lambda_l;
	geom_set.m_l = mass = af_params->m_l;
	geom_set.dt = dt;

	count_vertex_neighbors(&geom_set,sv);
	if (first)
	{
	    set_spring_vertex_memory(sv,size);
	    first = NO;
	}
	set_vertex_neighbors(&geom_set,sv,point_set);
        get_point_set_from(&geom_set,point_set);
	set_special_node_type(newc->start,size-2,start_type,sv,mass,payload,g);
	set_special_node_type(newc->end,size-1,end_type,sv,mass,payload,g);

	/* Start intensive computation */

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
	put_point_set_to(&geom_set,point_set);
	set_vertex_impulse(&geom_set,point_set);

	/* End intensive computation */
	
	if (debugging("trace"))
	    (void) printf("Leaving "
			"fourth_order_elastic_curve_propagate()\n");
}	/* end fourth_order_elastic_curve_propagate */

extern void fixed_length_tan_curve_propagate(
	Front           *fr,
        Front           *newfr,
        INTERFACE       *intfc,
        CURVE           *oldc,
        CURVE           *newc,
        double           dt)
{
	BOND *b,*bs,*bs2,*be;
	int nb,n;
	double seg_len,total_length,total_length0;
	static boolean first = YES;

	if (debugging("trace"))
	    (void) printf("Entering fixed_length_tan_curve_propagate()\n");
	if (fr->rect_grid->dim != 2) return;

	if (wave_type(newc) != ELASTIC_BOUNDARY &&
	    wave_type(newc) != ELASTIC_STRING)
	    return;

	if (first)
	{
	    first = NO;
	    nb = 0;
	    total_length  = 0.0;
	    for (b = newc->first; b != NULL; b = b->next)
	    {
	    	total_length += bond_length(b);
		nb++;
	    }
	    for (b = newc->first; b != NULL; b = b->next)
		b->length0 = total_length/(double)nb;
	}
	if (debugging("airfoil"))
	{
	    printf("Entering fixed_length_tan_curve_propagate()\n");
	    nb = 0;
	    total_length  = 0.0;
	    total_length0 = 0.0;
	    for (b = newc->first; b != NULL; b = b->next)
	    {
	    	total_length += bond_length(b);
	    	total_length0 += b->length0;
		nb++;
	    }
	    printf("Entering: total_length  = %16.12f\n", total_length);
	    printf("Entering: total_length0 = %16.12f\n", total_length0);
	}

	nb = 0;
	for (b = newc->first; b != NULL; b = b->next) nb++;
	if (nb%2 == 0)	n = nb/2;
	else n = nb/2 + 1;

	seg_len = 0.0;
	bs = newc->first;
	for (nb = 0, b = bs; nb < n; b = b->next) 
	{
	    nb++;
	    seg_len += b->length0;
	    be = b;
	}
	bs2 = be->next;
	FT_CurveSegLengthConstr(newc,bs,be,nb,seg_len,
				BACKWARD_REDISTRIBUTION);

	seg_len = 0.0;
	bs = bs2;
	for (nb = 0, b = bs; b != NULL; b = b->next) 
	{
	    nb++;
	    seg_len += b->length0;
	    be = b;
	}
	FT_CurveSegLengthConstr(newc,bs,be,nb,seg_len,
				FORWARD_REDISTRIBUTION);
	
	total_length = 0.0;
	for (b = newc->first; b != NULL; b = b->next)
	{
	    total_length += bond_length(b);
	    b->length0 = seg_len/(double)nb;
	}
	if (debugging("trace"))
	    (void) printf("Leaving fixed_length_tan_curve_propagate()\n");
}	/* end fixed_length_tan_curve_propagate */

extern void fourth_order_elastic_surf_propagate(
        Front           *newfr,
        double           fr_dt)
{
	static int size = 0;
	AF_PARAMS *af_params = (AF_PARAMS*)newfr->extra2;
        double *g = af_params->gravity;
	double mass;
	int i,j,num_pts,count;
	int n,n_sub = af_params->n_sub;
	double dt_tol,dt = fr_dt/(double)n_sub;
	ELASTIC_SET geom_set;
	CURVE **nc,*newc[MAX_SURF_CURVES];
	SURFACE *news,**s;
	NODE *newn[MAX_SURF_NODES];
	int num_nodes,num_curves;	/* Numbers of nodes and curves */
	boolean in_list;
	static SPRING_VERTEX *sv;
        static boolean first = YES;
	int countc[100],countn[100];
	int dim = newfr->rect_grid->dim;
	static GLOBAL_POINT **point_set;
        static GLOBAL_POINT *point_set_store;
        static GLOBAL_POINT *local_point_store;
	long max_point_gindex = newfr->interf->max_point_gindex;
	int owner[MAXD];

	if (pp_numnodes() > 1)
        {
            INTERFACE *elastic_intfc;
            owner[0] = 0;
            owner[1] = 0;
            owner[2] = 0;
            add_to_debug("collect_intfc");
            elastic_intfc = FT_CollectHypersurfFromSubdomains(newfr,owner,
                                ELASTIC_BOUNDARY);
            printf("Parallel code needed!\n");
            clean_up(0);
        }	
	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_surf_propagate()\n");

	start_clock("set_spring_model");
	for (s = newfr->interf->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
	    {
		news = *s;
		break;
	    }
	}

	dt_tol = sqrt((af_params->m_s)/(af_params->ks))/10.0;
        if (af_params->m_l != 0.0 &&
	    dt_tol > sqrt((af_params->m_l)/(af_params->kl))/10.0)
            dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
        if (af_params->m_g != 0.0 &&
	    dt_tol > sqrt((af_params->m_g)/(af_params->kg))/10.0)
            dt_tol = sqrt((af_params->m_g)/(af_params->kg))/10.0;
	if (dt > dt_tol)
        {
            n_sub = (int)(fr_dt/dt_tol);
            dt = fr_dt/(double)n_sub;
        }
	geom_set.ks = af_params->ks;
	geom_set.lambda_s = af_params->lambda_s;
	geom_set.m_s = mass = af_params->m_s;
	geom_set.kg = af_params->kg;
	geom_set.lambda_g = af_params->lambda_g;
	geom_set.m_g = af_params->m_g;
	geom_set.kl = af_params->kl;
	geom_set.lambda_l = af_params->lambda_l;
	geom_set.m_l = af_params->m_l;
	geom_set.front = newfr;
	geom_set.dt = dt;
	if (debugging("step_size"))
	{
	    (void) printf("n_sub = %d\n",n_sub);
	    (void) printf("ks = %f  kl = %f\n",geom_set.ks,geom_set.kl);
	    (void) printf("m_s = %f  m_l = %f\n",geom_set.m_s,geom_set.m_l);
	    (void) printf("lambda_s = %f  lambda_l = %f\n",
				geom_set.lambda_s,geom_set.lambda_l);
	}
	(void) printf("\nfr_dt = %f  dt_tol = %20.14f  dt = %20.14f\n",
                                fr_dt,dt_tol,dt);
        (void) printf("Number of interior sub-steps = %d\n\n",n_sub);

	/* Assume there is only one closed boundary curve */
	num_nodes = num_curves = 0;
	for (nc = news->pos_curves; nc && *nc; ++nc)
	{
	    //if (hsbdry_type(*nc) == FIXED_HSBDRY) continue;
	    if (!pointer_in_list((POINTER)(*nc),num_curves,
				(POINTER*)newc))
	    {
	    	newc[num_curves] = *nc;
	    	geom_set.curves[num_curves] = *nc;
	    	num_curves++;
	    }
	}
	for (nc = news->neg_curves; nc && *nc; ++nc)
	{
	    if (!pointer_in_list((POINTER)(*nc),num_curves,
				(POINTER*)newc))
	    {
	    	newc[num_curves] = *nc;
	    	geom_set.curves[num_curves] = *nc;
	    	num_curves++;
	    }
	}
	geom_set.num_curves = num_curves;
	for (i = 0; i < num_curves; ++i)
	{
	    if (!pointer_in_list((POINTER)newc[i]->start,num_nodes,
				(POINTER*)newn))
	    {
		newn[num_nodes] = newc[i]->start;
	    	geom_set.nodes[num_nodes] = newc[i]->start;
		num_nodes++;
	    }
	    if (is_closed_curve(newc[i])) continue;
	    if (!pointer_in_list((POINTER)newc[i]->end,num_nodes,
				(POINTER*)newn))
	    {
		newn[num_nodes] = newc[i]->end;
	    	geom_set.nodes[num_nodes] = newc[i]->end;
		num_nodes++;
	    }
	}
	geom_set.num_nodes = num_nodes;
	geom_set.num_surfs = 1;
	geom_set.surfs[0] = news;

	num_pts = I_NumOfSurfPoints(news);
	if (point_set == NULL)
        {
            FT_VectorMemoryAlloc((POINTER*)&point_set,max_point_gindex,
                                        sizeof(GLOBAL_POINT*));
            for (i = 0; i < max_point_gindex; ++i)
                point_set[i] = NULL;
        }
	if (size < num_pts && first)
	{
	    size = num_pts;
	    if (sv != NULL)
		FT_FreeThese(1,sv);
	    FT_VectorMemoryAlloc((POINTER*)&sv,size,sizeof(SPRING_VERTEX));
	    FT_VectorMemoryAlloc((POINTER*)&point_set_store,size,
                                        sizeof(GLOBAL_POINT));
            link_point_set(&geom_set,point_set,point_set_store);
	    FT_VectorMemoryAlloc((POINTER*)&local_point_store,size,
                                        sizeof(GLOBAL_POINT));
	}

	count_vertex_neighbors(&geom_set,sv);
	if (first)
	{
	    set_spring_vertex_memory(sv,size);
            first = NO;
	}
	set_vertex_neighbors(&geom_set,sv,point_set);
        get_point_set_from(&geom_set,point_set);

	stop_clock("set_spring_model");

        /* Start intensive computation */

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

	put_point_set_to(&geom_set,point_set);
	set_vertex_impulse(&geom_set,point_set);

	stop_clock("spring_model");

	if (debugging("trace"))
	    (void) printf("Leaving "
			"fourth_order_elastic_surf_propagate()\n");
}	/* end fourth_order_elastic_surf_propagate */

static void adjust_for_curve_type(
	CURVE *c,
	int index0,
	double **f,
	double **v,
	double m_l,
	double *g)
{
	C_PARAMS *c_params =  (C_PARAMS*)c->extra;
	int j,dir,dim = 3;
	double load_point_mass;
	BOND *b;
	int n,index = index0;
	double ave_accel;
	double *force;

	if (c_params == NULL) return;
	force = c_params->force;
	dir = c_params->dir;
	load_point_mass = c_params->point_mass;

	if (c_params->load_type == NO_LOAD) return;
	else if (c_params->load_type == FREE_LOAD)
	{
	    for (b = c->first; b != c->last; b = b->next)
	    {
	    	for (j = 0; j < dim; ++j)
		{
                    f[index][j] = f[index][j]*m_l/load_point_mass + g[j]
				+ force[j];
		}
		index++;
	    }
	}
	else if (c_params->load_type == RIGID_LOAD)
	{
	    ave_accel = 0.0;
	    n = 0;
	    for (b = c->first; b != c->last; b = b->next)
	    {
            	ave_accel += f[index][dir];
		n++;
		index++;
	    }
	    ave_accel /= n;
	    c_params->ave_accel = ave_accel;
	    index = index0;
	    for (b = c->first; b != c->last; b = b->next)
	    {
	    	for (j = 0; j < dim; ++j)
	    	{
		    if (j == dir)
		    {
            	    	f[index][dir] = ave_accel*m_l/load_point_mass + 
					g[dir] + force[dir];
		    }
		    else
		    	f[index][j] = v[index][j] = 0.0;
	    	}
		index++;
	    }
	}
}	/* end adjust_for_curve_type */

static void adjust_for_cnode_type(
	NODE *n,
	int index,
	double **f,
	double **v,
	double m_l,
	double *g)
{
	CURVE **c;
	C_PARAMS *c_params =  NULL;
	int j,dir,dim = 3;
	double load_point_mass;
	double *force;
	double ave_accel;

	for (c = n->in_curves; c && *c; ++c)
	{
	    if ((*c)->extra != NULL)
	    	c_params = (C_PARAMS*)(*c)->extra;
	}
	for (c = n->out_curves; c && *c; ++c)
	{
	    if ((*c)->extra != NULL)
	    	c_params = (C_PARAMS*)(*c)->extra;
	}
	if (c_params == NULL) return;

	dir = c_params->dir;
	force = c_params->force;
	load_point_mass = c_params->point_mass;
	if (c_params->load_type == NO_LOAD) return;
	else if (c_params->load_type == FREE_LOAD)
	{
	    for (j = 0; j < dim; ++j)
                f[index][j] = f[index][j]*m_l/load_point_mass + g[j]
				+ force[j];
	}
	else if (c_params->load_type == RIGID_LOAD)
	{
	    ave_accel = c_params->ave_accel;
	    for (j = 0; j < dim; ++j)
	    {
		if (j == dir)
            	    //f[index][dir] = f[index][dir]*m_l/load_point_mass + g[dir]
            	    f[index][dir] = ave_accel*m_l/load_point_mass + g[dir]
				+ force[dir];
		else
		    f[index][j] = v[index][j] = 0.0;
	    }
	}
}	/* end adjust_for_cnode_type */

static void adjust_for_node_type(
	NODE *n,
	int index,
	STRING_NODE_TYPE end_type,
	double **f,
	double **v,
	double mass,
	double payload,
	double *g)
{
	int j;
	int dim = n->interface->dim;

	if (end_type == FIXED_END)
	{
	    for (j = 0; j < dim; ++j)
		f[index][j] = v[index][j] = 0.0;
	}
	else if (end_type == LOADED_END)
	{
	    for (j = 0; j < dim; ++j)
		f[index][j] = f[index][j]*mass/payload + g[j];
	}
}	/* end adjust_for_node_type */

static void propagate_curve(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	double **x)
{
	int j;
	POINT *p;
	BOND *b;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	double *v;
	double dt = geom_set->dt;
	int n = 1;

	if (debugging("trace"))
	    (void) printf("Entering propagate_curve()\n");
	hs = Hyper_surf(curve);
	for (b = curve->first;  b != curve->last; b = b->next)
	{
	    hse = Hyper_surf_element(b);
	    p = b->end;
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    v = sl->vel;
	    for (j = 0; j < 2; ++j)
		x[n][j] += v[j]*dt;
	    ++n;
	}
}	/* end propagate_curve */

static void set_special_node_type(
	NODE *node,
	int n,
	STRING_NODE_TYPE end_type,
	SPRING_VERTEX *sv,
	double mass,
	double payload,
	double *g)
{
	int j;
	int dim = node->interface->dim;

	if (end_type == FIXED_END)
	{
	    for (j = 0; j < sv[n].num_nb; ++j)
	    	sv[n].k[j] = 0.0;
	    for (j = 0; j < dim; ++j)
		sv[n].ext_accel[j] = 0.0;
	    sv[n].lambda = 0.0;
	}
	else if (end_type == LOADED_END)
	{
	    sv[n].m = payload;
	    for (j = 0; j < dim; ++j)
		sv[n].ext_accel[j] = g[j];
	}
}	/* set_special_node_type */

extern void resolve_wall_collision(
        Front* front,
        SPRING_VERTEX* sv,
        int size)
{
        int i,j,dim;
        double *x, *v, crx[MAXD],center[MAXD];
        double *L, *U;
        boolean out_domain;
        L = front->rect_grid->L;
        U = front->rect_grid->U;
        dim = FT_Dimension();

        for (i = 0; i < size; i++)
        {
            x = sv[i].x;
            v = sv[i].v;
            out_domain = NO;
            for (j = 0; j < dim; j++)
            {
                center[j] = x[j];
                if (x[j] < L[j] || x[j] > U[j])
                {
                    center[j] = 0.5*(U[j] + L[j]);
                    out_domain = YES;
                }
            }
            if (out_domain == YES)
            {
                /*point is out of domain*/
                /*find crossing*/
                if(find_crx_between_box_line(L,U,x,center,0.1,dim,crx))
                {
                    for (j = 0; j < dim; j++)
                        x[j] = crx[j];
                }
                else
                {
                    printf("find crx failed, point is in the box\n");
                    clean_up(ERROR);
                }
            }
        }
}

static boolean find_crx_between_box_line(
        double* bmin, /*lower point of a box*/
        double* bmax, /*upper point of a box*/
        double* s,   /*start of line segment*/
        double* e,   /*end of line segment*/
        double  tol, /*shift tol to end*/
        int     dim, /*dimension*/
        double* ans)  /*intersection*/
{
        int i;
        double d, st, et, fst = 0.0, fet = 1.0;
        double times;
        for (i = 0; i < dim; i++)
        {
            if (s[i] < e[i])
            {
                if (s[i] > bmax[i] || e[i] < bmin[i])
                    return NO;
                d = e[i] - s[i];
                st = (s[i] < bmin[i]) ? (bmin[i] - s[i])/d : 0.0;
                et = (e[i] > bmax[i]) ? (bmax[i] - s[i])/d : 1.0;
            }
            else
            {
                if (e[i] > bmax[i] || s[i] < bmin[i])
                    return NO;
                d = e[i] - s[i];
                st = (s[i] > bmax[i]) ? (bmax[i] - s[i])/d : 0.0;
                et = (e[i] < bmin[i]) ? (bmin[i] - s[i])/d : 1.0;
            }
            if (st > fst) fst = st;
            if (et < fet) fet = et;
            if (fet < fst)
                return NO;
        }
        times = fst;
        for (i = 0; i < dim; i++)
            ans[i] = s[i] + (e[i] - s[i]) * times *(1.0 + tol);
        return YES;
}

