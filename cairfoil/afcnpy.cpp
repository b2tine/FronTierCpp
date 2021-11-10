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

#include "FabricManager.h"
#include "bending.h"


#define	MAX_NUM_RING1 30


static void compute_total_canopy_force2d(Front*,double*,double*);
static void compute_total_canopy_force3d(Front*,double*,double*);

static void setCurveVelocity(ELASTIC_SET*,CURVE*,double**,GLOBAL_POINT**);
static void setSurfVelocity(ELASTIC_SET*,SURFACE*,double**,GLOBAL_POINT**);
static void setNodeVelocity(ELASTIC_SET*,NODE*,double**,GLOBAL_POINT**);
static void setNodeVelocity2d(ELASTIC_SET*,NODE*,GLOBAL_POINT**);
static void setNodeVelocity3d(ELASTIC_SET*,NODE*,GLOBAL_POINT**);

static void reduce_high_freq_vel(Front*,SURFACE*);
static void smooth_vel(double*,POINT*,TRI*,SURFACE*);

static void compute_center_of_mass_velo(ELASTIC_SET*);

static void break_string_curve(CURVE*,double);
static void linkGlobalIndexToTri(INTERFACE*,TRI***);//nowhere used

static void setCollisionFreePoints3d(INTERFACE*);

static FABRIC_COLLISION_PARAMS getFabricCollisionParams(Front* front);
static void print_elastic_params(const ELASTIC_SET& geom_set);

static void fourth_order_elastic_set_propagate3d_serial(Front*,Front**,double);
static void elastic_set_propagate_serial(Front* fr, Front** newfront, double fr_dt);
static int elastic_set_propagate3d_serial(ELASTIC_SET*, Front*, Front**, int, double);
static void new_fourth_order_elastic_set_propagate3d_parallel_1(Front*,Front**,double);

static void print_max_fabric_speed(Front* fr);
static void print_max_string_speed(Front* fr);

//TODO: UNUSED FUNCTION
extern void compute_total_canopy_force(
	Front *front,
	double *pos_force,
	double *neg_force)
{
	int dim = front->rect_grid->dim;
	switch (dim)
	{
	case 2:
	    compute_total_canopy_force2d(front,pos_force,neg_force);
	    return;
	case 3:
	    compute_total_canopy_force3d(front,pos_force,neg_force);
	    return;
	}
}	/* end compute_total_canopy_force */

//TODO: UNUSED FUNCTION
static void compute_total_canopy_force2d(
	Front *front,
        double *pos_force,
        double *neg_force)
{
    printf("\nWARNING compute_total_canopy_force2d() not implemented\n");
}	/* end compute_total_canopy_force2d */

//TODO: UNUSED FUNCTION
static void compute_total_canopy_force3d(
	Front *front,
        double *pos_force,
        double *neg_force)
{
	TRI *tri;
	SURFACE *surf;
	INTERFACE *intfc = front->interf;
	POINT *p;
	STATE *sl,*sr;
	double pres_p,pres_m;
	double area[MAXD];
	int i;
	static FILE *pfile;

	if (debugging("trace"))
	    (void) printf("Entering compute_total_canopy_force3d()\n");
	if (pfile == NULL)
	{
	    pfile = fopen("payload","w");
	    fprintf(pfile,"\"Net lift vs time\"\n");
	}
	
    for (i = 0; i < 3; ++i)
    {
	    pos_force[i] = 0.0;
        neg_force[i] = 0.0;
    }

    next_tri(intfc,NULL,NULL);
	while (next_tri(intfc,&tri,&surf))
	{
	    if (wave_type(surf) != ELASTIC_BOUNDARY) continue; 
	    
        pres_p = 0.0;
        pres_m = 0.0;
	    
        for (i = 0; i < 3; ++i)
	    {
            p = Point_of_tri(tri)[i];
            FT_GetStatesAtPoint(p,Hyper_surf_element(tri),Hyper_surf(surf),
                    (POINTER*)&sl,(POINTER*)&sr);
            pres_m += sl->pres;
            pres_p += sr->pres;
            area[i] = Tri_normal(tri)[i];//TODO: NO -- see ifluid_compute_torque_and_force()
	    }

        //TODO: This does not look correct below
	    for (i = 0; i < 3; ++i)
	    {
            pos_force[i] -= pres_p*area[i]/3.0;
            neg_force[i] += pres_m*area[i]/3.0;
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving compute_total_canopy_force3d()\n");
}	/* end compute_total_canopy_force3d */

//GFM at ELASTIC_BOUNDARY:
//      Add artificial source terms at boundary and have
//      the fluid solver treat as NO_PDE_BOUNDARY
int af_find_state_at_crossing(
    Front *front,
    int *icoords,
    GRID_DIRECTION dir,
    int comp,
    POINTER *state,
    HYPER_SURF **hs,
    double *crx_coords)
{
    boolean status;
    HYPER_SURF_ELEMENT *hse;
    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

    status = FT_StateStructAtGridCrossing2(front,icoords,
            dir,comp,state,hs,&hse,crx_coords);
    
    if (status == NO)
        return NO_PDE_BOUNDARY;
    else
    {
        switch (wave_type(*hs))
        {
        case FIRST_PHYSICS_WAVE_TYPE:
        case ELASTIC_STRING:
        case ELASTIC_BOUNDARY:
            return NO_PDE_BOUNDARY;
        case MOVABLE_BODY_BOUNDARY:
            return CONST_V_PDE_BOUNDARY;
        case DIRICHLET_BOUNDARY:
	        if (boundary_state_function(*hs) &&
                strcmp(boundary_state_function_name(*hs),
                "flowThroughBoundaryState") == 0)
                return CONST_P_PDE_BOUNDARY;
            else
                return CONST_V_PDE_BOUNDARY;
        /*case ELASTIC_BOUNDARY:
            return POROUS_BOUNDARY;*/
        }
    }
    //TODO: this should return CONST_V_PDE_BOUNDARY,
    //      or MOVABLE_BODY_BOUNDAY should return
    //      NEUMANN_PDE_BOUNDARY in order to be consistent
    return NEUMANN_PDE_BOUNDARY;
}       /* af_find_state_at_crossing */

static void compute_center_of_mass_velo(
	ELASTIC_SET *geom_set)
{
	int i,j,n;
	TRI *tri;
	POINT *p;
	STATE *state;
	Front *front = geom_set->front;
	SURFACE **surfs = geom_set->surfs;
	SURFACE *canopy;
	NODE *node;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double area_dens = af_params->area_dens;
	double xt[MAXD],vt[MAXD],xcan[MAXD],vcan[MAXD],xload[MAXD],vload[MAXD];
	double area,mass_canopy,payload;
	double *xcom,*vcom;

	if (debugging("canopy"))
	    (void) printf("Entering compute_center_of_mass_velo()\n");

	for (n = 0; n < geom_set->num_surfs; ++n)
	{
	    canopy = geom_set->surfs[n];

	    for (j = 0; j < 3; ++j) vcan[j] = 0.0;
	    area = mass_canopy = 0.0;
	    
        surf_tri_loop(canopy,tri)
	    {
	    	for (j = 0; j < 3; ++j)
	    	{
                vt[j] = 0.0;
                xt[j] = 0.0;
	    	}
	    	
            for (i = 0; i < 3; ++i)
	    	{
                p = Point_of_tri(tri)[i];
                state = (STATE*)left_state(p);
		    
                for (j = 0; j < 3; ++j)
                {
                    vt[j] += state->vel[j]/3.0;
                    xt[j] += Coords(p)[j]/3.0;
                }
	    	}

	    	for (j = 0; j < 3; ++j)
	    	{
                vcan[j] += vt[j]*tri_area(tri);
                xcan[j] += xt[j]*tri_area(tri);
	    	}

	    	area += tri_area(tri);
	    }
	    
        mass_canopy += area_dens*area;
	    
        for (j = 0; j < 3; ++j)
	    {
	    	vcan[j] /= area;
	    	xcan[j] /= area;
	    }

        
        //TODO: Add string curves; need to specify linear density for the paracord

        //TODO: Add rigid body center of mass velocity if present


	    if (NULL != geom_set->load_node)
	    {
            node = geom_set->load_node;
            state = (STATE*)left_state(node->posn);
            for (j = 0; j < 3; ++j)
            {
                vload[j] = state->vel[j];
                xload[j] = Coords(node->posn)[j];
            }

            payload = af_params->payload;
            xcom = center_of_mass(Hyper_surf(canopy));
            vcom = center_of_mass_velo(Hyper_surf(canopy));
            for (j = 0; j < 3; ++j)
            {
                vcom[j] = (vcan[j]*mass_canopy + vload[j]*payload)/(mass_canopy + payload);
                xcom[j] = (xcan[j]*mass_canopy + xload[j]*payload)/(mass_canopy + payload);
            }
	    }
	    else
	    {
            xcom = center_of_mass(Hyper_surf(canopy));
            vcom = center_of_mass_velo(Hyper_surf(canopy));
	    }	
	}

	if (debugging("canopy"))
	    (void) printf("Leaving compute_center_of_mass_velo()\n");
}	/* end compute_center_of_mass_velo */

static void smooth_vel(
	double *vel,
	POINT *p,
	TRI *tri,
	SURFACE *surf)
{
	TRI *tris[20];
	HYPER_SURF_ELEMENT *hse = Hyper_surf_element(tri);
        HYPER_SURF         *hs = Hyper_surf(surf);
	int i,j,k,nt,np;
	POINT *pt_list[20],*pt;
	STATE *sl,*sr;
	static double max_speed = 0.0;
	
	PointAndFirstRingTris(p,hse,hs,&nt,tris);
	np = 0;
	
    for (k = 0; k < 3; ++k)
	    vel[k] = 0.0;

	for (i = 0; i < nt; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
            pt = Point_of_tri(tris[i])[j];
            if (pointer_in_list((POINTER)pt,np,(POINTER*)pt_list))
                continue;

            pt_list[np++] = pt;	
            hse = Hyper_surf_element(tris[i]);
            FT_GetStatesAtPoint(pt,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            for (k = 0; k < 3; ++k)
                vel[k] += sl->vel[k];
	    }
	}
	for (k = 0; k < 3; ++k)
	    vel[k] /= (double)np;
}	/* end smooth_vel */

extern boolean is_registered_point(
	SURFACE *surf,
	POINT *p)
{
	REGISTERED_PTS *rgp = (REGISTERED_PTS*)surf->extra;
	int i,num_pts;
	int *global_ids;
	
	if (rgp == NULL) return NO;

	num_pts = rgp->num_pts;
	global_ids = rgp->global_ids;
	for (i = 0; i < num_pts; ++i)
	{
	    if (Gindex(p) == global_ids[i])
		return YES;
	}
	return NO;
}	/* end is_registered_point */

extern void propagate_surface(
        ELASTIC_SET *geom_set,
        SURFACE *surf,
        double **x,
        int *n)
{
        int i,j;
        TRI *tri;
        POINT *p;
        STATE *sl,*sr;
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double *g = eqn_params->gravity;

	hs = Hyper_surf(surf);
	unsort_surf_point(surf);
        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
                        tri = tri->next)
        {
            hse = Hyper_surf_element(tri);
            for (i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                if (sorted(p) || Boundary_point(p)) continue;
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		if (is_registered_point(surf,p))
		{
                    for (j = 0; j < 3; ++j)
                    {
                        x[*n][j] += sl->impulse[j]*dt;
                        sr->impulse[j] = sl->impulse[j] = sl->impulse[j];
                    }
		}
		else
                {
                    for (j = 0; j < 3; ++j)
                    {
                        x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
                        sr->impulse[j] = sl->impulse[j] = 
					sl->impulse[j] + g[j]*dt;
                    }
                }
                sorted(p) = YES;
                ++(*n);
            }
        }
}       /* propagate_surface */

extern void propagate_node(
        ELASTIC_SET *geom_set,
	NODE *node,
        double **x,
        int *n)
{
        int i,j;
        POINT *p;
        STATE *sl,*sr;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double *g = eqn_params->gravity;
	int dim = front->rect_grid->dim;

        sl = (STATE*)left_state(node->posn);
        sr = (STATE*)right_state(node->posn);
        for (j = 0; j < dim; ++j)
        {
            x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
            sr->impulse[j] = sl->impulse[j] = sl->impulse[j] + g[j]*dt;
        }
        ++(*n);
}	/* end propagate_node */

extern void propagate_curve(
        ELASTIC_SET *geom_set,
	CURVE *curve,
        double **x,
        int *n)
{
        int i,j;
        POINT *p;
	BOND *b;
        STATE *sl,*sr;
	double dt = geom_set->dt;
	Front *front = geom_set->front;
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
	double *g = eqn_params->gravity;
	int dim = front->rect_grid->dim;

	for (b = curve->first; b != curve->last; b = b->next)
        {
            p = b->end;
            sl = (STATE*)left_state(p);
            sr = (STATE*)right_state(p);
            for (j = 0; j < dim; ++j)
            {
                x[*n][j] += (sl->impulse[j] + 0.5*g[j]*dt)*dt;
                sr->impulse[j] = sl->impulse[j] = sl->impulse[j] + g[j]*dt;
            }
            ++(*n);
        }
}	/* end propagate_curve */

//TODO: This function needs some refactoring/debugging
static void reduce_high_freq_vel(
	Front *front,
	SURFACE *canopy)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double **vv;
	int ncan;
	POINT *p;
	TRI *tri;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	double crds_max[MAXD];
	int i,j,gindex_max;
	int l,num_layers = af_params->num_smooth_layers;
	
    double max_speed;
    double max_nor_speed = 0.0;

	hs = Hyper_surf(canopy);

	ncan = 0;
	unsort_surf_point(canopy);
	for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p) || Boundary_point(p)) continue;
		sorted(p) = YES;
		ncan++;
	    }
	}

	FT_MatrixMemoryAlloc((POINTER*)&vv,ncan,3,sizeof(double));

	l = 0;
	while (l < num_layers)
	{
	    ncan = 0;
	    unsort_surf_point(canopy);
	    for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	    {
	    	for (i = 0; i < 3; ++i)
	    	{
		    p = Point_of_tri(tri)[i];
		    if (sorted(p) || Boundary_point(p)) continue;
		    smooth_vel(vv[ncan],p,tri,canopy);
		    ncan++;
		    sorted(p) = YES;
	    	}
	    }

	    ncan = 0;
	    max_speed = 0.0;
	    unsort_surf_point(canopy);
	    for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	    {
	    	hse = Hyper_surf_element(tri);
	    	for (i = 0; i < 3; ++i)
	    	{
		    p = Point_of_tri(tri)[i];
		    if (sorted(p) || Boundary_point(p)) continue;
		    
            double nor[MAXD];
            FT_NormalAtPoint(p,front,nor,NO_COMP);
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    double nor_speed = scalar_product(vv[ncan],nor,3);
		    
            for (j = 0; j < 3; ++j)
		    {
		    	sl->vel[j] = nor_speed*nor[j];
		    	sr->vel[j] = nor_speed*nor[j];
		    	    //sl->vel[j] = vv[ncan][j];
		    	    //sr->vel[j] = vv[ncan][j];
		    }

		    if (max_speed < Mag3d(sl->vel)) 
            {
                max_speed = Mag3d(sl->vel);
            }
             
            ncan++;
		    sorted(p) = YES;
	    	}
	    }
	    if (debugging("smooth_canopy_vel"))
		(void) printf("Max speed after smoothing round %d: %f\n",
					l,max_speed);
	    l++;
	}

	unsort_surf_point(canopy);
	max_speed = 0.0;
	for (tri = first_tri(canopy); !at_end_of_tri_list(tri,canopy); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
            p = Point_of_tri(tri)[i];
            if (sorted(p) || Boundary_point(p)) continue;
            
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            
            for (j = 0; j < 3; ++j)
            {
                FT_RecordMaxFrontSpeed(j,sl->vel[j],NULL,Coords(p),front);
            }    
	        
            //need to use normal veocity when setting for dir = dim (= 3)   
            FT_RecordMaxFrontSpeed(3,Mag3d(sl->vel),NULL,Coords(p),front);
		
        if (max_speed < Mag3d(sl->vel)) 
		{
		    max_speed = Mag3d(sl->vel);
		    gindex_max = Gindex(p);
		    for (j = 0; j < 3; ++j)
			    crds_max[j] = Coords(p)[j];
		}

		sorted(p) = YES;
	    }
	}
	FT_FreeThese(1,vv);
}	/* end reduce_high_freq_vel */

static void print_elastic_params(
	const ELASTIC_SET& geom_set)
{
	int i;
	double *spfr;
	Front *fr = geom_set.front;

	spfr = Spfr(fr);
    for (i = 0; i <= 3; ++i)
    {
        printf("Max front speed(%d) = %f\n",i,spfr[i]);
    }
    
    (void) printf("Input surface parameters:\n");
    (void) printf("ks = %f  m_s = %f  lambda_s = %f\n",
                    geom_set.ks,
                    geom_set.m_s,
                    geom_set.lambda_s);
    
    (void) printf("Input string parameters:\n");
    (void) printf("kl = %f  m_l = %f  lambda_l = %f\n",
                    geom_set.kl,
                    geom_set.m_l,
                    geom_set.lambda_l);
    
    (void) printf("Input gore parameters:\n");
    (void) printf("kg = %f  m_g = %f  lambda_g = %f\n",
                    geom_set.kg,
                    geom_set.m_g,
                    geom_set.lambda_g);
	
    (void) printf("\ndt_tol = %20.14f  dt = %20.14f\n",
                        geom_set.dt_tol,geom_set.dt);
}	/* end print_elastic_params */

static FABRIC_COLLISION_PARAMS getFabricCollisionParams(Front* front)
{
    FABRIC_COLLISION_PARAMS collsn_params;
    AF_PARAMS* af_params = (AF_PARAMS*)front->extra2;

    collsn_params.fabric_eps = af_params->fabric_eps;
    collsn_params.fabric_thickness = af_params->fabric_thickness;
    collsn_params.mu_s = af_params->mu_s;
    collsn_params.k_s = af_params->ks;
    collsn_params.m_s = af_params->m_s;

    collsn_params.string_eps = af_params->string_eps;
    collsn_params.string_thickness = af_params->string_thickness;
    collsn_params.mu_l = af_params->mu_l;
    collsn_params.k_l = af_params->kl;
    collsn_params.m_l = af_params->m_l;

    collsn_params.overlap_coefficient = af_params->overlap_coefficient;

    collsn_params.strain_limit = af_params->strain_limit;
    collsn_params.compressive_strain_limit = af_params->compressive_strain_limit;
    collsn_params.strainrate_limit = af_params->strainrate_limit;

    collsn_params.coefRestitution = 1.0;
            
    collsn_params.collision_off = false;
    if (debugging("collision_off")) //TODO: don't use debug string -- add input file option
    {
        collsn_params.collision_off = true;
    }

    return collsn_params;
}

void fourth_order_elastic_set_propagate(Front* fr, double fr_dt)
{
    Front* newfront;

    if (pp_numnodes() > 1 && !debugging("collision_off"))
    {
        new_fourth_order_elastic_set_propagate3d_parallel_1(fr,&newfront,fr_dt);
            //fourth_order_elastic_set_propagate_parallel(fr,fr_dt);
    }
    else
    {
        fourth_order_elastic_set_propagate3d_serial(fr,&newfront,fr_dt);
            //fourth_order_elastic_set_propagate_serial(fr,fr_dt);
    }

    assign_interface_and_free_front(fr,newfront);
}

//NEW
void elastic_set_propagate(Front* fr, double fr_dt)
{
    Front* newfront;
   
    if (pp_numnodes() > 1 && !debugging("collision_off"))
    {
        printf("\nERROR elastic_set_propagate(): \
                parallel code not yet implemented!\n");
        LOC(); clean_up(EXIT_FAILURE);
        //new_fourth_order_elastic_set_propagate3d_parallel_1(fr,&newfront,fr_dt);
            //fourth_order_elastic_set_propagate_parallel(fr,fr_dt);
    }
    else
    {
        elastic_set_propagate_serial(fr,&newfront,fr_dt);
    }

    assign_interface_and_free_front(fr,newfront);
} /*end elastic_set_propagate*/

static void elastic_set_propagate_serial(
        Front* fr,
        Front** newfront,
        double fr_dt)
{
	static ELASTIC_SET geom_set;
    AF_PARAMS* af_params = (AF_PARAMS*)fr->extra2;
    
    static boolean first = YES;
	if (first)
    {
        geom_set.front = fr;
        set_elastic_params(&geom_set,fr_dt);
        if (debugging("step_size"))
            print_elastic_params(geom_set);
    }

    int n_sub = af_params->n_sub;
    double ss_dt = fr_dt/((double)n_sub);

    double ss_dt_tol = af_params->ss_dt_relax*geom_set.dt_tol;
    if (fr_dt > ss_dt_tol)
    {
        n_sub = (int)(fr_dt/ss_dt_tol);
        ss_dt = fr_dt/((double)n_sub);
    }
    geom_set.n_sub = n_sub;
    geom_set.dt = ss_dt;
    
    printf("fr_dt = %f ss_dt_tol = %f n_sub = %d ss_dt = %f\n",
            fr_dt, ss_dt_tol, n_sub, ss_dt);

    int nsub_per_collsn_step = af_params->collsn_step_max_nsub;
    int collsn_nsub = nsub_per_collsn_step;
    
    int num_collsn_steps = n_sub/nsub_per_collsn_step;
    int remainder_num_steps = n_sub % nsub_per_collsn_step;
    if (num_collsn_steps == 0)
    {
        num_collsn_steps = 1;
        collsn_nsub = n_sub;
        remainder_num_steps = 0;
    }
    double collsn_dt = collsn_nsub*ss_dt;
    printf("num_collsn_steps = %d\n",num_collsn_steps);
    
    *newfront = copy_front(fr);
    set_size_of_intfc_state(size_of_state(fr->interf));
    set_copy_intfc_states(YES);
    (*newfront)->interf = pp_copy_interface(fr->interf);
    if ((*newfront)->interf == NULL)
    {
        (void) printf("WARNING in fourth_order_elastic_set_propagate3d_serial(), "
                      "unable to copy interface\n");
        LOC(); clean_up(EXIT_FAILURE);
            //return return_advance_front(fr,newfront,ERROR_IN_STEP,fname);
    }


    //TODO: Need the FabricManager here.
    //      Perform initialization using newfront...
    


    //TODO: In the collsn substepping loop below we advance forward in time
    //      using the spring solver along with applying repulsions, friction
    //      and strain limiting.

    int status;
    for (int n = 1; n <= num_collsn_steps; ++n)
    {
        int cd_nsub = collsn_nsub;
        
        //adjust collsn_dt for first collision substep that may
        //have a different number of n_sub (spring solver) substeps
        if (n == 1 && remainder_num_steps > 0)
        {
            cd_nsub = collsn_nsub + remainder_num_steps;
        }
            
        collsn_dt = cd_nsub*ss_dt;

        printf("\nCollision Sub Step #%d: ",n);
        printf("\tcollsn_nsub = %d  collsn_dt = %f\n\n",cd_nsub,collsn_dt);

        Front* sub_newfront;
        
        //NOTE: ignoring status for now since using fixed nsub_per_collsn_step
        //      instead of time step modification.
        status = elastic_set_propagate3d_serial(
                &geom_set,
                *newfront,
                &sub_newfront,
                cd_nsub,
                collsn_dt);
    
        assign_interface_and_free_front(*newfront,sub_newfront);
    }


    //TODO: Should check the interface for intersections/crossing
    //      here and attempt to correct before returning to the
    //      calling function.
    //
    //      Prevent remove_unphys_pair() from failing and crashing the run.
    //
        //CROSS  *cross;
        //boolean intersection_status;
        //intersection_status = intersections(newfront->interf,&cross,YES);

    
    //      The velocity constraint of strain limiting procedure could
    //      be applied here after the vertices have been moved to their
    //      final collision free position.


    //setSpecialNodeForce((*newfront)->interf,geom_set.kl);
    //compute_center_of_mass_velo(&geom_set);

} /*end elastic_set_propagate_serial*/

static int elastic_set_propagate3d_serial(
        ELASTIC_SET* geom_set,
        Front* fr,
        Front** newfront,
        int collsn_nsub,
        double collsn_dt)
{
    int myid = pp_mynode();

    static const char *fname = "elastic_set_propagate3d_serial";
	static int total_size = 0;
    static int size = 0;
    static int total_owner_size, owner_size, total_client_size, client_size;
	static int *client_size_old, *client_size_new;

    AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	int owner_id = af_params->node_id[0];
    int dim = FT_Dimension();

    long max_point_gindex = fr->interf->max_point_gindex;
	int gindex;
    
	double *L = fr->rect_grid->L;
	double *U = fr->rect_grid->U;
	double *h = fr->rect_grid->h;
	double client_L[MAXD],client_U[MAXD];

    static boolean first_break_strings = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;
        
	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate3d_serial()\n");

    int status; //return status

    //NOTE: newfront is actually the sub_newfront and fr is actually the newfront
    *newfront = copy_front(fr);
    set_size_of_intfc_state(size_of_state(fr->interf));
    set_copy_intfc_states(YES); //Why set to NO in struct_advance_front3d() ???
    (*newfront)->interf = pp_copy_interface(fr->interf);
    if ((*newfront)->interf == NULL)
    {
        (void) printf("WARNING in fourth_order_elastic_set_propagate3d_serial(), "
                      "unable to copy interface\n");
        LOC(); clean_up(EXIT_FAILURE);
            //return return_advance_front(fr,newfront,ERROR_IN_STEP,fname);
    }


    //TODO: should call this function, in case we need the previous interface states
    //
    //  set_correspondence_between_interfaces((*newfront)->interf,fr->interf);
    
    geom_set->front = *newfront;


    static SPRING_VERTEX* sv;
    static GLOBAL_POINT** point_set;
    static GLOBAL_POINT* point_set_store;
	static GLOBAL_POINT** client_point_set_store;
    
    INTERFACE *elastic_intfc = nullptr;

    ///////////////////////////////////////////////////////////////////////////////////////////
    //TODO: Move this into separate intitialization function that is called
    //      in the level above this function (in its calling function)/
    static boolean first = YES;
	if (first)
	{

	    int owner[MAXD];
        owner[0] = 0;
        owner[1] = 0;
        owner[2] = 0;
	    
        if (point_set != NULL)
		    FT_FreeThese(1, point_set);

	    FT_VectorMemoryAlloc((POINTER*)&point_set,max_point_gindex,sizeof(GLOBAL_POINT*));
	    
        for (int i = 0; i < max_point_gindex; ++i)
		    point_set[i] = NULL;

	    if (pp_numnodes() > 1)
	    {
            int w_type[3] = {ELASTIC_BOUNDARY,MOVABLE_BODY_BOUNDARY,NEUMANN_BOUNDARY};
            elastic_intfc = collect_hyper_surfaces(*newfront,owner,w_type,3);
            collectNodeExtra(*newfront,elastic_intfc,owner_id);
	    }
	    else
        {
            elastic_intfc = (*newfront)->interf;
        }
	    
        start_clock("set_data");
        if (myid == owner_id)
        {

            if (client_size_old != NULL)
                FT_FreeThese(3, client_size_old, client_size_new, client_point_set_store);

            FT_VectorMemoryAlloc((POINTER*)&client_size_old,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_size_new,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_point_set_store,
                    pp_numnodes(),sizeof(GLOBAL_POINT*));

            for (int i = 0; i < pp_numnodes(); i++)
            {
                client_size_old[i] = 0;
                client_size_new[i] = 0;
            }

            assembleParachuteSet(elastic_intfc,geom_set);
            
            total_owner_size = geom_set->total_num_verts; //fabric + rgb points
            owner_size = geom_set->elastic_num_verts; //just fabric points
            
            if (point_set_store != NULL) 
                FT_FreeThese(2,point_set_store, sv);
            
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,
                    total_owner_size,sizeof(GLOBAL_POINT)); //NOTE: uses total_num_verts here
            FT_VectorMemoryAlloc((POINTER*)&sv,
                    owner_size,sizeof(SPRING_VERTEX)); //NOTE: uses elastic_num_verts here
            
            link_point_set(geom_set,point_set,point_set_store);
	    	count_vertex_neighbors(geom_set,sv);
	    	set_spring_vertex_memory(sv,owner_size);
	    	set_vertex_neighbors(geom_set,sv,point_set);
		
            if (elastic_intfc != (*newfront)->interf)
            {
                delete_interface(elastic_intfc);
            }
	    }
	    
        stop_clock("set_data");
	    first = NO;
	}
    //end initialization of point_set and sv array etc.
    ///////////////////////////////////////////////////////////////////////////////////////////

	elastic_intfc = (*newfront)->interf;

    setCollisionFreePoints3d(elastic_intfc);
    
    resetBendingForce(elastic_intfc);
    if (!debugging("bendforce_off"))
    {
        double bends = af_params->kbs;
        double bendd = af_params->lambda_bs;
        computeSurfBendingForce(elastic_intfc,bends,bendd);//TODO: make function monadic
        computeStringBendingForce(elastic_intfc);
    }

	assembleParachuteSet(elastic_intfc,geom_set);
	
    if (myid != owner_id)
	{
	    total_client_size = geom_set->total_num_verts; //fabric + rigid body points
        client_size = geom_set->elastic_num_verts; //just fabric points
	    if (size < client_size || total_size < total_client_size)
	    {
	    	size = client_size;
	    	total_size = total_client_size;
	    	if (point_set_store != NULL)
            {
	    	    FT_FreeThese(2,point_set_store,sv);
            }
	
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,total_size,sizeof(GLOBAL_POINT));
            FT_VectorMemoryAlloc((POINTER*)&sv,size,sizeof(SPRING_VERTEX));
	    }

	    for (int i = 0; i < max_point_gindex; ++i)
        {
            point_set[i] = nullptr;
        }
	    
        link_point_set(geom_set,point_set,point_set_store);
	    count_vertex_neighbors(geom_set,sv);
	    set_spring_vertex_memory(sv,client_size);
	    set_vertex_neighbors(geom_set,sv,point_set);
	    get_point_set_from(geom_set,point_set);
	    
        pp_send(5,L,MAXD*sizeof(double),owner_id);
	    pp_send(6,U,MAXD*sizeof(double),owner_id);
	    pp_send(1,&(total_client_size),sizeof(int),owner_id);
        pp_send(2,point_set_store,total_client_size*sizeof(GLOBAL_POINT),owner_id);
	}
	else
    {
	    size = owner_size;
        total_size = total_owner_size;
    }

    
    //TODO: FabricManager should be passed into this function, or the function
    //      should be a member function of FabricManager.
    FabricManager fabric_manager(*newfront);
    FABRIC_COLLISION_PARAMS collsn_params = getFabricCollisionParams(*newfront);
    fabric_manager.setCollisionParams(collsn_params);
    fabric_manager.setCollisionTimeStep(collsn_dt);

    if (myid == owner_id)
	{
        fabric_manager.initializeSystem();

        //write to GLOBAL_POINT** point_set
        get_point_set_from(geom_set,point_set);
        
	    for (int i = 0; i < pp_numnodes(); i++)
	    {
            if (i == myid) continue;
            
            pp_recv(5,i,client_L,MAXD*sizeof(double));
            pp_recv(6,i,client_U,MAXD*sizeof(double));
            pp_recv(1,i,client_size_new + i,sizeof(int));
            
            if (client_size_new[i] > client_size_old[i])
            {
                client_size_old[i] = client_size_new[i];
                if (client_point_set_store[i] != NULL)
                    FT_FreeThese(1,client_point_set_store[i]);

                FT_VectorMemoryAlloc((POINTER*)&client_point_set_store[i],
                        client_size_new[i],sizeof(GLOBAL_POINT));
            }

            pp_recv(2,i,client_point_set_store[i],
                    client_size_new[i]*sizeof(GLOBAL_POINT));
            
            copy_from_client_point_set(point_set,client_point_set_store[i],
                    client_size_new[i],client_L,client_U);
	    }
    
        
        //compute candidate positions and velocities with spring solver
        double ss_dt = geom_set->dt;
        spring_solver_RK4(sv,dim,size,collsn_nsub,ss_dt);

        
        // Owner send and patch point_set_store from other processors
	    for (int i = 0; i < pp_numnodes(); i++)
        {
            if (i == myid) continue;
            copy_to_client_point_set(point_set, client_point_set_store[i], client_size_new[i]);
            pp_send(3, client_point_set_store[i], client_size_new[i]*sizeof(GLOBAL_POINT),i);
        }
	}//owner_id


    if (myid != owner_id)
    {
        pp_recv(3,owner_id,point_set_store,total_client_size*sizeof(GLOBAL_POINT));
    }
	
    //  write from point_set to geom_set
    put_point_set_to(geom_set,point_set);
	
    
    //TODO: What about the weight of the mass -- not accounted for in setSpecialNodeForce() ???
    //
    // Calculate the real force on load_node and rg_string_node
    setSpecialNodeForce(elastic_intfc,geom_set->kl);

	set_vertex_impulse(geom_set,point_set);
	set_geomset_velocity(geom_set,point_set);
	
        //compute_center_of_mass_velo(geom_set);

    if (myid == owner_id)
    {
        //TODO: Need resolveCollisionSubstep() to return an error code.
        //      Temporarily have it emit an exception instead for
        //      now while still prototyping new code.
        if (FT_Dimension() == 3)
        {
            try
            {
                fabric_manager.resolveCollisionSubstep();
            }
            catch (...)
            {
                free_front(*newfront);
                *newfront = nullptr;
                
                //TODO: Should we clear rest of geom_set also?
                //      could write the below function to so:
                //
                //          clear_geom_set(geom_set);
                geom_set->front = nullptr;
            
                //TEMPORARY: until we have fixed nsub_per_collsn_step
                //           working correctly, then we can try time step
                //           modification...
                //
                //TODO: Or compute intersection???
                printf("\nERROR elastic_set_propagate3d_serial(): \
                        resolve_collision() failed!\n");
                LOC(); clean_up(EXIT_FAILURE);
                    //status = MODIFY_TIME_STEP;
                    //return status;
            }
        }
    }
    
    setSpecialNodeForce(elastic_intfc,geom_set->kl);
        //compute_center_of_mass_velo(geom_set);


    //TODO: Sync interfaces after collision handling?
    //      Or call in interior_propagate()?
        //scatter_front(*newfront);


    if (debugging("max_speed"))
    {
        print_max_fabric_speed(*newfront);
        print_max_string_speed(*newfront);
        fflush(stdout);
    }

	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate3d_serial()\n");

    status = GOOD_STEP;
    return status;
        //return return_advance_front(fr,newfront,status,fname);

}	/* end elastic_set_propagate3d_serial() */

static void fourth_order_elastic_set_propagate3d_serial(
        Front* fr,
        Front** newfront,
        double fr_dt)
{
	static ELASTIC_SET geom_set;
	static int total_size = 0;
    static int size = 0;
    static int total_owner_size, owner_size, total_client_size, client_size;
	static int *client_size_old, *client_size_new;
    static SPRING_VERTEX *sv;
    static boolean first = YES;
    static GLOBAL_POINT **point_set;
    static GLOBAL_POINT *point_set_store;
	static GLOBAL_POINT **client_point_set_store;
    
    AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
    int i,j,k,n_sub;
    double dt;
    int dim = FT_Dimension();
    long max_point_gindex = fr->interf->max_point_gindex;
	int owner[MAXD];
	int owner_id = af_params->node_id[0];
    int myid = pp_mynode();
	int gindex;
    INTERFACE *elastic_intfc = nullptr;
	double *L = fr->rect_grid->L;
	double *U = fr->rect_grid->U;
	double *h = fr->rect_grid->h;
	double client_L[MAXD],client_U[MAXD];

    static boolean first_break_strings = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;
	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate3d_serial()\n");

    *newfront = copy_front(fr);
    set_size_of_intfc_state(size_of_state(fr->interf));
    set_copy_intfc_states(YES);
        //set_copy_intfc_states(NO); //TODO: Why was this set to NO?
    (*newfront)->interf = pp_copy_interface(fr->interf);
    if ((*newfront)->interf == NULL)
    {
        (void) printf("WARNING in fourth_order_elastic_set_propagate3d_serial(), "
                      "unable to copy interface\n");
        LOC(); clean_up(EXIT_FAILURE);
            //return return_advance_front(fr,newfront,ERROR_IN_STEP,fname);
    }

    geom_set.front = *newfront;

    //TODO: omitting for now -- needs to be tested.
	//if (first_break_strings && break_strings_num > 0 &&
	//    break_strings_time >= 0.0 && 
	//    fr->time + fr->dt >= break_strings_time)
	//{
	//    printf("Some strings break! Count and set spring vertex again.\n");
	//    first_break_strings = NO;
	//    first = YES;
	//}

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
    printf("fr_dt = %f geom_set.dt_tol = %f n_sub = %d dt = %f\n",
            fr_dt,geom_set.dt_tol,n_sub,dt);

	if (first)
	{
        owner[0] = 0;
        owner[1] = 0;
        owner[2] = 0;
	    
        if (point_set != NULL)
		    FT_FreeThese(1, point_set);

	    FT_VectorMemoryAlloc((POINTER*)&point_set,max_point_gindex,sizeof(GLOBAL_POINT*));
	    
        for (i = 0; i < max_point_gindex; ++i)
		    point_set[i] = NULL;

	    if (pp_numnodes() > 1)
	    {
            int w_type[3] = {ELASTIC_BOUNDARY,MOVABLE_BODY_BOUNDARY,NEUMANN_BOUNDARY};
            elastic_intfc = collect_hyper_surfaces(*newfront,owner,w_type,3);
            collectNodeExtra(*newfront,elastic_intfc,owner_id);
	    }
	    else
        {
            elastic_intfc = (*newfront)->interf;
        }
	    
        start_clock("set_data");
        if (myid == owner_id)
        {

            if (client_size_old != NULL)
                FT_FreeThese(3, client_size_old, client_size_new, client_point_set_store);

            FT_VectorMemoryAlloc((POINTER*)&client_size_old,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_size_new,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_point_set_store,
                    pp_numnodes(),sizeof(GLOBAL_POINT*));

            for (i = 0; i < pp_numnodes(); i++)
                client_size_old[i] = client_size_new[i] = 0;

            assembleParachuteSet(elastic_intfc,&geom_set);
            
            total_owner_size = geom_set.total_num_verts; //fabric + rgb points
            owner_size = geom_set.elastic_num_verts; //just fabric points
            
            if (point_set_store != NULL) 
                FT_FreeThese(2,point_set_store, sv);
            
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,
                    total_owner_size,sizeof(GLOBAL_POINT)); //NOTE: uses total_num_verts here
            FT_VectorMemoryAlloc((POINTER*)&sv,
                    owner_size,sizeof(SPRING_VERTEX)); //NOTE: uses elastic_num_verts here
            
            link_point_set(&geom_set,point_set,point_set_store);
	    	count_vertex_neighbors(&geom_set,sv);
	    	set_spring_vertex_memory(sv,owner_size);
	    	set_vertex_neighbors(&geom_set,sv,point_set);
		
            if (elastic_intfc != (*newfront)->interf)
                delete_interface(elastic_intfc);
	    }
	    
        stop_clock("set_data");
	    first = NO;
	}

	elastic_intfc = (*newfront)->interf;

    //TODO: Do we need to call setCollisionFreePoints3d() prior to the
    //      bending force computations? It appears that we may be using
    //      lagged values for the STATE::is_fixed etc. boolean flags from
    //      the previous time step. 
    setCollisionFreePoints3d(elastic_intfc);
    
    if (!debugging("bendforce_off"))
    {
        resetBendingForce(elastic_intfc);
        double bends = af_params->kbs;
        double bendd = af_params->lambda_bs;
        computeSurfBendingForce(elastic_intfc,bends,bendd);//TODO: make function monadic
        computeStringBendingForce(elastic_intfc);
    }

	assembleParachuteSet(elastic_intfc,&geom_set);
	
    if (myid != owner_id)
	{
	    total_client_size = geom_set.total_num_verts; //fabric + rigid body points
        client_size = geom_set.elastic_num_verts; //just fabric points
	    if (size < client_size || total_size < total_client_size)
	    {
	    	size = client_size;
	    	total_size = total_client_size;
	    	if (point_set_store != NULL)
            {
	    	    FT_FreeThese(2,point_set_store,sv);
            }
	
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,total_size,sizeof(GLOBAL_POINT));
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
	    pp_send(1,&(total_client_size),sizeof(int),owner_id);
        pp_send(2,point_set_store,total_client_size*sizeof(GLOBAL_POINT),owner_id);
	}
	else
    {
	    size = owner_size;
        total_size = total_owner_size;
    }

    CollisionSolver3d* collision_solver = nullptr;

    if (myid == owner_id)
	{
        if (!debugging("collision_off"))
        {
            collision_solver = new CollisionSolver3d();
            printf("COLLISION DETECTION ON\n");
            
            //TODO: Moved setCollisionFreePoints3d() up to be called
            //      before the bending force computations.
            //
            //setCollisionFreePoints3d((*newfront)->interf);

            collision_solver->initializeSystem(*newfront);
        
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

            collision_solver->gpoints = (*newfront)->gpoints;
            collision_solver->gtris = (*newfront)->gtris;
        }
        else
        {
            printf("COLLISION DETECTION OFF\n");
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
                        client_size_new[i],sizeof(GLOBAL_POINT));
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

	    // Owner send and patch point_set_store from other processors
	    for (i = 0; i < pp_numnodes(); i++)
        {
            if (i == myid) continue;
            copy_to_client_point_set(point_set,
                        client_point_set_store[i],client_size_new[i]);
            pp_send(3,client_point_set_store[i],
                        client_size_new[i]*sizeof(GLOBAL_POINT),i);
        }
	}

    if (myid != owner_id)
    {
        pp_recv(3,owner_id,point_set_store,total_client_size*sizeof(GLOBAL_POINT));
            //pp_recv(3,owner_id,point_set_store,client_size*sizeof(GLOBAL_POINT));
    }
	
    //  write from point_set to geom_set
    put_point_set_to(&geom_set,point_set);
	
    // Calculate the real force on load_node and rg_string_node
    setSpecialNodeForce(elastic_intfc,geom_set.kl);

	set_vertex_impulse(&geom_set,point_set);
	set_geomset_velocity(&geom_set,point_set);
	compute_center_of_mass_velo(&geom_set);

	if (!debugging("collision_off"))
    {
        if (myid == owner_id)
        {
            if (FT_Dimension() == 3)
                collision_solver->resolveCollisionSubstep();
            delete collision_solver;
        }
        
        setSpecialNodeForce(elastic_intfc,geom_set.kl);
        compute_center_of_mass_velo(&geom_set);
    }
    

    //sync interfaces after collision handling?
        //scatter_front(*newfront);


    if (debugging("max_speed"))
    {
        print_max_fabric_speed(*newfront);
        print_max_string_speed(*newfront);
    }

	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate3d_serial()\n");
}


/*
void fourth_order_elastic_set_propagate_serial(Front* fr, double fr_dt)
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
	double client_L[MAXD],client_U[MAXD];
	static boolean first_break_strings = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;


	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate_serial()\n");
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

	if (first)
	{
        owner[0] = 0;
        owner[1] = 0;
        owner[2] = 0;
	    
        if (point_set != NULL)
            FT_FreeThese(1, point_set);
	    
        FT_VectorMemoryAlloc((POINTER*)&point_set,max_point_gindex,
					sizeof(GLOBAL_POINT*));
	    
        // every processor has a point_set with global indexing
        for (i = 0; i < max_point_gindex; ++i)
            point_set[i] = NULL;

	    if (pp_numnodes() > 1)
	    {
            elastic_intfc = FT_CollectHypersurfFromSubdomains(fr,owner,
                    ELASTIC_BOUNDARY);
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

            FT_VectorMemoryAlloc((POINTER*)&client_size_old,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_size_new,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_point_set_store,
                    pp_numnodes(),sizeof(GLOBAL_POINT*));

            for (i = 0; i < pp_numnodes(); i++)
                client_size_old[i] = client_size_new[i] = 0;

            assembleParachuteSet(elastic_intfc,&geom_set);
            owner_size = geom_set.num_verts;
            
            if (point_set_store != NULL) 
                FT_FreeThese(2,point_set_store, sv);
            
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,
                    owner_size,sizeof(GLOBAL_POINT));
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

    //compute bending force
    double bends = af_params->kbs;
    double bendd = af_params->lambda_bs;
    computeBendingForce(elastic_intfc,bends,bendd);

	assembleParachuteSet(elastic_intfc,&geom_set);

	if (myid != owner_id)
	{
	    client_size = geom_set.num_verts;
	    if (size < client_size)
	    {
	    	size = client_size;
	    	if (point_set_store != NULL)
                FT_FreeThese(2,point_set_store,sv);

            FT_VectorMemoryAlloc((POINTER*)&point_set_store,
                    size,sizeof(GLOBAL_POINT));
            FT_VectorMemoryAlloc((POINTER*)&sv,size,
                    sizeof(SPRING_VERTEX));
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


    CollisionSolver3d* collision_solver = nullptr;

	if (myid == owner_id)
	{
	
    if (!debugging("collision_off"))
    {
        printf("COLLISION DETECTION ON\n");
        collision_solver = new CollisionSolver3d();

	    if (FT_Dimension() == 3)
        {
            setCollisionFreePoints3d(fr->interf);
            collision_solver->initializeSystem(fr);

            collision_solver->setRestitutionCoef(1.0);
            collision_solver->setVolumeDiff(0.0);

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
            //collision_solver->setCompressiveStrainLimit(af_params->compress_strain_limit);
            collision_solver->setStrainRateLimit(af_params->strainrate_limit);

            collision_solver->gpoints = fr->gpoints;
            collision_solver->gtris = fr->gtris;
        }
    }
    else
    {
        printf("COLLISION DETECTION OFF\n");
    }

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

        // owner writes and sends back the updated client point_set_stores
	    for (i = 0; i < pp_numnodes(); i++)
        {
            if (i == myid) continue;
            copy_to_client_point_set(point_set,
                    client_point_set_store[i], client_size_new[i]);
            pp_send(3,client_point_set_store[i],
                            client_size_new[i]*sizeof(GLOBAL_POINT),i);
        }

	}//myid == owner_id

    // clients receive the updated point_set_stores
    if (myid != owner_id)
    {
        pp_recv(3,owner_id,point_set_store,
            client_size*sizeof(GLOBAL_POINT));
    }

    // all processes write from their point_set to their geom_set
	put_point_set_to(&geom_set,point_set);
	
    //TODO: Does this need to be called before the collision
    //      solver has completed?
    //      Also, should be passing in the geom_set instead of
    //      looping over all nodes of fr->interf.
    //      ELASTIC_SET likely did not consider rigid bodies
    //      when it was created -- should make a new struct
    //      PARACHUTE_SET that does and leave ELASTIC_SET for
    //      pointmass runs and for debugging.
	setSpecialNodeForce(fr,geom_set.kl);

    //TODO: I think this function can be eliminated,
    //      see comments inside put_point_value_to()
    //      for required modification.
	set_vertex_impulse(&geom_set,point_set);
    
    //TODO: why only normal component of velocities retained?
    //      see set_geomset_velocity()
	set_geomset_velocity(&geom_set,point_set);
	
    //TODO: verify com_velo computation
    compute_center_of_mass_velo(&geom_set);

	if(!debugging("collision_off"))
    {
        if (myid == owner_id)
        {
            if (FT_Dimension() == 3)
            {
                start_clock("resolveCollision");
                collision_solver->resolveCollisionSubstep();
                stop_clock("resolveCollision");
            }
            delete collision_solver;
        }

        // Calculate the real force on load_node and rg_string_node
        setSpecialNodeForce(fr,geom_set.kl);
	    compute_center_of_mass_velo(&geom_set);
    }

    if (debugging("max_speed"))
    {
        print_max_fabric_speed(fr);
        print_max_string_speed(fr);
    }

	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate_serial()\n");
}*/	/* end fourth_order_elastic_set_propagate_serial() */


//TODO: memory/synchronization bugs
//
//uses a collision solver on each processor
void new_fourth_order_elastic_set_propagate3d_parallel_1(
        Front* fr,
        Front** newfront,
        double fr_dt)
{
	static ELASTIC_SET geom_set;
	static int total_size = 0;
    static int size = 0;
    static int total_owner_size, owner_size, total_client_size, client_size;
	static int *client_size_old, *client_size_new;
    static SPRING_VERTEX *sv;
    static boolean first = YES;
    static GLOBAL_POINT **point_set;
    static GLOBAL_POINT *point_set_store;
	static GLOBAL_POINT **client_point_set_store;

    AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
    int i,j,k,n_sub;
    double dt;
    int dim = FT_Dimension();
    long max_point_gindex = fr->interf->max_point_gindex;
	int owner[MAXD];
	int owner_id = af_params->node_id[0];
    int myid = pp_mynode();
	int gindex;
    
    INTERFACE *elastic_intfc = NULL;
	double *L = fr->rect_grid->L;
	double *U = fr->rect_grid->U;
	double client_L[MAXD],client_U[MAXD];
	
    static boolean first_break_strings = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;

    static const char *fname = "fourth_order_elastic_set_propagate3d_parallel_1";

	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate()\n");
	
    *newfront = copy_front(fr);
    set_size_of_intfc_state(size_of_state(fr->interf));
    set_copy_intfc_states(YES);
        //set_copy_intfc_states(NO); //TODO: Why was this set to NO?
    (*newfront)->interf = pp_copy_interface(fr->interf);
    if ((*newfront)->interf == NULL)
    {
        (void) printf("WARNING in fourth_order_elastic_set_propagate3d_parallel_1(), "
                      "unable to copy interface\n");
        LOC(); clean_up(EXIT_FAILURE);
            //return return_advance_front(fr,newfront,ERROR_IN_STEP,fname);
    }

    geom_set.front = *newfront;

    /*
    //TODO: omitting for now -- needs to be tested.
	if (first_break_strings && break_strings_num > 0 &&
	    break_strings_time >= 0.0 && 
	    fr->time + fr->dt >= break_strings_time)
	{
	    printf("Some strings break! Count and set spring vertex again.\n");
	    first_break_strings = NO;
	    first = YES;
	}
    */

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
            int w_type[3] = {ELASTIC_BOUNDARY,MOVABLE_BODY_BOUNDARY,NEUMANN_BOUNDARY};
            elastic_intfc = collect_hyper_surfaces(*newfront,owner,w_type,3);
            collectNodeExtra(*newfront,elastic_intfc,owner_id);
	    }
	    else
        {
            elastic_intfc = (*newfront)->interf;
        }

        start_clock("set_data");
	    if (myid == owner_id)
        {
            if (client_size_old != NULL)
            {
                FT_FreeThese(3, client_size_old, client_size_new, 
                        client_point_set_store);
            }

            FT_VectorMemoryAlloc((POINTER*)&client_size_old,pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_size_new,pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_point_set_store,pp_numnodes(),sizeof(GLOBAL_POINT*));
            
            for (i = 0; i < pp_numnodes(); i++)
            {
                client_size_old[i] = 0;
                client_size_new[i] = 0;
            }

            assembleParachuteSet(elastic_intfc,&geom_set);

            owner_size = geom_set.elastic_num_verts; //just fabric points
            total_owner_size = geom_set.total_num_verts; //fabric + rgb points
            
            if (point_set_store != NULL) 
            {
                FT_FreeThese(2,point_set_store, sv);
            }

            FT_VectorMemoryAlloc((POINTER*)&sv,owner_size,sizeof(SPRING_VERTEX));
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,total_owner_size,sizeof(GLOBAL_POINT));
            
            //allocate mem for point_set via point_set_store, and store
            //gindex values of elastic_intfc points in array positions
            link_point_set(&geom_set,point_set,point_set_store);

            count_vertex_neighbors(&geom_set,sv);
            set_spring_vertex_memory(sv,owner_size);

            //links sv to point_set
            set_vertex_neighbors(&geom_set,sv,point_set);

            if (elastic_intfc != (*newfront)->interf)
                delete_interface(elastic_intfc);
        }

	    stop_clock("set_data");
	    first = NO;
	}

	
    elastic_intfc = (*newfront)->interf;
    
    //compute bending force
    double bends = af_params->kbs;
    double bendd = af_params->lambda_bs;
    computeSurfBendingForce(elastic_intfc,bends,bendd);
        //computeStringBendingForce(elastic_intfc); //TODO: finish implementation

    assembleParachuteSet(elastic_intfc,&geom_set);


	if (myid != owner_id)
	{
        client_size = geom_set.elastic_num_verts;
        total_client_size = geom_set.total_num_verts;

	    if (size < client_size || total_size < total_client_size)
	    {
	    	size = client_size;
	    	total_size = total_client_size;

	    	if (point_set_store != NULL)
            {
                FT_FreeThese(2,point_set_store,sv);
            }

            FT_VectorMemoryAlloc((POINTER*)&sv,size,sizeof(SPRING_VERTEX));
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,total_size,sizeof(GLOBAL_POINT));
	    }

        for (i = 0; i < max_point_gindex; ++i)
            point_set[i] = nullptr;
        
        //allocate mem for point_set via point_set_store, and store
        //gindex values of elastic_intfc points in array positions
	    link_point_set(&geom_set,point_set,point_set_store);
	    
        count_vertex_neighbors(&geom_set,sv);
	    set_spring_vertex_memory(sv,client_size);
	    
        //links sv to point_set
        set_vertex_neighbors(&geom_set,sv,point_set);

        //Write from client geom_set to client point_set
        // (which sv has pointers to)
	    get_point_set_from(&geom_set,point_set);

        //Send client point_sets to owner
	    pp_send(5,L,MAXD*sizeof(double),owner_id);
	    pp_send(6,U,MAXD*sizeof(double),owner_id);
	    pp_send(1,&(total_client_size),sizeof(int),owner_id);
        pp_send(2,point_set_store,total_client_size*sizeof(GLOBAL_POINT),owner_id);
	}
	else
    {
	    size = owner_size;
	    total_size = total_owner_size;
    }


    CollisionSolver3d* collision_solver = nullptr;

    if (!debugging("collision_off") && FT_Dimension() == 3) 
    {
        collision_solver = new CollisionSolver3d();
        printf("COLLISION DETECTION ON\n");

        setCollisionFreePoints3d(elastic_intfc);
        
        //collision_solver->initializeSystem(*newfront); //TODO: may use this now
        CollisionSolver3d::setStep((*newfront)->step);
        CollisionSolver3d::setTimeStepSize((*newfront)->dt);
        CollisionSolver3d::setOutputDirectory(OutName(*newfront));    

        collision_solver->assembleFromInterface(elastic_intfc);
        collision_solver->recordOriginalPosition();
        collision_solver->setHseTypeLists();
        collision_solver->initializeImpactZones();

        collision_solver->setRestitutionCoef(1.0);
        collision_solver->setVolumeDiff(0.0);

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
        collision_solver->setCompressiveStrainLimit(af_params->compressive_strain_limit);
        collision_solver->setStrainRateLimit(af_params->strainrate_limit);

        collision_solver->gpoints = (*newfront)->gpoints;
        collision_solver->gtris = (*newfront)->gtris;
    }
    else
    {
        printf("COLLISION DETECTION OFF\n");
    }


	if (myid == owner_id)
	{
        //Write from owner geom_set to owner point_set
        get_point_set_from(&geom_set,point_set);

        //Write from client point_sets to owner point_set
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
                        client_size_new[i],sizeof(GLOBAL_POINT));
            }

            pp_recv(2,i,client_point_set_store[i],
                    client_size_new[i]*sizeof(GLOBAL_POINT));

            //performs the actual write
            copy_from_client_point_set(point_set,client_point_set_store[i],
                    client_size_new[i],client_L,client_U);
	    } 

        //Call spring solver now that owner holds all the point set data
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

        //Owner writes back to client point sets and sends to clients
	    for (i = 0; i < pp_numnodes(); i++)
        {
            if (i == myid) continue;
            copy_to_client_point_set(point_set,
                    client_point_set_store[i],client_size_new[i]);
            pp_send(3,client_point_set_store[i],
                    client_size_new[i]*sizeof(GLOBAL_POINT),i);
        }
	
    }//end myid == owner_id


    //Clients receive their point sets back from owner
    if (myid != owner_id)
    {
        pp_recv(3,owner_id,point_set_store,total_client_size*sizeof(GLOBAL_POINT));
    }

	//All processes write from point sets to their geom_sets
	put_point_set_to(&geom_set,point_set);

	// Calculate the real force on load_node and rg_string_node
    setSpecialNodeForce(elastic_intfc,geom_set.kl);
        //setSpecialNodeForce(fr,geom_set.kl);

	set_vertex_impulse(&geom_set,point_set);
	set_geomset_velocity(&geom_set,point_set);
	compute_center_of_mass_velo(&geom_set);

    if (!debugging("collision_off") && FT_Dimension() == 3) 
    {
        start_clock("resolveCollision");
        collision_solver->resolveCollisionSubstep();
        stop_clock("resolveCollision");
        delete collision_solver;

        setSpecialNodeForce(elastic_intfc,geom_set.kl);
        compute_center_of_mass_velo(&geom_set);
    }


    //sync interfaces after collision handling
    pp_gsync();
	exchange_curve_gindex(*newfront);
	exchange_surf_gindex(*newfront);
    scatter_front(*newfront);

    if (debugging("max_speed"))
    {
        print_max_fabric_speed(*newfront);
        print_max_string_speed(*newfront);
    }

	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate_parallel()\n");
}	/* end new_fourth_order_elastic_set_propagate_parallel_1() */


/*
void fourth_order_elastic_set_propagate_parallel(Front* fr, double fr_dt)
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
	double client_L[MAXD],client_U[MAXD];
	static boolean first_break_strings = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;


	if (debugging("trace"))
	    (void) printf("Entering fourth_order_elastic_set_propagate()\n");
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
            //TODO: Including NEUMANN_BOUNDARY in the w_type array can
            //      cause problems with rectangular domain numann boundarys.
            //      Need to put a continue statement inside a domain boundary
            //      check of the hypersurface in collect_hyper_surfaces()
            int w_type[3] = {ELASTIC_BOUNDARY,MOVABLE_BODY_BOUNDARY,NEUMANN_BOUNDARY};
	        elastic_intfc = collect_hyper_surfaces(fr,owner,w_type,3);
            //TODO: Or use this one???
                //elastic_intfc = FT_CollectHypersurfFromSubdomains(fr,owner,ELASTIC_BOUNDARY);
            
            collectNodeExtra(fr,elastic_intfc,owner_id);
	    }
	    else
        {
            elastic_intfc = fr->interf;
        }

        start_clock("set_data");
	    if (myid == owner_id)
        {
            if (client_size_old != NULL)
            {
                FT_FreeThese(3, client_size_old, client_size_new, 
                        client_point_set_store);
            }

            FT_VectorMemoryAlloc((POINTER*)&client_size_old,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_size_new,
                    pp_numnodes(),sizeof(int));
            FT_VectorMemoryAlloc((POINTER*)&client_point_set_store,
                    pp_numnodes(),sizeof(GLOBAL_POINT*));
            
            for (i = 0; i < pp_numnodes(); i++)
            {
                client_size_old[i] = 0;
                client_size_new[i] = 0;
            }

            assembleParachuteSet(elastic_intfc,&geom_set);
            owner_size = geom_set.num_verts;
            
            if (point_set_store != NULL) 
                FT_FreeThese(2,point_set_store, sv);
            
            FT_VectorMemoryAlloc((POINTER*)&point_set_store,
                    owner_size,sizeof(GLOBAL_POINT));
            FT_VectorMemoryAlloc((POINTER*)&sv,owner_size,
                    sizeof(SPRING_VERTEX));
            
            //allocate mem for point_set via point_set_store, and store
            //elastic_intfc points in the point_set array mapped to their
            //respective gindex values.
            link_point_set(&geom_set,point_set,point_set_store);

            count_vertex_neighbors(&geom_set,sv);
            set_spring_vertex_memory(sv,owner_size);

            //links sv to the point_set
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
        
        //allocate mem for point_set via point_set_store, and store
        //elastic_intfc points in the point_set array mapped to their
        //respective gindex values.
	    link_point_set(&geom_set,point_set,point_set_store);
	    
        count_vertex_neighbors(&geom_set,sv);
	    set_spring_vertex_memory(sv,client_size);
	    
        //links sv to the point_set
        set_vertex_neighbors(&geom_set,sv,point_set);

        //Write from client's geom_set to client's point_set
        //(which sv has pointers to)
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


    CollisionSolver3d* collision_solver = nullptr;


    //Collect all point_set data at owner node then
    //run fabric solver and collision solver. Then return to clients
	if (myid == owner_id)
	{
        if (!debugging("collision_off") && FT_Dimension() == 3) 
        {
            printf("COLLISION DETECTION ON\n");
            collision_solver = new CollisionSolver3d();

            //NOTE: These two functions need to be called before the spring
            //      solver runs.
            setCollisionFreePoints3d(fr->interf);
            collision_solver->initializeSystem(fr);

            //These two also must be set before??
            collision_solver->gpoints = fr->gpoints;
            collision_solver->gtris = fr->gtris;

            //TODO: consolidate the following into a setDefaultCollisionParams() function
            collision_solver->setRestitutionCoef(1.0);
            collision_solver->setVolumeDiff(0.0);

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

            //collision_solver->setStrainLimit(af_params->strain_limit);
            //collision_solver->setStrainRateLimit(af_params->strainrate_limit);
        }
        else
        {
            printf("COLLISION DETECTION OFF\n");
        }


        //Write from owner geom_set to owner point_set
	    get_point_set_from(&geom_set,point_set);

        //Write from client point_sets to owner point_set
	    for (i = 0; i < pp_numnodes(); i++)
	    {
            //if (i == owner_id) continue;
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

            //receive point_set_stores from clients
            pp_recv(2,i,client_point_set_store[i],
                    client_size_new[i]*sizeof(GLOBAL_POINT));

            //performs the actual write to owner point_set
            copy_from_client_point_set(point_set,client_point_set_store[i],
                    client_size_new[i],client_L,client_U);
	    } 

        //Call spring solver now that owner holds all the point set data
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

        //Write back to owner geomset from owner point set
        //(which was updated through sv in the spring solver)
        put_point_set_to(&geom_set,point_set);
        set_vertex_impulse(&geom_set,point_set);
	    set_geomset_velocity(&geom_set,point_set);


        if (!debugging("collision_off") && FT_Dimension() == 3) 
        {
            start_clock("resolveCollision");
            collision_solver->resolveCollisionSubstep();
            stop_clock("resolveCollision");
            delete collision_solver;
        }

        //Write from owner geom_set to owner point_set
	    get_point_set_from(&geom_set,point_set);


        //Owner writes back to client point sets and sends
	    for (i = 0; i < pp_numnodes(); i++)
        {
            if (i == myid) continue;
            copy_to_client_point_set(point_set,
                    client_point_set_store[i], client_size_new[i]);
            pp_send(3,client_point_set_store[i],
                            client_size_new[i]*sizeof(GLOBAL_POINT),i);
        }
	
    }//end myid == owner_id


    //Clients recevie their point sets back from owner
    if (myid != owner_id)
    {
        pp_recv(3,owner_id,point_set_store,
            client_size*sizeof(GLOBAL_POINT));
    }

	//All processes write from point sets to their geom_sets
	put_point_set_to(&geom_set,point_set);

	// Calculate the real force on load_node and rg_string_node
	setSpecialNodeForce(fr,geom_set.kl);

	set_vertex_impulse(&geom_set,point_set);
	set_geomset_velocity(&geom_set,point_set);
	compute_center_of_mass_velo(&geom_set);

    if (debugging("max_speed"))
    {
        print_max_fabric_speed(fr);
        print_max_string_speed(fr);
    }

	if (debugging("trace"))
	    (void) printf("Leaving fourth_order_elastic_set_propagate_parallel()\n");
}*/	/* end fourth_order_elastic_set_propagate_parallel() */

static void print_max_fabric_speed(Front* fr)
{
    SURFACE **s;
    TRI *tri;
    POINT *pt;
    STATE *state;
    
    double max_speed = -HUGE;
    POINT* max_pt = nullptr;

    if (!FT_FrontContainWaveType(fr,ELASTIC_BOUNDARY)) return;

    intfc_surface_loop(fr->interf,s)
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
        surf_tri_loop(*s,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                pt = Point_of_tri(tri)[i];
	            
                double nor[MAXD];
                FT_NormalAtPoint(pt,fr,nor,NO_COMP);
                
                state = (STATE*)left_state(pt);
                double* vel = state->vel;
                double nor_speed = scalar_product(vel,nor,3);
                double abs_nor_speed = std::abs(nor_speed);
                    //double speed = sqrt(sqr(state->vel[0]) + sqr(state->vel[1]) + sqr(state->vel[2]));
                
                if (max_speed < abs_nor_speed)
                {
                    max_speed = abs_nor_speed;
                    max_pt = pt;
                }
            }
        }
    }
    
    if (max_pt != nullptr)
    {
        printf("\n");
        state = (STATE*)left_state(max_pt);
        printf("max normal speed of fabric/canopy: %g\n",max_speed);
        printf("Velocity: %g %g %g\n",
                state->vel[0],state->vel[1],state->vel[2]);
        printf("coords = %f %f %f    Point Gindex: %d\n",
                Coords(max_pt)[0],Coords(max_pt)[1],Coords(max_pt)[2],
                Gindex(max_pt));
    }

    if (max_speed > 1000)
    {
        printf("\n\n\tWARNING: max normal speed of fabric/canopy exceeds 1000 m/s\n\n");
            //LOC(); clean_up(EXIT_FAILURE);
    }
}

static void print_max_string_speed(Front* fr)
{
    CURVE **c;
    CURVE *curve;
    BOND *b;
    POINT *pt;
    STATE *state;
    
    double speed;
    double max_speed = -HUGE;
    POINT* max_pt = nullptr;

    if (!FT_FrontContainHsbdryType(fr,STRING_HSBDRY)) return;

    intfc_curve_loop(fr->interf,c)
    {
        if (hsbdry_type(*c) != STRING_HSBDRY) continue;

        //TODO: Add nodes etc. below is not traversing every point of the curve.
        curve = *c;
        for (b = curve->first; b != curve->last; b = b->next)
        {
            pt = b->end;
            state = (STATE*)left_state(pt);
            speed = sqrt(sqr(state->vel[0]) + sqr(state->vel[1])
                        + sqr(state->vel[2]));
            if (max_speed < speed)
            {
                max_speed = speed;
                max_pt = pt;
            }
        }
    }
    
    if (max_pt != nullptr)
    {
        printf("\n");
        state = (STATE*)left_state(max_pt);
        printf("max speed of elastic strings: %g\n",max_speed);
        printf("Velocity: %g %g %g\n",
                state->vel[0],state->vel[1],state->vel[2]);
        printf("coords = %f %f %f    Point Gindex: %d\n",
                Coords(max_pt)[0],Coords(max_pt)[1],Coords(max_pt)[2],
                Gindex(max_pt));
    }

    if (max_speed > 1000)
    {
        printf("\n\n\tWARNING: max speed of elastic strings exceeds 1000 m/s\n\n");
            //LOC(); clean_up(EXIT_FAILURE);
    }
}

static void setSurfVelocity(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	GLOBAL_POINT **point_set)
{
	int i,j;
	TRI *tri;
	POINT *p;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
	Front *front = geom_set->front;
	double nor[MAXD];
	double *vel = nullptr;
	int gindex_max;
	long gindex;

    int dim = front->rect_grid->dim;
    double nor_speed = 0.0;
    double max_nor_speed = 0.0;
    double *max_coords = nullptr;

    //temporary variable to test reduce_high_freq_vel() function 
    static double prev_max_nor_speed = HUGE;

    /*
    static STATE* max_nor_speed_state;
    if (max_nor_speed_state == nullptr)
    {
        FT_ScalarMemoryAlloc((POINTER*)&max_nor_speed_state,sizeof(STATE));
    }
    */

	unsort_surf_point(surf);
	HYPER_SURF* hs = Hyper_surf(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
        {
            p = Point_of_tri(tri)[i];
            if (sorted(p) || Boundary_point(p)) continue;
            
            gindex = Gindex(p);
            FT_NormalAtPoint(p,front,nor,NO_COMP);
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            vel = point_set[gindex]->v;
            nor_speed = scalar_product(vel,nor,3);
            
            for (j = 0; j < 3; ++j)
            {
                sl->vel[j] = vel[j];
                sr->vel[j] = vel[j];
                    //sl->vel[j] = nor_speed*nor[j];
                    //sr->vel[j] = nor_speed*nor[j];
            }
            
            if (max_nor_speed < fabs(nor_speed))
            {
                max_nor_speed = fabs(nor_speed);
                max_coords = Coords(p);
                    //ft_assign(max_nor_speed_state,sl,fr->sizest);
            }

            sorted(p) = YES;
        }
	}

    set_max_front_speed(dim,max_nor_speed,NULL,max_coords,front);
    
    //TODO: add a switch for using reduce_high_freq_vel() when max_speed spikes abruptly
    if (debugging("reduce_highfreq_vel"))
    {
        //Hardcode 3x prev_max_nor_speed to test
        if (max_nor_speed > 3.0*prev_max_nor_speed)
        {
            printf("\nEntering reduce_high_freq_vel()\n");
            reduce_high_freq_vel(front,surf);
            prev_max_nor_speed = Spfr(front)[3];
        }
        else
        {
            prev_max_nor_speed = max_nor_speed;
        }
    }
}	/* end setSurfVelocity */

static void setCurveVelocity(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	GLOBAL_POINT **point_set)
{
	int i,j;
	BOND *b;
	POINT *p;
	BOND_TRI **btris;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	Front *front = geom_set->front;
	double nor[MAXD];
	double *vel = nullptr;
	double *crds_max = nullptr;
	int gindex_max;
	long gindex;
	int dim = FT_Dimension();

    if (hsbdry_type(curve) == STRING_HSBDRY)
    {
        double max_vel[3] = {0};
        double max_speed = 0.0;

        /*
        //TODO: Nodes set in other functions, and in this case if max speed is set
        //      then it would be overwriting the value set when traversing the mono_comp_hsbry
        //      reducing the time step unnecessarily.
        
        NODE* string_nodes[2];
        string_nodes[0] = curve->start;
        string_nodes[1] = curve->end;

        for (int i = 0; i < 2; ++i)
        {
            if (!is_string_node(string_nodes[i])) continue;

            p = curve->start->posn;
            gindex = Gindex(p);
            sl = (STATE*)left_state(p);
            sr = (STATE*)right_state(p);
            vel = point_set[gindex]->v;
            
            double speed = Mag3d(vel);
            if (max_speed < speed)
            {
                max_speed = speed;
                crds_max = Coords(p);
            }

            for (j = 0; j < 3; ++j)
            {
                sl->vel[j] = vel[j];
                sr->vel[j] = vel[j];
            }
        }
        */

        for (b = curve->first; b != curve->last; b = b->next)
        {
            p = b->end;
            gindex = Gindex(p);
            sl = (STATE*)left_state(p);
            sr = (STATE*)right_state(p);
            vel = point_set[gindex]->v;
            
            double speed = Mag3d(vel);
            if (max_speed < speed)
            {
                max_speed = speed;
                crds_max = Coords(p);
            }

            for (j = 0; j < 3; ++j)
            {
                sl->vel[j] = vel[j];
                sr->vel[j] = vel[j];

                if (std::abs(max_vel[j]) < std::abs(vel[j]))
                {
                    max_vel[j] = vel[j];
                }
            }
        }

        for (int k = 0; k < dim; ++k)
        {
            set_max_front_speed(k,max_vel[k],nullptr,nullptr,front);
        }
    }
    else
    {
        double max_nor_speed = 0.0;

        for (b = curve->first; b != curve->last; b = b->next)
        {
            p = b->end;
            for (btris = Btris(b); btris && *btris; ++btris)
            {
                p->hse = hse = Hyper_surf_element((*btris)->tri);
                p->hs = hs = Hyper_surf((*btris)->surface);
                gindex = Gindex(p);
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                FT_NormalAtPoint(p,front,nor,NO_COMP);
                vel = point_set[gindex]->v;

                double nor_speed = scalar_product(vel,nor,3);
                if (max_nor_speed < fabs(nor_speed))
                {
                    max_nor_speed = fabs(nor_speed);
                    crds_max = Coords(p);
                }
                
                for (j = 0; j < 3; ++j)
                {
                    sl->vel[j] = vel[j];
                    sr->vel[j] = vel[j];
                        //sl->vel[j] = nor_speed*nor[j];
                        //sr->vel[j] = nor_speed*nor[j];
                }
            }
        }
    
        set_max_front_speed(dim,max_nor_speed,NULL,crds_max,front);
    }

    /*
    //TODO: How significant is this????
    for (b = curve->first; b != NULL; b = b->next)
    {
	    set_bond_length(b,dim);
    }
    */
}	/* end setCurveVelocity */

static void setNodeVelocity(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set)
{
	int dim = FT_Dimension();
	switch(dim)
	{
	case 2:
	    setNodeVelocity2d(geom_set,node,point_set);
	    return;
	case 3:
	    setNodeVelocity3d(geom_set,node,point_set);
	    return;
	}
}	/* end setNodeVelocity */

static void setNodeVelocity2d(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set)
{
	int i,j;
	BOND *b;
	POINT *p;
	STATE *sl,*sr;
    Front *front = geom_set->front;
	double *vel = nullptr;
	long gindex;

    double max_speed = 0.0;

	if (is_load_node(node))
	{
	    sl = (STATE*)left_state(node->posn);
        sr = (STATE*)right_state(node->posn);
	    gindex = Gindex(node->posn);
	    vel = point_set[gindex]->v;
        
        for (j = 0; j < 2; ++j)
        {
            sl->vel[j] = vel[j];
            sr->vel[j] = vel[j];
            max_speed += sqr(vel[j]);
        }

        max_speed = sqrt(max_speed);
        set_max_front_speed(2,max_speed,NULL,Coords(node->posn),front);
	}
}	/* end setNodeVelocity2d */

static void setNodeVelocity3d(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set)
{
	int i,j;
	BOND *b;
	POINT *p;
	BOND_TRI **btris;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs;
	Front *front = geom_set->front;
	CURVE **c;
	double nor[MAXD];
	double *vel = nullptr;
	int gindex_max;
	long gindex;

    double nor_speed = 0.0;
    double max_nor_speed = 0.0;
    double *crds_max = nullptr;

	for (c = node->out_curves; c && *c; ++c)
    {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY && 
		    hsbdry_type(*c) != PASSIVE_HSBDRY) continue;

        b = (*c)->first;
        p = b->start;
		for (btris = Btris(b); btris && *btris; ++btris)
		{
            p->hse = hse = Hyper_surf_element((*btris)->tri);
            p->hs = hs = Hyper_surf((*btris)->surface);
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
		    FT_NormalAtPoint(p,front,nor,NO_COMP);
		    gindex = Gindex(p);
		    vel = point_set[gindex]->v;
		    
            if (hsbdry_type(*c) == PASSIVE_HSBDRY)
		    {
                //This sets the rg_string_node vel
                for (j = 0; j < 3; ++j)
                {
                    sl->vel[j] = vel[j];
                    sr->vel[j] = vel[j];
                }
    			continue;
		    }

		    nor_speed = scalar_product(vel,nor,3);
		    if (max_nor_speed < fabs(nor_speed)) 
		    {
		    	max_nor_speed = fabs(nor_speed);
		    	gindex_max = Gindex(p);
                crds_max = Coords(p);
		    }

            for (j = 0; j < 3; ++j)
            {
                sl->vel[j] = vel[j];
                sr->vel[j] = vel[j];
                    //sl->vel[j] = nor_speed*nor[j];
                    //sr->vel[j] = nor_speed*nor[j];
            }
		}
        
    }

    for (c = node->in_curves; c && *c; ++c)
    {
        if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
            hsbdry_type(*c) != GORE_HSBDRY && 
            hsbdry_type(*c) != PASSIVE_HSBDRY) continue;

        b = (*c)->last;
        p = b->end;

        for (btris = Btris(b); btris && *btris; ++btris)
        {
            p->hse = hse = Hyper_surf_element((*btris)->tri);
            p->hs = hs = Hyper_surf((*btris)->surface);
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            FT_NormalAtPoint(p,front,nor,NO_COMP);
            gindex = Gindex(p);
            vel = point_set[gindex]->v;
        
            if (hsbdry_type(*c) == PASSIVE_HSBDRY)
            {
                for (j = 0; j < 3; ++j)
                {
                    sl->vel[j] = vel[j];
                    sr->vel[j] = vel[j];
                }
                continue;
            }

            nor_speed = scalar_product(vel,nor,3);
            if (max_nor_speed < fabs(nor_speed))
            {
                max_nor_speed = fabs(nor_speed);
                gindex_max = Gindex(p);
                crds_max = Coords(p);
            }

            for (j = 0; j < 3; ++j)
            {
                sl->vel[j] = vel[j];
                sr->vel[j] = vel[j];
                    //sl->vel[j] = nor_speed*nor[j];
                    //sr->vel[j] = nor_speed*nor[j];
            }
        }
    }

    set_max_front_speed(3,max_nor_speed,NULL,crds_max,front);
}	/* end setNodeVelocity3d */

extern void set_geomset_velocity(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set)
{
	int ns = geom_set->num_surfs;
	int nc = geom_set->num_curves;
	int nn = geom_set->num_nodes;

	for (int i = 0; i < ns; ++i)
    {
	    setSurfVelocity(geom_set,geom_set->surfs[i],point_set);
    }
    for (int i = 0; i < nc; ++i)
    {
        setCurveVelocity(geom_set,geom_set->curves[i],point_set);
    }
    for (int i = 0; i < nn; ++i)
	{
        if (is_load_node(geom_set->nodes[i]))
        {
            continue;
        }
        setNodeVelocity(geom_set,geom_set->nodes[i],point_set);
	}

    //TODO: add rgb_surfs ?

}	/* end set_geomset_velocity */

extern void collectNodeExtra(
	Front *front,
	INTERFACE *host_intfc,
	int owner_id)
{
	NODE **n;
	NODE **node_with_extra;
	int i,j,k,num_nodes;
	INTERFACE *intfc = front->interf;
	long global_index;
	AF_NODE_EXTRA extra_recv, *extra;
	RECT_GRID *gr = front->rect_grid;
        double *L = gr->L;
        double *U = gr->U;
        int dim = gr->dim;

	if (pp_mynode() != owner_id)
	{
	    num_nodes = 0;
	    intfc_node_loop(intfc,n)
	    {
		for (k = 0; k < dim; ++k)
		{
		    if (Coords((*n)->posn)[k] <= L[k] ||
		        Coords((*n)->posn)[k] > U[k])
		        break;
		}
		if (k != dim || (*n)->extra == NULL) continue;
		num_nodes++;
	    }
	    pp_send(10,&num_nodes,sizeof(int),owner_id);
	    intfc_node_loop(intfc,n)
	    {
		for (k = 0; k < dim; ++k)
		{
		    if (Coords((*n)->posn)[k] <= L[k] ||
		        Coords((*n)->posn)[k] > U[k])
		        break;
		}
		if (k != dim || (*n)->extra == NULL) continue;
	    	pp_send(11,&(Gindex((*n)->posn)),sizeof(long),owner_id);
	    	pp_send(12,(*n)->extra,sizeof(AF_NODE_EXTRA),owner_id);
	    }
	}
	else
	{
	    for (i = 0; i < pp_numnodes(); ++i)
	    {
		if (i == owner_id) continue;
		pp_recv(10,i,&num_nodes,sizeof(int));
		for (j = 0; j < num_nodes; ++j)
		{
		    pp_recv(11,i,&global_index,sizeof(long));
		    pp_recv(12,i,&extra_recv,sizeof(AF_NODE_EXTRA));
		    intfc_node_loop(host_intfc,n)
		    {
			if (Gindex((*n)->posn) != global_index)
			    continue;
		        FT_ScalarMemoryAlloc((POINTER*)&extra,
					sizeof(AF_NODE_EXTRA));
			*extra = extra_recv;
			(*n)->extra = (POINTER)extra;
			(*n)->size_of_extra = sizeof(AF_NODE_EXTRA);
		    }
		}
	    }
	}
}	/* end collectNodeExtra */

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

    //TODO: is this loop more efficient than using the intfc_surface_loop()?
    next_point(intfc,NULL,NULL,NULL);
    while(next_point(intfc,&p,&hse,&hs))
    {
        STATE* sl = (STATE*)left_state(p);
        sl->is_fixed = false;
        sl->is_movableRG = false;
        sl->is_registeredpt = false;
        
        if ((surf = Surface_of_hs(hs)))
        {
            if (wave_type(hs) == NEUMANN_BOUNDARY)
                sl->is_fixed = true;
            if (wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                sl->is_movableRG = true;
            if (is_registered_point(surf,p))
            {
                sl->is_registeredpt = true;
                sl->is_fixed = true;
            }
        }
    }

    CURVE **c;
    BOND* b;
    intfc_curve_loop(intfc,c)
    {
        if (hsbdry_type(*c) != FIXED_HSBDRY &&
            hsbdry_type(*c) != STRING_HSBDRY) continue;

        for (b = (*c)->first; b != (*c)->last; b = b->next)
        {
            STATE* sl = (STATE*)left_state(b->end);
            if (hsbdry_type(*c) == FIXED_HSBDRY)
                sl->is_fixed = true;
            /*
            else if (is_registered_point(*c,b->end))
                sl->is_registeredpt = true;
            */ //TODO: define registered points on curves!
        }
    }

    NODE** n;
    intfc_node_loop(intfc,n)
    {
        STATE* sl = (STATE*)left_state((*n)->posn);
        sl->is_fixed = false;
        sl->is_movableRG = false;
        sl->is_load_node = false;

        AF_NODE_EXTRA* extra = (AF_NODE_EXTRA*)(*n)->extra;
        if (extra)
        {
            if (extra->af_node_type == RG_STRING_NODE)
            {
                node_out_curve_loop(*n,c)
                {
                    if (hsbdry_type(*c) == PASSIVE_HSBDRY)
                    {
                        b = (*c)->first;
                        hs = Hyper_surf(b->_btris[0]->surface);
                        break;
                    }
                }

                node_in_curve_loop(*n,c)
                {
                    if (hsbdry_type(*c) == PASSIVE_HSBDRY)
                    {
                        b = (*c)->last;
                        hs = Hyper_surf(b->_btris[0]->surface);
                        break;
                    }
                }

                if (hs == NULL)
                {
                    printf("ERROR in setCollisionFreePoints3d() : ");
                    printf("No related hs found\n");
                    LOC(); clean_up(ERROR);
                }

                if (wave_type(hs) == NEUMANN_BOUNDARY)
                    sl->is_fixed = true;
                if (wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                    sl->is_movableRG = true;
            }
            else if (extra->af_node_type == LOAD_NODE)
            {
                sl->is_load_node = true;
            }
            else if (extra->af_node_type == PRESET_NODE)
            {
                sl->is_fixed = true;
            }
        }
        else if ((*n)->hsb && is_fixed_node(*n))
        {
            sl->is_fixed = true;
        }
    }
}       /* setCollisionFreePoints3d() */

extern void scatterAirfoilExtra(
	Front *front)
{
	NODE **n;
	NODE **node_with_extra;
	int i,j,k,num_nodes;
	INTERFACE *intfc = front->interf;
	long global_index;
	AF_NODE_EXTRA extra_recv,*extra;
	RECT_GRID *gr = front->rect_grid;
	double *L = gr->L;
	double *U = gr->U;
	int dim = gr->dim;

	num_nodes = 0;
        intfc_node_loop(intfc,n)
        {
            for (k = 0; k < dim; ++k)
            {
                if (Coords((*n)->posn)[k] <= L[k] ||
                    Coords((*n)->posn)[k] > U[k])
                    break;
            }
            if (k != dim || (*n)->extra == NULL) continue;
            num_nodes++;
        }
        for (i = 0; i < pp_numnodes(); ++i)
        {
            if (i == pp_mynode()) continue;
            pp_send(30,&num_nodes,sizeof(int),i);
            intfc_node_loop(intfc,n)
            {
                for (k = 0; k < dim; ++k)
                {
                    if (Coords((*n)->posn)[k] <= L[k] ||
                        Coords((*n)->posn)[k] > U[k])
                        break;
                }
                if (k != dim || (*n)->extra == NULL) continue;
                pp_send(31,&(Gindex((*n)->posn)),sizeof(long),i);
                pp_send(32,(*n)->extra,sizeof(AF_NODE_EXTRA),i);
            }
        }
	pp_gsync();
        for (i = 0; i < pp_numnodes(); ++i)
        {
            if (i == pp_mynode()) continue;
            pp_recv(30,i,&num_nodes,sizeof(int));
            for (j = 0; j < num_nodes; ++j)
            {
                pp_recv(31,i,&global_index,sizeof(long));
                pp_recv(32,i,&extra_recv,sizeof(AF_NODE_EXTRA));
                intfc_node_loop(intfc,n)
                {
                    if (Gindex((*n)->posn) != global_index)
                        continue;
                    FT_ScalarMemoryAlloc((POINTER*)&extra,
                                sizeof(AF_NODE_EXTRA));
                    *extra = extra_recv;
                    (*n)->extra = (POINTER)extra;
                    (*n)->size_of_extra = sizeof(AF_NODE_EXTRA);
                }
            }
	}
}	/* end scatterAirfoilExtra */

extern void setSpecialNodeForce(
	INTERFACE* intfc,
	double kl)
{
	int i, k;
	double f[MAXD], vec[MAXD];
	NODE **n;
	CURVE **c;
	BOND *b;
	RECT_GRID *gr = &(intfc->table->rect_grid);
	double *L = gr->L;
	double *U = gr->U;
	int dim = gr->dim;
	
	if (debugging("trace"))
	    printf("Entering setSpecialNodeForce() \n");

	intfc_node_loop(intfc, n)
	{
	    if ((!is_load_node(*n)) && (!is_rg_string_node(*n))) continue;

        for (k = 0; k < dim; ++k)
        {
            if (Coords((*n)->posn)[k] <= L[k] ||
                Coords((*n)->posn)[k] > U[k]) break;
        }
        if (k != dim) continue;

	    for (i = 0; i < dim; ++i) f[i] = 0.0;
	    
        node_out_curve_loop(*n,c)
	    {
            b = (*c)->first;
            set_bond_length(b,dim);
	    }
	    
        node_in_curve_loop(*n,c)
	    {
            b = (*c)->last;
            set_bond_length(b,dim);
	    }
	    
        node_out_curve_loop(*n,c)
	    {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
            
            b = (*c)->first;
            double dL = bond_length(b) - bond_length0(b);
            
            for (i = 0; i < dim; ++i)
            {
                vec[i] = Coords(b->end)[i] - Coords(b->start)[i];
                vec[i] /= bond_length(b);
                f[i] += kl*dL*vec[i];
            }
	    }
	    
        node_in_curve_loop(*n,c)
	    {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
            
            b = (*c)->last;
            double dL = bond_length(b) - bond_length0(b);
            
            for (i = 0; i < dim; ++i)
            {
                vec[i] = Coords(b->start)[i] - Coords(b->end)[i];
                vec[i] /= bond_length(b);
                f[i] += kl*dL*vec[i];
            }
	    }

        //TODO: What about the weight of the mass?
	    
        //The update
        for (i = 0; i < dim; ++i)
        {
            (*n)->posn->force[i] = f[i];
        }

	    if (debugging("rigid_body"))
	    {
            printf("special node coords = %f %f %f \n", 
                Coords((*n)->posn)[0], Coords((*n)->posn)[1], 
                Coords((*n)->posn)[2]);
            printf("force on the node = %f %f %f \n", f[0], f[1], f[2]);
            printf("velo of the node = %f %f %f \n", (*n)->posn->vel[0], 
                    (*n)->posn->vel[1], (*n)->posn->vel[2]);
	    }
	}

	if (debugging("trace"))
	    printf("Leaving setSpecialNodeForce() \n");
}	/* end setSpecialNodeForce */

extern void break_strings(Front *front)
{
        INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
        static boolean first = YES;
	static double break_strings_time = af_params->break_strings_time;
	static int break_strings_num = af_params->break_strings_num;
	static int *bs_gindex = af_params->break_strings_gindex;
        CURVE **c, **curves;
        int i = 0;

	if (FT_Dimension() != 3) return;
	if (break_strings_time < 0.0 || break_strings_num <= 0) return;
        if ((!first) || (front->time + front->dt < break_strings_time)) return;

	if (debugging("trace"))
	    printf("Entering break_strings() \n");

        first = NO;
	double z_max = -HUGE, z_min = HUGE;
	intfc_curve_loop(intfc, c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	    if (z_max < Coords((*c)->start->posn)[2])
		z_max = Coords((*c)->start->posn)[2];
	    if (z_max < Coords((*c)->end->posn)[2])
		z_max = Coords((*c)->end->posn)[2];
	    if (z_min > Coords((*c)->start->posn)[2])
		z_min = Coords((*c)->start->posn)[2];
	    if (z_min > Coords((*c)->end->posn)[2])
		z_min = Coords((*c)->end->posn)[2];
	}
	pp_gsync();
	pp_global_max(&z_max, 1);
	pp_global_min(&z_min, 1);
	double break_height = 0.5 * (z_max + z_min);
	for (int i = 0; i < break_strings_num; ++i)
	{
	    intfc_curve_loop(intfc, c)
	    {
		if (Gindex(*c) == bs_gindex[i]) 
		    break;
	    }
	    if (c && *c)
		break_string_curve(*c, break_height);
	    pp_gsync();
	    pp_global_lmax(&(current_interface()->max_point_gindex),1);
	    pp_global_lmax(&(current_interface()->max_curve_gindex),1);
	    exchange_curve_gindex(front);
	    scatter_front(front);
	}

	if (debugging("trace"))
	    printf("Leaving break_strings() \n");
}	/* end break_strings */

static void break_string_curve(CURVE* c, double height)
{
	static CURVE *curves[2] = {NULL, NULL};

	if (debugging("trace"))
	    printf("Entering split_string_curve()\n");
	INTERFACE *cur_intfc = current_interface();
	RECT_GRID *gr = computational_grid(cur_intfc);
	BOND *b = NULL;

	set_current_interface(c->interface);
	curve_bond_loop(c, b)
	{
	    double tmp = (Coords(b->start)[2] - height) * 
				(Coords(b->end)[2] - height);
	    if (tmp <= 0.0) break;
	}
	if (b == c->first) b = b->next;
	if (b == NULL || c->num_points <= 2)
	{
	    set_current_interface(cur_intfc);
	    return;
	}

	POINT *p[2];
	p[0] = b->start;
	p[1] = Point(Coords(p[0]));
	/* for consistency and uniqueness of global index */
	Gindex(p[1]) = cur_intfc->max_point_gindex + 1;
	cur_intfc->max_point_gindex += 1;

	NODE *n[2];
	n[0] = make_node(p[0]);
	n[1] = make_node(p[1]);
	if (n[0] == NULL || n[1] == NULL)
	{
	    printf("ERROR: split_string_curve returning NULL, "
	    	"cannot split at point p\n");
	    clean_up(ERROR);
	}
	curves[0] = make_curve(0,0,c->start,n[0]);
	curves[1] = make_curve(0,0,n[1],c->end);
	if (curves[0] == NULL || curves[1] == NULL)
	{
	    printf("ERROR: split_string_curve returning NULL, "
	    	"cannot make curve\n");
	    clean_up(ERROR);
	}
	curves[0]->first = c->first;
	curves[0]->last = b->prev;
	curves[0]->first->prev = NULL;
	curves[0]->last->next = NULL;
	curves[1]->first = b;
	curves[1]->last = c->last;
	curves[1]->first->prev = NULL;
	curves[1]->last->next = NULL;

	/* for consistency and uniqueness of global index */
	Gindex(curves[1]) = Gindex(c);
	Gindex(curves[0]) = cur_intfc->max_curve_gindex + 1;
	cur_intfc->max_curve_gindex += 1;

	hsbdry_type(curves[0]) = hsbdry_type(c);
	hsbdry_type(curves[1]) = hsbdry_type(c);
	curve_tangent_function(curves[0]) = curve_tangent_function(c);
	curve_tangent_function(curves[1]) = curve_tangent_function(c);

	delete_curve(c);

	curves[0]->num_points = num_points_on_curve(curves[0]);
	curves[1]->num_points = num_points_on_curve(curves[1]);

	AF_NODE_EXTRA *extra[2];
	FT_ScalarMemoryAlloc((POINTER*)&extra[0],sizeof(AF_NODE_EXTRA));
	extra[0]->af_node_type = STRING_NODE;
	n[0]->extra = (POINTER)extra[0];
	n[0]->size_of_extra = sizeof(AF_NODE_EXTRA);
	FT_ScalarMemoryAlloc((POINTER*)&extra[1],sizeof(AF_NODE_EXTRA));
	extra[1]->af_node_type = STRING_NODE;
	n[1]->extra = (POINTER)extra[1];
	n[1]->size_of_extra = sizeof(AF_NODE_EXTRA);

	set_current_interface(cur_intfc);
	if (debugging("trace"))
	    printf("Leaving split_string_curve()\n");
}	/* end split_string_curve */

extern void record_break_strings_gindex(Front *front)
{
	INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int *bs_gindex = af_params->break_strings_gindex;
	CURVE **c;

	if(FT_Dimension() != 3) return;

	/* build up the map from string gindex to string id */
	/* to make the input and output more intuitive to users */
	std::vector<CURVE*> string_curves = af_params->string_curves;
	int *array, n = string_curves.size();
	FT_ScalarMemoryAlloc((POINTER*)&array, sizeof(int) * n);
	for (int i = 0; i < n; ++i)
	    array[i] = Gindex(string_curves[i]);
	pp_gsync();
	pp_global_imax(array, n);
	for (int i = 0; i < n; ++i)
	    af_params->string_hash[array[i]] = i;

	if (af_params->break_strings_num <= 0) return;

	if (debugging("trace"))
	    printf("Entering record_break_strings_gindex()\n");

	for (int i = 0; i < af_params->break_strings_num; ++i)
	    bs_gindex[i] = array[bs_gindex[i]];

	FT_FreeThese(1,array);
	if (debugging("trace"))
	    printf("Leaving record_break_strings_gindex()\n");
}	/* end record_break_strings_gindex */

extern void set_unequal_strings(Front *front)
{
	INTERFACE *intfc = front->interf;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int *us_gindex = af_params->unequal_strings_gindex;
	double tol = MACH_EPS;
	CURVE **c;

	if(FT_Dimension() != 3) return;
	if (af_params->unequal_strings_num <= 0) return;
	if (fabs(af_params->unequal_coeff - 1.0) < tol) return;

	if (debugging("trace"))
	    printf("Entering set_unequal_strings\n");

	std::vector<CURVE*> string_curves = af_params->string_curves;
	for (int i = 0; i < af_params->unequal_strings_num; ++i)
	{
	    CURVE *curve = string_curves[us_gindex[i]];
	    BOND *bond;
	    curve_bond_loop(curve, bond)
		bond->length0 *= af_params->unequal_coeff;
	}

	/* become invalid after first step and no longer used */
	af_params->string_curves.clear();

	if (debugging("trace"))
	    printf("Leaving set_unequal_strings()\n");
}	/* end set_unequal_strings */

//NOTE: nowhere used
static void linkGlobalIndexToTri(
        INTERFACE *intfc,
        TRI ***gtri)
{
        SURFACE **s;
        TRI *tri;
        int i;

        if (*gtri == NULL)
        {
            FT_VectorMemoryAlloc((POINTER*)gtri,intfc->max_tri_gindex,
                                    sizeof(TRI*));
        }
        for (i = 0; i < intfc->max_tri_gindex; ++i)
            (*gtri)[i] = NULL;
        intfc_surface_loop(intfc,s)
        {
            surf_tri_loop(*s,tri)
            {
                if (Gindex(tri) >= intfc->max_tri_gindex)
                {
                    (void) printf("ERROR: Tri global index exceeds "
                                  "max_tri_gindex!\n");
                    (void) printf("Tri gindex: %d  max_tri_gindex: %d\n",
                                Gindex(tri),intfc->max_tri_gindex);
                    clean_up(ERROR);
                }
                (*gtri)[Gindex(tri)] = tri;
            }
        }
}       /* end linkGlobalIndexToTri */

