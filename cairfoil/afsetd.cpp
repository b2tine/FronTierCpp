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

#define MAX_NUM_RING1 30

static void count_surf_neighbors(SURFACE*,SPRING_VERTEX*,int*);
static void count_curve_neighbors(CURVE*,SPRING_VERTEX*,int*);
static void count_node_neighbors(NODE*,SPRING_VERTEX*,int*);

static void link_surf_point_set(ELASTIC_SET*,SURFACE*,GLOBAL_POINT**,
				GLOBAL_POINT*,int*);
static void link_curve_point_set(ELASTIC_SET*,CURVE*,GLOBAL_POINT**,
				GLOBAL_POINT*,int*);
static void link_node_point_set(ELASTIC_SET*,NODE*,GLOBAL_POINT**,
				GLOBAL_POINT*,int*);

static void surf_get_point_set_from(SURFACE*,GLOBAL_POINT**);
static void curve_get_point_set_from(CURVE*,GLOBAL_POINT**);
static void node_get_point_set_from(NODE*,GLOBAL_POINT**);

static void surf_put_point_set_to(SURFACE*,GLOBAL_POINT**);
static void curve_put_point_set_to(CURVE*,GLOBAL_POINT**);
static void node_put_point_set_to(NODE*,GLOBAL_POINT**);

static void get_point_value_from(POINT*,GLOBAL_POINT**);
static void put_point_value_to(POINT*,GLOBAL_POINT**);
static void set_node_impulse(ELASTIC_SET*,NODE*,GLOBAL_POINT**);
static void set_curve_impulse(ELASTIC_SET*,CURVE*,GLOBAL_POINT**);
static void set_surf_impulse(ELASTIC_SET*,SURFACE*,GLOBAL_POINT**);

static void reorder_string_curves(NODE*);
static void assembleParachuteSet3d(INTERFACE*,ELASTIC_SET*);

static void set_equilibrium_mesh2d(Front*);
static void set_equilibrium_mesh3d(Front*);


static void count_node_neighbors(
	NODE *node,
	SPRING_VERTEX *sv,
	int *n)
{
    int num_nb = 0;
	int dim = Dimension(node->interface);

	CURVE **c;
	for (c = node->out_curves; c && *c; ++c)
	{
	    if (dim == 3 && hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
	    num_nb++;
	}
	
    for (c = node->in_curves; c && *c; ++c)
	{
	    if (dim == 3 && hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
	    num_nb++;
	}

	sv[*n].ix = node->posn->indx = *n;
	
    if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    TRI *tri;

	    int num_tris = 0;
	    POINT *p = node->posn;

	    for (c = node->out_curves; c && *c; ++c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;

            BOND* b = (*c)->first;
            for (btris = Btris(b); btris && *btris; ++btris)
            {
                int nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
                for (int j = 0; j < nt; ++j)
                {
                    if (!pointer_in_list((POINTER)tris[j],num_tris,(POINTER*)tri_list))
                        tri_list[num_tris++] = tris[j];
                }
            }
        }
	    
        for (c = node->in_curves; c && *c; ++c)
        {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
    
            BOND* b = (*c)->last;
            for (btris = Btris(b); btris && *btris; ++btris)
            {
                int nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
                for (int j = 0; j < nt; ++j)
                {
                    if (!pointer_in_list((POINTER)tris[j],num_tris,(POINTER*)tri_list))
                        tri_list[num_tris++] = tris[j];
                }
            }
        }
	    
        for (int i = 0; i < num_tris; ++i)
        {
            tri = tri_list[i];
            for (int side = 0; side < 3; ++side)
            {
                if (p == Point_of_tri(tri)[side])
                {
                    if (is_side_bdry(tri,side)) continue;
                    num_nb++;
                }
            }
        }
	}

	sv[*n].num_nb = num_nb;
	(*n)++;
}	/* end count_node_neighbors */

static void count_curve_neighbors(
	CURVE *curve,
	SPRING_VERTEX *sv,
	int *n)
{
	if (curve->first == curve->last) return; // no interior point

	int i = *n;
	
    BOND *b;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    sv[i].num_nb = 2;
	    i++;
	}

	int dim = Dimension(curve->interface);
	if (dim == 2)
	{
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
            POINT *p = b->end;
            sv[i].ix = p->indx = i;
            i++;
	    }
	}

	if (dim == 3)
	{
	    i = *n;
	    for (b = curve->first; b != curve->last; b = b->next)
        {
            POINT* p = b->end;
            for (BOND_TRI** btris = Btris(b); btris && *btris; ++btris)
            {
	            TRI **tris;
                int nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
                for (int j = 0; j < nt; ++j)
                {
                    for (int side = 0; side < 3; ++side)
                    {
                        if (p == Point_of_tri(tris[j])[side])
                        {
                            if (is_side_bdry(tris[j],side)) continue;
                            sv[i].num_nb++;
                        }
                    }
                }
            }

            sv[i].ix = p->indx = i;
            i++;
        }
	}

	*n = i;
}	/* end count_curve_neighbors */

static void count_surf_neighbors(
	SURFACE *surf,
	SPRING_VERTEX *sv,
	int *n)
{
    int i = *n;
    TRI *tris[MAX_NUM_RING1];
	
    unsort_surf_point(surf);
	for (TRI* tri = first_tri(surf); !at_end_of_tri_list(tri,surf); tri = tri->next)
	{
	    for (int j = 0; j < 3; ++j)
	    {
            POINT* p = Point_of_tri(tri)[j];
            if (sorted(p) || Boundary_point(p)) continue;
            sorted(p) = YES;
            
            int nt;
            PointAndFirstRingTris(p,Hyper_surf_element(tri),Hyper_surf(surf),&nt,tris);
            sv[i].num_nb = nt;
            sv[i].ix = p->indx = i;
            i++;
	    }
	}

	*n = i;
}	/* end count_surf_neighbors */

extern void set_spring_vertex_memory(
	SPRING_VERTEX *sv,
	int size)
{
	for (int i = 0; i < size; ++i)
	{
	    if (sv[i].x_nb != nullptr)
        {
		    FT_FreeThese(5, sv[i].x_nb, sv[i].v_nb, sv[i].k,
                    sv[i].len0, sv[i].ix_nb);
        }

	    int num_nb = sv[i].num_nb;
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].x_nb,num_nb,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].v_nb,num_nb,sizeof(double*));
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].k,num_nb,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].len0,num_nb,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&sv[i].ix_nb,num_nb,sizeof(int));
	    
	    // reset external acceleration
        for (int j = 0; j < MAXD; ++j)
        {
		    sv[i].ext_accel[j] = 0.0;
        }
	}
}	/* end set_spring_vertex_memory */

extern void compute_spring_accel1(
	SPRING_VERTEX *sv,
	double *accel,
	int dim)
{
	for (int k = 0; k < dim; ++k)
	    accel[k] = 0.0;

	//acceleration due to elastic stretching force
    for (int j = 0; j < sv->num_nb; ++j)
	{
	    double v_rel[MAXD];
	    double vec[MAXD];
	    
        double len = 0.0;
	    for (int k = 0; k < dim; ++k)
	    {
            v_rel[k] = sv->v_nb[j][k] - sv->v[k];
            vec[k] = sv->x_nb[j][k] - sv->x[k];
            len += vec[k]*vec[k];
	    }
	    len = sqrt(len);

        double dL = len - sv->len0[j];
        
        //zero compressive stress
        if (dL <= 0.0) continue;

	    for (int k = 0; k < dim; ++k)
	    {
            accel[k] += sv->k[j]*dL*vec[k]/len/sv->m;
	    }
	}

    //acceleration due to elastic bending and damping forces
    for (int k = 0; k < dim; ++k)
    {
        accel[k] += sv->bendforce[k]/sv->m;
        accel[k] -= sv->lambda*(sv->v[k] - sv->ext_impul[k])/sv->m;
    }

    for (int k = 0; k < dim; ++k)
    {
        sv->f[k] = accel[k]*sv->m;
    }

    for (int k = 0; k < dim; ++k)
    {
        accel[k] += sv->ext_accel[k] + sv->fluid_accel[k] + sv->other_accel[k];
    }

}	/* end compute_spring_accel */


//RK4
void generic_spring_solver(
	SPRING_VERTEX *sv,
	int dim,
	int size,
	int n_loop,
	double dt)
{
	static double **x_old,**x_new,**v_old,**v_new,**accel;
	static int old_size = size;
	int i,j,n;

	if (debugging("trace"))
	    (void) printf("Entering generic_spring_solver()\n");
	if (x_old == NULL)
	{
	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&accel,size,3,sizeof(double));
	}
	if (size > old_size)
	{
	    FT_FreeThese(5, x_old, v_old, x_new, v_new, accel);
	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&accel,size,3,sizeof(double));
	    printf("size = %d, old_size = %d\n", size, old_size);
	}
	old_size = size;

    /*
//TEMP DEBUG
///////////////////////////////////////////////////////////////////////
for (i = 0; i < size; ++i)
{
    char msg[25];
    sprintf(msg,"sv[%d].bendforce",i);
    fprint_general_vector(stdout,msg,sv[i].bendforce,dim,"\n");
}
///////////////////////////////////////////////////////////////////////
    */	

    for (i = 0; i < size; ++i)
	{
	    compute_spring_accel1(&sv[i],accel[i],dim);
	}
	for (i = 0; i < size; ++i)
	for (j = 0; j < dim; ++j)
	{
	    x_old[i][j] = sv[i].x[j];
	    v_old[i][j] = sv[i].v[j];
	}
	
    /*
//TEMP DEBUG
///////////////////////////////////////////////////////////////////////
for (i = 0; i < size; ++i)
{
    char msg[25];
    sprintf(msg,"accel[%d]",i);
    fprint_general_vector(stdout,msg,accel[i],dim,"\n");
}
///////////////////////////////////////////////////////////////////////
    */
    
    for (n = 0; n < n_loop; ++n)
	{
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
            x_new[i][j] = x_old[i][j] + dt*v_old[i][j]/6.0;
            v_new[i][j] = v_old[i][j] + dt*accel[i][j]/6.0;
	    	sv[i].x[j] = x_old[i][j] + 0.5*v_old[i][j]*dt;
	    	sv[i].v[j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }
	    
        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/6.0;
	        }

	    for (i = 0; i < size; ++i)
	    {
		    compute_spring_accel1(&sv[i],accel[i],dim);
	    }
	    
        for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		    x_new[i][j] += dt*sv[i].v[j]/3.0;
            v_new[i][j] += dt*accel[i][j]/3.0;
	    	sv[i].x[j] = x_old[i][j] + 0.5*sv[i].v[j]*dt;
	    	sv[i].v[j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }
	    
        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/3.0;
	        }
	
	    for (i = 0; i < size; ++i)
	    {
		    compute_spring_accel1(&sv[i],accel[i],dim);
	    }

	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		    x_new[i][j] += dt*sv[i].v[j]/3.0;
            v_new[i][j] += dt*accel[i][j]/3.0;
	    	sv[i].x[j] = x_old[i][j] + sv[i].v[j]*dt;
	    	sv[i].v[j] = v_old[i][j] + accel[i][j]*dt; 
	    }

        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/3.0;
	        }

	    for (i = 0; i < size; ++i)
	    {
            compute_spring_accel1(&sv[i],accel[i],dim);
	    }

        for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		    x_new[i][j] += dt*sv[i].v[j]/6.0;
            v_new[i][j] += dt*accel[i][j]/6.0;
	    }

	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
            sv[i].x[j] = x_new[i][j];
            sv[i].v[j] = v_new[i][j];
	    }

        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/6.0;
            }

	    if (n != n_loop-1)
	    {
            for (i = 0; i < size; ++i)
                for (j = 0; j < 3; ++j)
                {
                    x_old[i][j] = sv[i].x[j];
                    v_old[i][j] = sv[i].v[j];
                    if (std::isnan(x_old[i][j]) || std::isinf(x_old[i][j]))
                    {
                        printf("After loop %d: x_old[%d][%d] = %f\n",n,i,j,x_old[i][j]);
                        LOC(); clean_up(ERROR);
                    }

                }
	    	for (i = 0; i < size; ++i)
            {
                compute_spring_accel1(&sv[i],accel[i],dim);
            }
	    }
	}
	if (debugging("trace"))
	    (void) printf("Leaving generic_spring_solver()\n");
}	/* end generic_spring_solver */

//TODO: Same as generic_spring_solver() -- after changing all instances of generic_spring_solver()
//      to spring_solver_RK4() can delete generic_spring_solver().
void spring_solver_RK4(
	SPRING_VERTEX *sv,
	int dim,
	int size,
	int n_loop,
	double dt)
{
	static double **x_old,**x_new,**v_old,**v_new,**accel;
	static int old_size = size;
	int i,j,n;

	if (debugging("trace"))
	    (void) printf("Entering spring_solver_RK4()\n");

	if (x_old == NULL)
	{
	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&accel,size,3,sizeof(double));
	}
	if (size > old_size)
	{
	    FT_FreeThese(5, x_old, v_old, x_new, v_new, accel);
	    FT_MatrixMemoryAlloc((POINTER*)&x_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_old,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&x_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&v_new,size,3,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&accel,size,3,sizeof(double));
	    printf("size = %d, old_size = %d\n", size, old_size);
	}
	old_size = size;

    for (i = 0; i < size; ++i)
	{
	    compute_spring_accel1(&sv[i],accel[i],dim);
	}

	for (i = 0; i < size; ++i)
	for (j = 0; j < dim; ++j)
	{
	    x_old[i][j] = sv[i].x[j];
	    v_old[i][j] = sv[i].v[j];
	}
    
    for (n = 0; n < n_loop; ++n)
	{
	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
            x_new[i][j] = x_old[i][j] + dt*v_old[i][j]/6.0;
            v_new[i][j] = v_old[i][j] + dt*accel[i][j]/6.0;
	    	sv[i].x[j] = x_old[i][j] + 0.5*v_old[i][j]*dt;
	    	sv[i].v[j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }
	    
        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/6.0;
	        }

	    for (i = 0; i < size; ++i)
	    {
		    compute_spring_accel1(&sv[i],accel[i],dim);
	    }
	    
        for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		    x_new[i][j] += dt*sv[i].v[j]/3.0;
            v_new[i][j] += dt*accel[i][j]/3.0;
	    	sv[i].x[j] = x_old[i][j] + 0.5*sv[i].v[j]*dt;
	    	sv[i].v[j] = v_old[i][j] + 0.5*accel[i][j]*dt;
	    }
	    
        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/3.0;
	        }
	
	    for (i = 0; i < size; ++i)
	    {
		    compute_spring_accel1(&sv[i],accel[i],dim);
	    }

	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		    x_new[i][j] += dt*sv[i].v[j]/3.0;
            v_new[i][j] += dt*accel[i][j]/3.0;
	    	sv[i].x[j] = x_old[i][j] + sv[i].v[j]*dt;
	    	sv[i].v[j] = v_old[i][j] + accel[i][j]*dt; 
	    }

        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/3.0;
	        }

	    for (i = 0; i < size; ++i)
	    {
            compute_spring_accel1(&sv[i],accel[i],dim);
	    }

        for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
		    x_new[i][j] += dt*sv[i].v[j]/6.0;
            v_new[i][j] += dt*accel[i][j]/6.0;
	    }

	    for (i = 0; i < size; ++i)
	    for (j = 0; j < dim; ++j)
	    {
            sv[i].x[j] = x_new[i][j];
            sv[i].v[j] = v_new[i][j];
	    }

        for (i = 0; i < size; ++i)
            for (j = 0; j < 3; ++j)
            {
                sv[i].ext_impul[j] += (sv[i].ext_accel[j] + 
					sv[i].fluid_accel[j])*dt/6.0;
            }

	    if (n != n_loop-1)
	    {
            for (i = 0; i < size; ++i)
                for (j = 0; j < 3; ++j)
                {
                    x_old[i][j] = sv[i].x[j];
                    v_old[i][j] = sv[i].v[j];
                    if (std::isnan(x_old[i][j]) || std::isinf(x_old[i][j]))
                    {
                        printf("After loop %d: x_old[%d][%d] = %f\n",n,i,j,x_old[i][j]);
                        LOC(); clean_up(ERROR);
                    }

                }
	    	for (i = 0; i < size; ++i)
            {
                compute_spring_accel1(&sv[i],accel[i],dim);
            }
	    }

	}//end n loop

	if (debugging("trace"))
	    (void) printf("Leaving spring_solver_RK4()\n");
}	/* end spring_solver_RK4 */

extern void count_vertex_neighbors(
	ELASTIC_SET *geom_set,
	SPRING_VERTEX *sv)
{
	if (debugging("canopy"))
	    (void) printf("Entering count_vertex_neighbors()\n");

	int ns = geom_set->num_surfs;
	int nc = geom_set->num_curves;
	int nn = geom_set->num_nodes;
	int n = 0;

	for (int i = 0; i < ns; ++i)
    {
	    count_surf_neighbors(geom_set->surfs[i],sv,&n);
    }
	for (int i = 0; i < nc; ++i)
    {
	    count_curve_neighbors(geom_set->curves[i],sv,&n);
    }
	for (int i = 0; i < nn; ++i)
    {
	    count_node_neighbors(geom_set->nodes[i],sv,&n);	
    }

	if (debugging("canopy"))
	    (void) printf("Leaving count_vertex_neighbors()\n");
}	/* end  count_vertex_neighbors */

extern void link_point_set(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store)
{
	if (debugging("canopy"))
	    (void) printf("Entering link_point_set()\n");

	int ns = geom_set->num_surfs;
	int nc = geom_set->num_curves;
	int nn = geom_set->num_nodes;
	int n = 0;

	for (int i = 0; i < ns; ++i)
    {
	    link_surf_point_set(geom_set,geom_set->surfs[i],point_set,
				point_set_store,&n);
    }
	for (int i = 0; i < nc; ++i)
	{
	    link_curve_point_set(geom_set,geom_set->curves[i],point_set,
				point_set_store,&n);
	}
	for (int i = 0; i < nn; ++i)
    {
	    link_node_point_set(geom_set,geom_set->nodes[i],point_set,
				point_set_store,&n);
    }

    int nrgbs = geom_set->num_rgb_surfs;
    for (int i = 0; i < nrgbs; ++i)
    {
        link_surf_point_set(geom_set,geom_set->rgb_surfs[i],point_set,
                point_set_store,&n);
    }

	if (debugging("canopy"))
	{
	    (void) printf("Final n = %d\n",n);
	    (void) printf("Leaving link_point_set()\n");
	}
}	/* end link_point_set */

static void link_surf_point_set(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store,
	int *n)
{
	TRI *tri;
	POINT *p;
	long gindex;
	int i,j;

	unsort_surf_point(surf);
	i = *n;
	surf_tri_loop(surf,tri)
	{
	    for (j = 0; j < 3; ++j)
	    {
            p = Point_of_tri(tri)[j];
            if (sorted(p) || Boundary_point(p)) continue;
            gindex = Gindex(p);
            point_set[gindex] = point_set_store + i;
            point_set[gindex]->gindex = gindex;
            sorted(p) = YES;
            i++;
	    }
	}
	*n = i;
}	/* end link_surf_point_set */

static void link_curve_point_set(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store,
	int *n)
{
	BOND *b;
	POINT *p;
	long gindex;
	int i = *n;

	if (curve->first == curve->last) return;	// Treated as nodes
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    p = b->end;
	    gindex = Gindex(p);
	    point_set[gindex] = point_set_store + i;
	    point_set[gindex]->gindex = gindex;
	    i++;
	}
	*n = i;
}	/* end link_curve_point_set */

static void link_node_point_set(
	ELASTIC_SET *geom_set,
	NODE *node,
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *point_set_store,
	int *n)
{
	long gindex;
	POINT *p = node->posn;

	gindex = Gindex(p);
	point_set[gindex] = point_set_store + (*n);
	point_set[gindex]->gindex = gindex;
	(*n)++;
}	/* end link_node_point_set */

extern void set_vertex_neighbors(
	ELASTIC_SET *geom_set,
	SPRING_VERTEX *sv,
	GLOBAL_POINT **point_set)
{
	if (debugging("canopy"))
	    (void) printf("Entering set_vertex_neighbors()\n");

	int ns = geom_set->num_surfs;
	int nc = geom_set->num_curves;
	int nn = geom_set->num_nodes;
	int n = 0;
	
    for (int i = 0; i < ns; ++i)
    {
	    set_surf_spring_vertex(geom_set,geom_set->surfs[i],sv,&n,
					point_set);
    }
	for (int i = 0; i < nc; ++i)
    {
	    set_curve_spring_vertex(geom_set,geom_set->curves[i],sv,&n,
					point_set);
    }
	for (int i = 0; i < nn; ++i)
    {
	    set_node_spring_vertex(geom_set,geom_set->nodes[i],sv,&n,
					point_set);
    }

	if (debugging("canopy"))
	    (void) printf("Leaving set_vertex_neighbors()\n");
}	/* end  set_vertex_neighbors */

extern void set_node_spring_vertex(
	ELASTIC_SET *geom_set,
	NODE *node,
	SPRING_VERTEX *sv,
	int *n,
	GLOBAL_POINT **point_set)
{
	CURVE **c;
	BOND *b;
	POINT *p,*p_nb;
	Front *front = geom_set->front;
	INTERFACE *intfc = front->interf;
	int i,j,nn,dim = Dimension(intfc);
	double ks = geom_set->ks;
	double ks_band = geom_set->ks_band;
	double kl = geom_set->kl;
	double kl_band = geom_set->kl_band;
	double kg = geom_set->kg;
	double mass = 0.0;
	double lambda_s = geom_set->lambda_s;
	double lambda_l = geom_set->lambda_l;
	double lambda_g = geom_set->lambda_g;
	boolean is_fixed = NO;
	STATE *sl,*sr;
	long gindex,gindex_nb;
	
    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double *g = nullptr;
	if (af_params)
 	    g = af_params->gravity;

    bool on_diskgap = false;
	if (dim == 3)
	{
	    AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    if (extra != NULL)
        {
            //TODO: Need to deal with case of FIXED_LOAD_NODE AND FIXED_RG_STRING_NODE
            if (extra->af_node_type == PRESET_NODE)
            {
                mass = geom_set->m_s;
                is_fixed = YES;
            }
            else if (extra->af_node_type == LOAD_NODE || 
                     extra->af_node_type == RG_STRING_NODE)
            {
                //TODO: Is this being set correctly?
                Front *front = geom_set->front;
                AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
                    
                mass = af_params->payload;

                if (extra->af_node_type == RG_STRING_NODE)
                {
                    mass = af_params->payload/af_params->num_rg_string_nodes;
                }
            }
            else if (extra->af_node_type == STRING_NODE)
                mass = geom_set->m_l;
            else if (extra->af_node_type == CANOPY_STRING_NODE ||
                     extra->af_node_type == DISKGAP_STRING_NODE)
            {
                mass = geom_set->m_s;
                if (extra->af_node_type == DISKGAP_STRING_NODE)
                    on_diskgap = true;
            }
            else if (extra->af_node_type == THR_LOAD_NODE)
                mass = geom_set->m_l;
            else if (extra->af_node_type == SEC_LOAD_NODE)
                mass = geom_set->m_l;
            else if (extra->af_node_type == GORE_NODE)
                mass = geom_set->m_g;
        }
	    else
        {
            mass = geom_set->m_s;
        }
	}
	else
	{//dim == 2
        mass = geom_set->m_l;
	    boolean on_canopy = NO;
	    node_out_curve_loop(node, c)
	    {
		    if (wave_type(*c) == ELASTIC_BOUNDARY)
		        on_canopy = YES;
	    }
	    
        node_in_curve_loop(node, c)
	    {
		    if (wave_type(*c) == ELASTIC_BOUNDARY)
		        on_canopy = YES;
	    }
	    
        if (on_canopy)
		    mass = geom_set->m_s;
	    
        AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
	    Front *front = geom_set->front;
        AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	    
        if (extra != NULL && extra->af_node_type == PRESET_NODE)
		    is_fixed = YES;
	    
        if (extra != NULL && extra->af_node_type == LOAD_NODE)
	    	mass = af_params->payload;
	}

	if (mass == 0.0)
	{
	    printf("\n\nERROR: mass is not set for some node\n\n");
        //////////////////
        AF_NODE_EXTRA *extra = (AF_NODE_EXTRA*)node->extra;
        printf("node->extra = %s\n",(extra == nullptr) ? "nullptr": "non-nullptr");
        print_node(node);
        //////////////////
	    LOC(); clean_up(ERROR);
	}

	nn = 0;
	gindex = Gindex(node->posn);
	sv[*n].x = point_set[gindex]->x;
	sv[*n].v = point_set[gindex]->v;
	sv[*n].f = point_set[gindex]->f;
    sv[*n].bendforce = point_set[gindex]->bendforce;
	sv[*n].ext_impul = point_set[gindex]->impuls;
	sv[*n].fluid_accel = point_set[gindex]->fluid_accel;
	sv[*n].other_accel = point_set[gindex]->other_accel;

	for (c = node->out_curves; c && *c; ++c)
	{
	    if (dim == 2 && wave_type(*c) == PASSIVE_HSBDRY) continue;
	    if (dim == 3 && hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
	    b = (*c)->first;
	    gindex_nb = Gindex(b->end);
	    sv[*n].x_nb[nn] = point_set[gindex_nb]->x;
	    sv[*n].v_nb[nn] = point_set[gindex_nb]->v;
	    sv[*n].len0[nn] = bond_length0(b);
	    sv[*n].m = mass;
	    sv[*n].ix_nb[nn] = b->end->indx; //important for gpu
	    
        if (dim == 3)
	    {
            if (is_fixed || is_load_node(node) || is_rg_string_node(node))
                sv[*n].k[nn] = 0.0;
            else if (hsbdry_type(*c) == STRING_HSBDRY)
                sv[*n].k[nn] = kl;
            else if (hsbdry_type(*c) == DISKGAP_STRING_HSBDRY)
                sv[*n].k[nn] = kl_band;
            else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
            {
                if (on_diskgap)
                    sv[*n].k[nn] = ks_band;
                else
                    sv[*n].k[nn] = ks;
            }
            else if (hsbdry_type(*c) == GORE_HSBDRY)
                sv[*n].k[nn] = kg;
            else if (hsbdry_type(*c) == FIXED_HSBDRY)
                is_fixed = YES;
	    }
	    else
        {
		    sv[*n].k[nn] = kl;
        }
	    ++nn;
	}

	for (c = node->in_curves; c && *c; ++c)
	{
	    if (dim == 2 && wave_type(*c) == PASSIVE_HSBDRY) continue;
	    if (dim == 3 && hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
	    b = (*c)->last;
	    gindex_nb = Gindex(b->start);
	    sv[*n].x_nb[nn] = point_set[gindex_nb]->x;
	    sv[*n].v_nb[nn] = point_set[gindex_nb]->v;
	    sv[*n].len0[nn] = bond_length0(b);
	    sv[*n].m = mass;
	    sv[*n].ix_nb[nn] = b->start->indx; //important for gpu
	    
        if (dim == 3)
	    {
            if (is_fixed || is_load_node(node) || is_rg_string_node(node))
                sv[*n].k[nn] = 0.0;
            else if (hsbdry_type(*c) == STRING_HSBDRY)
                sv[*n].k[nn] = kl;
            else if (hsbdry_type(*c) == DISKGAP_STRING_HSBDRY)
                sv[*n].k[nn] = kl_band;
            else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
            {
                if (on_diskgap)
                    sv[*n].k[nn] = ks_band;
                else
                    sv[*n].k[nn] = ks;
            }
            else if (hsbdry_type(*c) == GORE_HSBDRY)
                sv[*n].k[nn] = kg;
            else if (hsbdry_type(*c) == FIXED_HSBDRY)
                is_fixed = YES;
	    }
	    else
        {
            sv[*n].k[nn] = kl;
        }
	    ++nn;
	}
	
    if (dim == 3)
	{
	    BOND_TRI **btris;
	    TRI **tris,*tri_list[500];
	    int k,side,nt,num_tris;
	    TRI *tri;

	    num_tris = 0;
	    p = node->posn;
	    for (c = node->out_curves; c && *c; ++c)
	    {
            if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
            b = (*c)->first;
            for (btris = Btris(b); btris && *btris; ++btris)
            {
                nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
                for (j = 0; j < nt; ++j)
                {
                    if (!pointer_in_list((POINTER)tris[j],num_tris, (POINTER*)tri_list))
                        tri_list[num_tris++] = tris[j];
                }
            }
	    }

	    for (c = node->in_curves; c && *c; ++c)
	    {
		    if (hsbdry_type(*c) == PASSIVE_HSBDRY) continue;
		    b = (*c)->last;
            for (btris = Btris(b); btris && *btris; ++btris)
            {
                nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
                for (j = 0; j < nt; ++j)
                {
                    if (!pointer_in_list((POINTER)tris[j],num_tris, (POINTER*)tri_list))
                        tri_list[num_tris++] = tris[j];
                }
            }
	    }

	    for (i = 0; i < num_tris; ++i)
        {
            tri = tri_list[i];
            for (side = 0; side < 3; ++side)
            {
                if (p == Point_of_tri(tri)[side])
                {
                    if (is_side_bdry(tri,side)) continue;
                    p_nb = Point_of_tri(tri)[(side+1)%3];
                    gindex_nb = Gindex(p_nb);
                    sv[*n].x_nb[nn] = point_set[gindex_nb]->x;
                    sv[*n].v_nb[nn] = point_set[gindex_nb]->v;
                    sv[*n].ix_nb[nn] = p_nb->indx;
                    if (on_diskgap)
                        sv[*n].k[nn] = ks_band;
                    else
                        sv[*n].k[nn] = ks;
                    if (is_fixed) sv[*n].k[nn] = 0.0;
                    sv[*n].len0[nn] = tri->side_length0[side];
                    ++nn;
                }
            }
        }
	    
        //TODO: should load_nodes and rg_string_nodes be handled differently than a fixed node?
        sv[*n].lambda = lambda_s;
        if (is_fixed || is_load_node(node) || is_rg_string_node(node)) 
	    {
		    sv[*n].lambda = 0.0;
	    	for (i = 0; i < sv[*n].num_nb; ++i)
		        sv[*n].k[i] = 0.0;
	    }
	
    }
	else
	{
	    sv[*n].lambda = lambda_l;
	    if (is_fixed)
        {
            sv[*n].lambda = 0.0;
            for (i = 0; i < sv[*n].num_nb; ++i)
                sv[*n].k[i] = 0.0;
        }
	}
        
    for (i = 0; i < dim; ++i)
    {
	    if (is_fixed || g == NULL)
	    	sv[*n].ext_accel[i] = 0;
	    else
	    	sv[*n].ext_accel[i] = g[i];
	}

	(*n)++;
}	/* end set_node_spring_vertex */

extern void set_curve_spring_vertex(
	ELASTIC_SET *geom_set,
	CURVE *curve,
	SPRING_VERTEX *sv,
	int *n,
	GLOBAL_POINT **point_set)
{
	Front *front = geom_set->front;
	int i,j,nn;
	BOND *b;
	double kl,m_l,lambda_l;
	int dim = front->rect_grid->dim;
	long gindex,gindex_nb;

    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double *g = nullptr;
	if (af_params)
 	    g = af_params->gravity;

	if (dim == 3)
	{
	    if (hsbdry_type(curve) == STRING_HSBDRY)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
        else if (hsbdry_type(curve) == DISKGAP_STRING_HSBDRY)
        {
	    	kl = geom_set->kl_band;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
        }
	    else if (hsbdry_type(curve) == GORE_HSBDRY)
	    {
	    	kl = geom_set->kg;
	    	m_l = geom_set->m_g;
	    	lambda_l = geom_set->lambda_g;
	    }
	    else if (hsbdry_type(curve) == FIXED_HSBDRY)
	    {
	    	kl = 0.0;
	    	m_l = geom_set->m_l;
	    	lambda_l = 0.0;
	    }
	    else //MONO_COMP_HSBDRY
	    {
	    	kl = geom_set->ks;
	    	m_l = geom_set->m_s;
	    	lambda_l = geom_set->lambda_s;

            b = curve->first;
            BOND_TRI **btris = Btris(b);
            HYPER_SURF* hs = Hyper_surf((*btris)->surface);
            if (wave_type(hs) == ELASTIC_BAND_BOUNDARY)
                kl = geom_set->ks_band;
	    }
	}
	else if (dim == 2)
	{
	    if (wave_type(curve) == ELASTIC_STRING)
	    {
	    	kl = geom_set->kl;
	    	m_l = geom_set->m_l;
	    	lambda_l = geom_set->lambda_l;
	    }
	    else if (wave_type(curve) == ELASTIC_BOUNDARY)
	    {
		    kl = geom_set->ks;
            m_l = geom_set->m_s;
            lambda_l = geom_set->lambda_s;
	    }
	}

	i = *n;
	for (b = curve->first; b != curve->last; b = b->next)
	{
	    gindex = Gindex(b->end);
	    sv[i].x = point_set[gindex]->x;
	    sv[i].v = point_set[gindex]->v;
	    sv[i].f = point_set[gindex]->f;
        sv[i].bendforce = point_set[gindex]->bendforce;
	    sv[i].ext_impul = point_set[gindex]->impuls;
	    sv[i].fluid_accel = point_set[gindex]->fluid_accel;
	    sv[i].other_accel = point_set[gindex]->other_accel;
	    
        gindex_nb = Gindex(b->start);
	    sv[i].x_nb[0] = point_set[gindex_nb]->x;
	    sv[i].v_nb[0] = point_set[gindex_nb]->v;
	    
        gindex_nb = Gindex(b->next->end);
	    sv[i].x_nb[1] = point_set[gindex_nb]->x;
	    sv[i].v_nb[1] = point_set[gindex_nb]->v;
	    
        sv[i].ix_nb[0] = b->start->indx;
	    sv[i].ix_nb[1] = b->next->end->indx;
	    sv[i].len0[0] = bond_length0(b);
	    sv[i].len0[1] = bond_length0(b->next);
	    sv[i].k[0] = sv[i].k[1] = kl;
	    sv[i].m = m_l;
	    sv[i].num_nb = 2;
	    sv[i].lambda = lambda_l;
	
        if (dim == 3)
	    {
	    	if (hsbdry_type(curve) == FIXED_HSBDRY || g == NULL)
	    	{
                for (j = 0; j < dim; ++j)
                    sv[i].ext_accel[j] = 0;
	    	}
            else
            {
                for (j = 0; j < dim; ++j)
                    sv[i].ext_accel[j] = g[j];
            }
	    }
	    else if (dim == 2 && g)
	    {
                for (j = 0; j < dim; ++j)
                	sv[i].ext_accel[j] = g[j];
	    }

	    SURFACE** surf;
	    boolean is_stationary_point = NO;
	    intfc_surface_loop(front->interf, surf)
	    {
		    is_stationary_point = is_registered_point(*surf, b->end);
		    if (is_stationary_point) break;
	    }
	    
        if (is_stationary_point)
	    {
		    sv[i].k[0] = sv[i].k[1] = 0.0;
		    sv[i].lambda = 0.0;
		    for (j = 0; j < dim; ++j)
		        sv[i].ext_accel[j] = 0.0;
	    }
	    
        ++i;
	}

	if (dim == 3)
	{
	    POINT *p,*p_nb;
	    BOND_TRI **btris;
	    TRI **tris;
	    int j,k,side,nt;
	    double ks;

	    if (hsbdry_type(curve) == FIXED_HSBDRY)
		    ks = 0.0;
	    else
		    ks = geom_set->ks;
	    
        i = *n;
	    
        for (b = curve->first; b != curve->last; b = b->next)
	    {
		    SURFACE** surf;
		    boolean is_stationary_point = NO;
            intfc_surface_loop(front->interf, surf)
            {
                is_stationary_point = is_registered_point(*surf, b->end);
                if (is_stationary_point) break;
            }
            
            p = b->end;
            nn = sv[i].num_nb;
            sv[i].m = m_l;
            
            for (btris = Btris(b); btris && *btris; ++btris)
            {
                HYPER_SURF* hs = Hyper_surf((*btris)->surface);
                if (wave_type(hs) == ELASTIC_BAND_BOUNDARY)
                    ks = geom_set->ks_band;

                nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
                for (j = 0; j < nt; ++j)
                {
                    for (side = 0; side < 3; ++side)
                    {
                        if (p == Point_of_tri(tris[j])[side])
                        {
                            if (is_side_bdry(tris[j],side)) continue;
                            p_nb = Point_of_tri(tris[j])[(side+1)%3];
                            gindex_nb = Gindex(p_nb);
                            sv[i].x_nb[nn] = point_set[gindex_nb]->x;
                            sv[i].v_nb[nn] = point_set[gindex_nb]->v;
                            sv[i].ix_nb[nn] = p_nb->indx;
                            sv[i].k[nn] = ks;
                            sv[i].len0[nn] = tris[j]->side_length0[side];
                            
                            if (is_stationary_point)
                            {
                                sv[i].k[nn] = 0.0;
                            }
                            
                            ++nn;
                        }
                    }
                }
            }

            sv[i].num_nb = nn;
            i++;
	    }
	}
	
    *n = i;
}	/* end set_curve_spring_vertex */

extern void set_surf_spring_vertex(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
	SPRING_VERTEX *sv,
	int *n,
	GLOBAL_POINT **point_set)
{
	Front *front = geom_set->front;
	int i,j,k,l,nt;
	TRI *tri;
	TRI *tris[MAX_NUM_RING1];
	POINT *p,*p_nb;
	double ks = geom_set->ks;
	double ks_band = geom_set->ks_band;
	double m_s = geom_set->m_s;
	double lambda_s = geom_set->lambda_s;
	boolean is_stationary_point;
	int dim = front->rect_grid->dim;
	STATE *sl,*sr;
	HYPER_SURF_ELEMENT *hse;
        HYPER_SURF         *hs = Hyper_surf(surf);
	long gindex,gindex_nb;

    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double *g = nullptr;
	if (af_params)
 	    g = af_params->gravity;

    if (wave_type(surf) == ELASTIC_BAND_BOUNDARY)
    {
        ks = ks_band;
    }

	unsort_surf_point(surf);
	i = *n;
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		
        sv[i].m = m_s;
		sv[i].lambda = lambda_s;
		
        is_stationary_point = is_registered_point(surf,p);
		if (is_stationary_point == YES)
		    sv[i].lambda = 0.0;

        for (k = 0; k < dim; ++k)
        {
		    if (is_stationary_point == YES || g == NULL)
	    	    	sv[i].ext_accel[k] = 0.0;
		    else
			    sv[i].ext_accel[k] = g[k];
        }

		gindex = Gindex(p);
		sv[i].x = point_set[gindex]->x;
		sv[i].v = point_set[gindex]->v;
		sv[i].f = point_set[gindex]->f;

        sv[i].bendforce = point_set[gindex]->bendforce;

        sv[i].ext_impul = point_set[gindex]->impuls;
        sv[i].fluid_accel = point_set[gindex]->fluid_accel;
        sv[i].other_accel = point_set[gindex]->other_accel;
		
        PointAndFirstRingTris(p,Hyper_surf_element(tri),
				Hyper_surf(surf),&nt,tris);
		sv[i].num_nb = nt;
		for (k = 0; k < nt; ++k)
		for (l = 0; l < 3; ++l)
		if (Point_of_tri(tris[k])[l] == p)
		{
		    p_nb = Point_of_tri(tris[k])[(l+1)%3];
		    gindex_nb = Gindex(p_nb);
		    sv[i].x_nb[k] = point_set[gindex_nb]->x;
		    sv[i].v_nb[k] = point_set[gindex_nb]->v;
		    sv[i].ix_nb[k] = p_nb->indx;

		    if (is_stationary_point == YES)
		    	sv[i].k[k] = 0.0;
		    else
            {
                sv[i].k[k] = ks;
            }

            sv[i].len0[k] = tris[k]->side_length0[l];
		}
		sorted(p) = YES;
	    	++i;
	    }
	}
	*n = i;
}	/* end set_surf_spring_vertex */

static void get_point_value_from(
	POINT *p,
	GLOBAL_POINT **point_set)
{
	int i;
	STATE *state = (STATE*)left_state(p);
	long gindex = Gindex(p);

	point_set[gindex]->gindex = gindex;
	for (i = 0; i < 3; ++i)
	{
	    point_set[gindex]->x[i] = Coords(p)[i] - p->pshift[i];
	    point_set[gindex]->v[i] = p->vel[i];
	    point_set[gindex]->f[i] = p->force[i];
        point_set[gindex]->bendforce[i] = state->bendforce[i];
	    point_set[gindex]->impuls[i] = state->impulse[i];
	    point_set[gindex]->fluid_accel[i] = state->fluid_accel[i];
	    point_set[gindex]->other_accel[i] = state->other_accel[i];
	}
}	/* end get_point_value_from */
	
static void put_point_value_to(
	POINT *p,
	GLOBAL_POINT **point_set)
{
	int i;
	long gindex = Gindex(p);
    
    //TODO: Can we set both the left and right states here? 
    //      If so we can eliminate the call to set_vertex_impulse(),
    //      since all it does is set the left and right state's impulse.
	/*
    STATE *sl = (STATE*)left_state(p);
	STATE *sr = (STATE*)right_state(p);
	*/

    STATE *state = (STATE*)left_state(p);

	for (i = 0; i < 3; ++i)
	{
	    Coords(p)[i] = point_set[gindex]->x[i] + p->pshift[i];
	    p->vel[i] = point_set[gindex]->v[i];
	    p->force[i] = point_set[gindex]->f[i];
	    state->impulse[i] = point_set[gindex]->impuls[i];
	    /*
        sl->impulse[i] = point_set[gindex]->impuls[i];
	    sr->impulse[i] = point_set[gindex]->impuls[i];
        */
	}
}	/* end put_point_value_to */
	
extern void get_point_set_from(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set)
{
	if (debugging("canopy"))
	    (void) printf("Entering get_point_set_from()\n");

	int ns = geom_set->num_surfs;
	int nc = geom_set->num_curves;
	int nn = geom_set->num_nodes;

	for (int i = 0; i < ns; ++i)
	    surf_get_point_set_from(geom_set->surfs[i],point_set);
	for (int i = 0; i < nc; ++i)
	    curve_get_point_set_from(geom_set->curves[i],point_set);
	for (int i = 0; i < nn; ++i)
	    node_get_point_set_from(geom_set->nodes[i],point_set);

    int nrgbs = geom_set->num_rgb_surfs;
    for (int i = 0; i < nrgbs; ++i)
        surf_get_point_set_from(geom_set->rgb_surfs[i],point_set);

	if (debugging("canopy"))
	    (void) printf("Leaving get_point_set_from()\n");
}	/* end  get_point_set_from */

extern void put_point_set_to(
	ELASTIC_SET *geom_set,
	GLOBAL_POINT **point_set)
{
	if (debugging("canopy"))
	    (void) printf("Entering put_point_set_to()\n");

	int ns = geom_set->num_surfs;
	int nc = geom_set->num_curves;
	int nn = geom_set->num_nodes;

	for (int i = 0; i < ns; ++i)
	    surf_put_point_set_to(geom_set->surfs[i],point_set);
	
    for (int i = 0; i < nc; ++i)
	    curve_put_point_set_to(geom_set->curves[i],point_set);
	
    for (int i = 0; i < nn; ++i)
	    node_put_point_set_to(geom_set->nodes[i],point_set);

    int nrgbs = geom_set->num_rgb_surfs;
    for (int i = 0; i < nrgbs; ++i)
        surf_put_point_set_to(geom_set->rgb_surfs[i],point_set);

	if (debugging("canopy"))
	    (void) printf("Leaving put_point_set_to()\n");
}	/* end  put_point_set_to */

static void surf_get_point_set_from(
	SURFACE *surf,
	GLOBAL_POINT **point_set)
{
	TRI *tri;
	POINT *p;
	int j;

	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		if (sorted(p) || Boundary_point(p)) continue;
		get_point_value_from(p,point_set);
		sorted(p) = YES;
	    }
	}
}	/* end surf_get_point_set_from */

static void curve_get_point_set_from(
	CURVE *curve,
	GLOBAL_POINT **point_set)
{
	BOND *b;
	POINT *p;

	for (b = curve->first; b != curve->last; b = b->next)
	{
	    p = b->end;
	    get_point_value_from(p,point_set);
	}
}	/* end curve_get_point_set_from */

static void node_get_point_set_from(
	NODE *node,
	GLOBAL_POINT **point_set)
{
	POINT *p = node->posn;
	get_point_value_from(p,point_set);
}	/* end node_get_point_set_from */

static void surf_put_point_set_to(
	SURFACE *surf,
	GLOBAL_POINT **point_set)
{
	TRI *tri;
	POINT *p;
	int j;

	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
	{
	    for (j = 0; j < 3; ++j)
	    {
            p = Point_of_tri(tri)[j];
            if (sorted(p) || Boundary_point(p)) continue;
            put_point_value_to(p,point_set);
            sorted(p) = YES;
	    }
	}
}	/* end surf_put_point_set_to */

static void curve_put_point_set_to(
	CURVE *curve,
	GLOBAL_POINT **point_set)
{
	BOND *b;
	POINT *p;

	for (b = curve->first; b != curve->last; b = b->next)
	{
	    p = b->end;
	    put_point_value_to(p,point_set);
	}
}	/* end curve_put_point_set_to */

static void node_put_point_set_to(
	NODE *node,
	GLOBAL_POINT **point_set)
{
	POINT *p = node->posn;
	put_point_value_to(p,point_set);
}	/* end node_put_point_set_to */

extern void set_elastic_params(
	ELASTIC_SET *geom_set,
	double fr_dt)
{
	Front *front = geom_set->front;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	/* Set elastic set kinetic parameters */
    geom_set->ks = af_params->ks;
    geom_set->lambda_s = af_params->lambda_s;
    geom_set->m_s = af_params->m_s;
    geom_set->kl = af_params->kl;
    geom_set->lambda_l = af_params->lambda_l;
    geom_set->m_l = af_params->m_l;
    geom_set->kg = af_params->kg;
    geom_set->lambda_g = af_params->lambda_g;
    geom_set->m_g = af_params->m_g;

    if (debugging("rigid_canopy"))
    {
	    geom_set->dt_tol = HUGE;
        return;
    }

	/* Set elastic set time step */
    double dt_tol = sqrt((af_params->m_s)/(af_params->ks))/10.0;
    
    /*
    const int Ms = 6;
    dt_tol = 0.1*sqrt(0.5*af_params->m_s/(Ms*af_params->ks));
    */

    if (af_params->strings_present &&
        dt_tol > sqrt((af_params->m_l)/(af_params->kl))/10.0)
    {
        dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
    }
    
    if (af_params->gores_present &&
        dt_tol > sqrt((af_params->m_g)/(af_params->kg))/10.0)
    {
        dt_tol = sqrt((af_params->m_g)/(af_params->kg))/10.0;
    }

	pp_global_min(&dt_tol,1);
	geom_set->dt_tol = dt_tol;
}	/* end set_elastic_params */

extern void set_elastic_params(
	ELASTIC_SET *geom_set,
    AF_PARAMS* af_params,
	double fr_dt)
{
	/* Set elastic set kinetic parameters */
    geom_set->ks = af_params->ks;
    geom_set->lambda_s = af_params->lambda_s;
    geom_set->m_s = af_params->m_s;
    geom_set->kl = af_params->kl;
    geom_set->lambda_l = af_params->lambda_l;
    geom_set->m_l = af_params->m_l;
    geom_set->kg = af_params->kg;
    geom_set->lambda_g = af_params->lambda_g;
    geom_set->m_g = af_params->m_g;

    geom_set->ks_band = af_params->ks_band;
    geom_set->kl_band = af_params->kl_band;

	/* Set elastic set time step */
	double dt_tol = sqrt((af_params->m_s)/(af_params->ks))/10.0;
    
    if (af_params->strings_present &&
        dt_tol > sqrt((af_params->m_l)/(af_params->kl))/10.0)
    {
        dt_tol = sqrt((af_params->m_l)/(af_params->kl))/10.0;
    }
    
    if (af_params->gores_present &&
        dt_tol > sqrt((af_params->m_g)/(af_params->kg))/10.0)
    {
        dt_tol = sqrt((af_params->m_g)/(af_params->kg))/10.0;
    }

	pp_global_min(&dt_tol,1);
	geom_set->dt_tol = dt_tol;
}	/* end set_elastic_params */

extern void merge_global_point_set(
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *gpoint_store,
	int num_gpoint)
{
	int i,k,gindex;
	for (i = 0; i < num_gpoint; ++i)
	{
	    gindex = gpoint_store[i].gindex;
	    *point_set[gindex] = gpoint_store[i];
	}
}	/* end merge_global_point_set */

extern void assembleParachuteSet(
	INTERFACE *intfc,
	ELASTIC_SET *geom_set)
{
    assembleParachuteSet3d(intfc,geom_set);
}	/* end assembleParachuteSet */

static void assembleParachuteSet3d(
	INTERFACE *intfc,
	ELASTIC_SET *geom_set)
{
	SURFACE **s = NULL;
	CURVE **c = NULL;
	NODE **n = NULL;

    SURFACE **rgb_surfs = geom_set->rgb_surfs;
    SURFACE **surfs = geom_set->surfs;
	CURVE **curves = geom_set->curves;
	NODE **nodes = geom_set->nodes;

    INTERFACE *cur_intfc = current_interface();
    set_current_interface(intfc);

	/* Assemble canopy surfaces */
	
    int nrgbs = 0;
    int ns = 0;
    int nc = 0;
    int nn = 0;
	
    intfc_surface_loop(intfc,s)
	{
        if (is_bdry(*s)) continue;

        if (wave_type(*s) == NEUMANN_BOUNDARY ||
            wave_type(*s) == MOVABLE_BODY_BOUNDARY)
        {
            rgb_surfs[nrgbs++] = *s;
        }

	    if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_BAND_BOUNDARY) continue;

	    surfs[ns++] = *s;
	    surf_pos_curve_loop(*s,c)
	    {
		if (hsbdry_type(*c) == SUBDOMAIN_HSBDRY) continue;
	    	if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    	{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
	    	}
	    }
	    surf_neg_curve_loop(*s,c)
	    {
		if (hsbdry_type(*c) == SUBDOMAIN_HSBDRY) continue;
	    	if (!pointer_in_list(*c,nc,(POINTER*)curves))
	    	{
		    curves[nc++] = *c;
		    if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->start;
		    if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	nodes[nn++] = (*c)->end;
	    	}
	    }
	}
	
    //TODO: below if finished?
    /* Change for dealing the cases where there is both canopy surface
	   and isolated 3d curves */
	intfc_curve_loop(intfc,c)
	{
	    if (pointer_in_list(*c,nc,(POINTER*)curves)) continue;
	    if (hsbdry_type(*c) == STRING_HSBDRY ||
            hsbdry_type(*c) == DISKGAP_STRING_HSBDRY ||
	        hsbdry_type(*c) == MONO_COMP_HSBDRY ||
	        hsbdry_type(*c) == GORE_HSBDRY)
	    {
	        curves[nc++] = *c;
	        if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
                	nodes[nn++] = (*c)->start;
                if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
                	nodes[nn++] = (*c)->end;
	    }
	}	


	/* Assemble curves and nodes */
	
    int num_layers = 3; //TODO: what does this mean???
	for (int l = 0; l < num_layers; ++l)
	{
	    for (int i = 0; i < nn; ++i)
	    {
	    	node_in_curve_loop(nodes[i],c)
	    	{
		    if (hsbdry_type(*c) == PASSIVE_HSBDRY ||
			hsbdry_type(*c) == SUBDOMAIN_HSBDRY) 
			continue;
		    if (!pointer_in_list(*c,nc,(POINTER*)curves))
		    {
		    	curves[nc++] = *c;
		    	if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->start;
		    	if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->end;
		    }
	    	}
	    	node_out_curve_loop(nodes[i],c)
	    	{
		    if (hsbdry_type(*c) == PASSIVE_HSBDRY ||
			hsbdry_type(*c) == SUBDOMAIN_HSBDRY) 
			continue;
		    if (!pointer_in_list(*c,nc,(POINTER*)curves))
		    {
		    	curves[nc++] = *c;
		    	if (!pointer_in_list((*c)->start,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->start;
		    	if (!pointer_in_list((*c)->end,nn,(POINTER*)nodes))
		    	    nodes[nn++] = (*c)->end;
		    }
	    	}
	    }
	}

	geom_set->num_surfs = ns;
	geom_set->num_curves = nc;
	geom_set->num_nodes = nn;
	geom_set->elastic_num_verts = 0;
	
    for (int i = 0; i < ns; ++i)
	    geom_set->elastic_num_verts += I_NumOfSurfInteriorPoints(surfs[i]);
	for (int i = 0; i < nc; ++i)
	    geom_set->elastic_num_verts += I_NumOfCurveInteriorPoints(curves[i]);
	geom_set->elastic_num_verts += nn;
	
 
    geom_set->num_rgb_surfs = nrgbs;
    geom_set->total_num_verts = geom_set->elastic_num_verts;

    for (int i = 0; i < nrgbs; ++i)
        geom_set->total_num_verts += I_NumOfSurfInteriorPoints(rgb_surfs[i]);


    int n_rgbsn = 0;
    geom_set->load_node = NULL;

	for (int i = 0; i < nn; ++i)
	{
	    if (is_load_node(nodes[i]) || is_rg_string_node(nodes[i]))
	    {
            if (is_load_node(nodes[i]))
                geom_set->load_node = nodes[i];
            else
                geom_set->rg_string_nodes[n_rgbsn++] = nodes[i];
            reorder_string_curves(nodes[i]);
            //TODO: Does the node pointer become invalid after reorder string curves?
            //      Or at least an incomplete copy of the original node?
	    }
	}

    if (debugging("intfc_assembly"))
    {
        printf("ns = %d, nc = %d, nn = %d, elastic_num_verts = %d\n",
                ns, nc, nn, geom_set->elastic_num_verts);
        printf("nrgbs = %d, total_num_verts = %d\n", geom_set->total_num_verts);
    }

    set_current_interface(cur_intfc);
}	/* end assembleParachuteSet3d */

extern void copy_from_client_point_set(
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *client_point_set,
	int client_size,
	double *client_L,
	double *client_U)
{
	for (int j = 0; j < client_size; ++j)
    {
        GLOBAL_POINT cpoint = client_point_set[j];

        bool skip_copy = false;
        for (int k = 0; k < 3; ++k)
        {
            if (cpoint.x[k] < client_L[k] || cpoint.x[k] >= client_U[k])
            {
                skip_copy = true;
                break;
            }
        }    
        if (skip_copy) continue;

        long gindex = cpoint.gindex;
        for (int k = 0; k < 3; ++k)
        {
            point_set[gindex]->x[k] = cpoint.x[k];
            point_set[gindex]->v[k] = cpoint.v[k];
            point_set[gindex]->f[k] = cpoint.f[k];
            point_set[gindex]->bendforce[k] = cpoint.bendforce[k];
            point_set[gindex]->impuls[k] = cpoint.impuls[k];
            point_set[gindex]->fluid_accel[k] = cpoint.fluid_accel[k];
            point_set[gindex]->other_accel[k] = cpoint.other_accel[k];
        }
    }
}	/* end copy_from_client_point_set */

extern void copy_to_client_point_set(
	GLOBAL_POINT **point_set,
	GLOBAL_POINT *client_point_set,
	int client_size)
{
	int j,k;
	long gindex;
	for (j = 0; j < client_size; j++)
        {
            gindex = client_point_set[j].gindex;
            for (k = 0; k < 3; k++)
            {
                client_point_set[j].x[k] = point_set[gindex]->x[k];
                client_point_set[j].v[k] = point_set[gindex]->v[k];
                client_point_set[j].f[k] = point_set[gindex]->f[k];
                client_point_set[j].bendforce[k] = point_set[gindex]->bendforce[k];
                client_point_set[j].impuls[k] = point_set[gindex]->impuls[k];
                client_point_set[j].fluid_accel[k] = point_set[gindex]->fluid_accel[k];
                client_point_set[j].other_accel[k] = point_set[gindex]->other_accel[k];
            }
        }
}	/* end copy_to_client_point_set */

static void reorder_string_curves(NODE *node)
{
	CURVE **c,**string_curves,*c_tmp;
	int i,j,num_curves;
	POINT **nb_points,*p_tmp;
    
    INTERFACE *save_intfc = current_interface();
    set_current_interface(node->interface);

	num_curves = I_NumOfNodeCurves(node);
	FT_VectorMemoryAlloc((POINTER*)&string_curves,num_curves,
				sizeof(CURVE*));
	FT_VectorMemoryAlloc((POINTER*)&nb_points,num_curves,
				sizeof(POINT*));
	i = 0;	
	node_in_curve_loop(node,c)
	{
	    string_curves[i] = *c;
	    nb_points[i++] = (*c)->last->start;
	}
	node_out_curve_loop(node,c)
	{
	    string_curves[i] = *c;
	    nb_points[i++] = (*c)->first->end;
	}
	for (i = 0; i < num_curves; ++i)
	{
	    if (string_curves[i]->end == node)
		delete_from_pointers(string_curves[i],&node->in_curves);
	    if (string_curves[i]->start == node)
		delete_from_pointers(string_curves[i],&node->out_curves);
	}
	for (i = 0; i < num_curves-1; ++i)
	for (j = i+1; j < num_curves; ++j)
	{
	    if (Gindex(nb_points[i]) > Gindex(nb_points[j]))
	    {
		c_tmp = string_curves[i];
		string_curves[i] = string_curves[j];
		string_curves[j] = c_tmp;
		p_tmp = nb_points[i];
		nb_points[i] = nb_points[j];
		nb_points[j] = p_tmp;
	    }	
	}
	for (i = 0; i < num_curves; ++i)
	{
	    if (string_curves[i]->end == node)
		unique_add_to_pointers(string_curves[i],&node->in_curves);
	    if (string_curves[i]->start == node)
		unique_add_to_pointers(string_curves[i],&node->out_curves);
	}

	FT_FreeThese(2,string_curves,nb_points);
    set_current_interface(save_intfc);
}	/* end reorder_string_curves */

//TODO: it seems that we should be able to perform the
//      assignments of this function inside put_point_value_to()
//      called by the XXX_put_point_set_to() functions
extern void set_vertex_impulse(
        ELASTIC_SET *geom_set,
        GLOBAL_POINT **point_set)
{
        int i,ns,nc,nn;

        ns = geom_set->num_surfs;
        nc = geom_set->num_curves;
        nn = geom_set->num_nodes;
        for (i = 0; i < ns; ++i)
            set_surf_impulse(geom_set,geom_set->surfs[i],point_set);
        for (i = 0; i < nc; ++i)
            set_curve_impulse(geom_set,geom_set->curves[i],point_set);
        for (i = 0; i < nn; ++i)
            set_node_impulse(geom_set,geom_set->nodes[i],point_set);

}       /* end set_vertex_impulse */

static void set_node_impulse(
	ELASTIC_SET *geom_set,
	NODE *node,
        GLOBAL_POINT **point_set)
{
	int i,dim;
	STATE *sl,*sr;
	long gindex = Gindex(node->posn);

	dim = FT_Dimension();
	sl = (STATE*)left_state(node->posn);
	sr = (STATE*)right_state(node->posn);

	for (i = 0; i < dim; ++i)
    {
	    sl->impulse[i] = point_set[gindex]->impuls[i];
        sr->impulse[i] = point_set[gindex]->impuls[i];
    }
}	/* end set_node_impulse */

static void set_curve_impulse(
	ELASTIC_SET *geom_set,
	CURVE *curve,
        GLOBAL_POINT **point_set)
{
	int j,dim;
	STATE *sl,*sr;
	BOND *b;
	long gindex;

	dim = FT_Dimension();

	for (b = curve->first; b != curve->last; b = b->next)
	{
	    gindex = Gindex(b->end);
	    sl = (STATE*)left_state(b->end);
	    sr = (STATE*)right_state(b->end);

        for (j = 0; j < dim; ++j)
        {
	    	sl->impulse[j] = point_set[gindex]->impuls[j]; 
            sr->impulse[j] = point_set[gindex]->impuls[j];
	    }
	}
}	/* end set_curve_impulse */

static void set_surf_impulse(
	ELASTIC_SET *geom_set,
	SURFACE *surf,
        GLOBAL_POINT **point_set)
{
	int j,k;
	TRI *tri;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	long gindex;

	unsort_surf_point(surf);
	hs = Hyper_surf(surf);
	for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
			tri = tri->next)
	{
	    hse = Hyper_surf_element(tri);
	    for (j = 0; j < 3; ++j)
	    {
            p = Point_of_tri(tri)[j];
            if (sorted(p) || Boundary_point(p)) continue;
            sorted(p) = YES;
            gindex = Gindex(p);
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
    
            for (k = 0; k < 3; ++k)
            {
                sl->impulse[k] = point_set[gindex]->impuls[k];
                sr->impulse[k] = point_set[gindex]->impuls[k];
            }
	    }
	}
}	/* end set_surf_impulse */

extern void set_equilibrium_mesh(
	Front *front)
{
	switch (front->rect_grid->dim)
	{
	case 3:
	    set_equilibrium_mesh3d(front);
	    return;
    default:
        printf("ERROR: dim must equal 3\n");
        LOC(); clean_up(EXIT_FAILURE);
	}
}	/* end set_equilibrium_mesh */


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
        //TODO: THIS MAKES NO SENSE
        if ((hsbdry_type(*c) != STRING_HSBDRY && hsbdry_type(*c) != DISKGAP_STRING_HSBDRY)
                || hsbdry_type(*c) != GORE_HSBDRY) continue;
        
        //printf("hsbdry_type(*c): %d\n",hsbdry_type(*c));

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
    
        //TODO: UNDEFINED IF DISKGAP_STRING_HSBDRY ...
        never_redistribute(Hyper_surf(*c)) = YES;
	}

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_BAND_BOUNDARY) continue;
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

    /*
	printf("Original length:\n");
	printf("min_len = %16.12f\n",min_len);
	printf("max_len = %16.12f\n",max_len);
	printf("ave_len = %16.12f\n",ave_len/count);
    */

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY && 
            wave_type(*s) != ELASTIC_BAND_BOUNDARY) continue;
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
	    if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_BAND_BOUNDARY) continue;
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
    /*
	printf("Perturbed length:\n");
	printf("min_len = %16.12f\n",min_len);
	printf("max_len = %16.12f\n",max_len);
	printf("ave_len = %16.12f\n",ave_len/count);
    */
}	/* end set_equilibrium_mesh3d */

