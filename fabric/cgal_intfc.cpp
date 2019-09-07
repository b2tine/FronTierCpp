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


/*
*				cgal_compare.cpp:
*
*				Author: Brandon Ballentine
*
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/


#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Origin.h>

#include <FronTier.h> //This needs to be below CGAL includes for some reason

#include<stdio.h>
#include<iostream>
#include<fstream>

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT,Function> Surface_3;

typedef C2t3::Vertex_handle Vertex_handle;
typedef C2t3::Vertex_iterator Vertex_iterator; 
typedef C2t3::Facet_iterator Facet_iterator; 


FT sphere_function(Point_3 p)
{
    const FT x = p.x();
    const FT y = p.y();
    const FT z = p.z();
    return sqrt(x*x + y*y + z*z) - 1.0;
}


void Cgal_InitIntfc(Front*,Function);
void GenerateCgalSurface(Front*, SURFACE**, C2t3*);
//double LinfIntfcPoints(INTERFACE*);
static void computeError_surface(INTERFACE*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;


/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

double sphere_func(POINTER,double*);

typedef struct {
    double center[3];
	double radius;
} TEST_SPHERE_PARAMS;


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	TEST_SPHERE_PARAMS s_params;
	static LEVEL_FUNC_PACK level_func_pack;
	char dname[100];

	f_basic.dim = 3;
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = -1.0;	f_basic.L[1] = -1.0; 	f_basic.L[2] = -1.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; 	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 32;	f_basic.gmax[1] = 32; f_basic.gmax[2] = 32;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;


        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
                        right_flush(RestartStep,7));
        if (pp_numnodes() > 1)
            sprintf(restart_name,"%s-nd%s",restart_name, 
                                right_flush(pp_mynode(),4));

	FT_StartUp(&front,&f_basic);
    front.rect_grid->h;

	if (!RestartRun)
	{
	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
        
	    FT_InitIntfc(&front,&level_func_pack);
        Cgal_InitIntfc(&front,sphere_function);
	}

    /*
    FILE* efile;
    sprintf(ename,"%s/errorfile.txt",out_name);
    efile = fopen(ename,"w");
    fprintf(efile,"CGAL_Linfty = %24.16g\n\n",LinfIntfcPoints(front.interf));
    */

    computeError_surface(front.interf);

	sprintf(dname,"%s/cgal_intfc",out_name);
	gview_plot_interface(dname,front.interf);


    delete_interface(front.interf);

     
	FT_StartUp(&front,&f_basic);

	if (!RestartRun)
	{
        s_params.center[0] = 0.0;
        s_params.center[1] = 0.0;
        s_params.center[2] = 0.0;
        s_params.radius = 1.0;
    
	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
        level_func_pack.func_params = (POINTER) &s_params;
        level_func_pack.func = sphere_func; 
        level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
        
	    FT_InitIntfc(&front,&level_func_pack);
	}


    computeError_surface(front.interf);
//    fprintf(efile,"FT_Linfty = %24.16g\n\n",LinfIntfcPoints(front.interf));

	sprintf(dname,"%s/FT_intfc",out_name);
	gview_plot_interface(dname,front.interf);


	clean_up(0);

}


/********************************************************************
 *                 	Function Definitions                            *
 ********************************************************************/


double sphere_func(POINTER func_params, double *coords)
{
    TEST_SPHERE_PARAMS* s_params = (TEST_SPHERE_PARAMS*) func_params;
	double x0,y0,z0,R;
	double distance;

    x0 = s_params->center[0];
    y0 = s_params->center[1];
    z0 = s_params->center[2];
	R = s_params->radius;

	distance = sqrt(sqr(coords[0] - x0) +
               sqr(coords[1] - y0) + sqr(coords[2] - z0)) - R;

    return distance;

}       /* end sphere_func */


void Cgal_InitIntfc(Front* front, Function func)
{
    Tr tr;
    C2t3 c2t3(tr);
    SURFACE* surf;

    double epsilon = 0.0425;
    //double epsilon = 0.03;
    double lfs = 1.0;
    double ub_rad = epsilon*lfs;
    double ub_dist = 4.5*sqr(epsilon)*lfs;
    //double ub_dist = sqr(epsilon)*lfs;
    double at = 30.0;

    Surface_3 surface(func, Sphere_3(CGAL::ORIGIN,2.0));
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(at,ub_rad,ub_dist);
    CGAL::make_surface_mesh(c2t3,surface,criteria,
                                CGAL::Non_manifold_tag());

    GenerateCgalSurface(front,&surf,&c2t3);
    wave_type(surf) = FIRST_PHYSICS_WAVE_TYPE;
    FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);
    
    if( consistent_interface(front->interf) == NO )
        clean_up(ERROR);
}


void GenerateCgalSurface(Front* front, SURFACE** surf, C2t3* c2t3)
{
    COMPONENT amb_comp = front->interf->default_comp;
    COMPONENT neg_comp, pos_comp;
    INTERFACE* intfc;
    SURFACE* newsurf;
    INTERFACE* sav_intfc;

    /*
    neg_comp = amb_comp; //not true for closed interface e.g. sphere
    pos_comp = amb_comp;
    */
    neg_comp = amb_comp-1;
    pos_comp = amb_comp;

    intfc = front->interf;
    sav_intfc = current_interface();
    set_current_interface(intfc);

    newsurf = make_surface(neg_comp,pos_comp,NULL,NULL);

    std::map<Vertex_handle,int> V;

    int num_vtx = 0;
    Vertex_iterator vit;
    vit = c2t3->vertices_begin();

    //Create vertex map and find number of vertices
    while( vit != c2t3->vertices_end() )
    {
        V[vit] = num_vtx++;
        vit++;
    }
    printf("NumMeshVertices = %d\n",num_vtx);

    POINT** points;
    uni_array(&points,num_vtx,sizeof(POINT*));
    
    int ii = 0;
    vit = c2t3->vertices_begin();

    //printf("\t\t\t\t Vertex Coordinates \t\t\t\t Magnitude \n\n");
    while( vit != c2t3->vertices_end() )
    {
        double* pcoords;
        uni_array(&pcoords,3,sizeof(double)); 

        pcoords[0] = (double) vit->point()[0];
        pcoords[1] = (double) vit->point()[1];
        pcoords[2] = (double) vit->point()[2];
        
        /*
        double r = sqrt( sqr(pcoords[0]) + sqr(pcoords[1]) + sqr(pcoords[2]) );
        printf("%12.3g \t %12.3g \t %12.3g \t %12.5g\n",pcoords[0],
                                            pcoords[1],pcoords[2],r);
    Surface_3 surface(cgal_lf_pack->Cgal_func, Sphere_3(CGAL::ORIGIN,2.0));
        */

        points[ii] = Point(pcoords);
        points[ii]->num_tris = 0;

        FT_FreeThese(1,pcoords);

        ii++;
        vit++;

    }


    ii = 0;
    int num_tris = c2t3->number_of_facets();
    printf("NumMeshFacets = %d\n",num_tris);

    TRI** tris;
    uni_array(&tris,num_tris,sizeof(TRI*));
    Facet_iterator fit;
    fit = c2t3->facets_begin();
    Tr tr1 = c2t3->triangulation();

    while( fit != c2t3->facets_end() )
    {
        typename C2t3::Cell_handle cell = fit->first;
        int index = fit->second;

        int i1 = V[cell->vertex(tr1.vertex_triple_index(index,0))];
        int i2 = V[cell->vertex(tr1.vertex_triple_index(index,1))];
        int i3 = V[cell->vertex(tr1.vertex_triple_index(index,2))];

        tris[ii] = make_tri(points[i1],points[i2],points[i3],
                                               NULL,NULL,NULL,NO);
        
        tris[ii]->surf = newsurf;

        points[i1]->num_tris++;
        points[i2]->num_tris++;
        points[i3]->num_tris++;

        ii++;
        fit++;
    }


    //From here this is the same as lines 1264-1319 of ../airfoil/cgal.cpp
    
    int num_point_tris = 3*num_tris;
    intfc->point_tri_store = (TRI**) store(num_point_tris*sizeof(TRI*));
    TRI** ptris = intfc->point_tri_store;
    POINT* p;

    for( int i = 0; i < num_vtx; i++ )
    {
        points[i]->tris = ptris;
        ptris += points[i]->num_tris;
        points[i]->num_tris = 0;
    }

    for( int i = 0; i < num_tris; i++ )
    {
        if( i != 0 )
        {
            tris[i]->prev = tris[i-1];
            tris[i-1]->next = tris[i];
        }

        for( int j = 0; j < 3; j++ )
        {
            p = Point_of_tri(tris[i])[j];
            p->tris[p->num_tris++] = tris[i];
        }
    }

    for( int i = 0; i < num_vtx; i++ )
    {
        ptris = points[i]->tris;
        int num_ptris = points[i]->num_tris;
        for( int j = 0; j < num_ptris; j++ )
        {
            for( int k = 0; k < j; k++ )
            {
                TRI* tri1 = ptris[j];
                TRI* tri2 = ptris[k];
                for( int m = 0; m < 3; m++ )
                {
                    for( int l = 0; l < 3; l++ )
                    {
                        if( Point_of_tri(tri1)[m] == Point_of_tri(tri2)[(l+1)%3] &&
                                Point_of_tri(tri1)[(m+1)%3] == Point_of_tri(tri2)[l] )
                        {
                            Tri_on_side(tri1,m) = tri2;
                            Tri_on_side(tri2,l) = tri1;
                        }
                    }
                }
            }
        }
    }

    newsurf->num_tri = num_tris;

    first_tri(newsurf) = tris[0];
    last_tri(newsurf) = tris[num_tris-1];
    last_tri(newsurf)->next = tail_of_tri_list(newsurf);
    first_tri(newsurf)->prev = head_of_tri_list(newsurf);
    reset_intfc_num_points(newsurf->interface);

    *surf = newsurf;
    set_current_interface(sav_intfc);
}
// end GenerateCgalSurface


//static void computeError_surface(INTERFACE* intfc)
void computeError_surface(INTERFACE* intfc)
{
    int dim = intfc->dim;

    double x0 = 0.0; // Get these from level func params
    double y0 = 0.0;
    double z0 = 0.0;
    double R = 1.0;

    SURFACE** s;
    int num_surfaces = I_NumOfSurfaces(intfc);
    uni_array(&s,num_surfaces,sizeof(SURFACE*));
    I_ArrayOfSurfaces(intfc,s);

    POINT** points;
    int num_points;
    int num_tris;

    for( int i = 0; i < num_surfaces; i++ )
    {
        if( wave_type(s[i]) < FIRST_PHYSICS_WAVE_TYPE )
            continue;

        num_points = I_NumOfSurfPoints(s[i]);
        num_tris = I_NumOfSurfTris(s[i]);
        printf("\nnumber of points: %d\n",num_points);
        printf("number of triangles: %d\n",num_tris);
        uni_array(&points,num_points,sizeof(POINT*));
        I_ArrayOfSurfPoints(s[i],points);
        break;
    }

    double L1_error = 0.0;
    double L2_error = 0.0;
    double Li_error = 0.0;

    printf("\n\n\t\t\t\t Vertex Coordinates \t\t\t\t Magnitude \n\n");
    for( int i = 0; i < num_points; i++ )
    {

        double xi = static_cast<double>(Coords(points[i])[0]);
        double yi = static_cast<double>(Coords(points[i])[1]);
        double zi = static_cast<double>(Coords(points[i])[2]);

        double r = sqrt( sqr(xi-x0) + sqr(yi-y0) + sqr(zi-z0) );
        double delta = fabs(r-R);

        printf("%12.6g \t %12.6g \t %12.6g \t %24.16g\n",xi,yi,zi,r);

        L1_error += delta;
        L2_error += delta;

        if( delta > Li_error )
            Li_error = delta;
    }

    L1_error /= (double) num_points;
    L2_error /= (double) num_points;

    FT_FreeThese(2,s,points);

    char ename[100];
    sprintf(ename,"%s/errorfile.txt",out_name);
    static FILE* efile;
    if( efile == NULL )
        efile = fopen(ename,"w");

    fprintf(efile,"L1:%24.18g \t L2:%24.18g \t Linfty:%24.18g\n\n",
                                            L1_error,L2_error,Li_error);

    fflush(efile);

}







