#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include "FronTier.h"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <limits>
#include <unordered_map>


using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point3>;
using Vertex_index = Mesh::Vertex_index;
using Face_index = Mesh::Face_index;

void TriMeshOFF2MonoCompSurf(Front*, Mesh*);
void TriMeshOFF2Surf(INTERFACE*, COMPONENT,
        COMPONENT, Mesh*, SURFACE**);

char *in_name, *out_name;

int main(int argc, char* argv[])
{
    static Front front;
    static RECT_GRID comp_grid;
    static F_BASIC_DATA f_basic;
    static LEVEL_FUNC_PACK level_func_pack;

    f_basic.dim = 3;
    FT_Init(argc,argv,&f_basic);

    //TODO: Set domain based on input mesh max and min coords.
    f_basic.L[0] = -4.0;    f_basic.L[1] = -4.0;    f_basic.L[2] = -1.0;
    f_basic.U[0] = 4.0;     f_basic.U[1] = 4.0;     f_basic.U[2] = 4.0;
    f_basic.gmax[0] = 32;   f_basic.gmax[1] = 32;   f_basic.gmax[2] = 32;
    f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
    f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
    f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
    
    f_basic.size_of_intfc_state = 0;

    in_name = f_basic.in_name;
    out_name = f_basic.out_name;
    FT_StartUp(&front,&f_basic);
   
    Mesh mesh;
    std::ifstream input(in_name);
    if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Error: input file must be a triangular mesh OFF file" << std::endl;
        return 1;
    }

    //TODO: Add boolean to LEVEL_FUNC_PACK and an execution branch
    //      to FT_InitIntfc() that allows the interface to be read in
    //      from a surface mesh OFF file. Would eliminate the need to
    //      set the LEVEL_FUNC_PACK::pos_component value erroneously.

    level_func_pack.pos_component = 1;
    FT_InitIntfc(&front,&level_func_pack);
    
    TriMeshOFF2MonoCompSurf(&front,&mesh);
    
    char dname[100];
    sprintf(dname,"%s/off_read",out_name);
    gview_plot_interface(dname,front.interf);
    print_interface(front.interf);

    clean_up(0);

}


void TriMeshOFF2MonoCompSurf(Front* front, Mesh* mesh)
{
    SURFACE* surf;
    INTERFACE* intfc = front->interf;
    COMPONENT amb_comp = intfc->default_comp;
    TriMeshOFF2Surf(intfc,amb_comp,amb_comp,mesh,&surf);
    wave_type(surf) = ELASTIC_BOUNDARY;
    FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);
    if( consistent_interface(front->interf) == NO )
        clean_up(ERROR);
}

//TODO: Rewrite as a CGALSurfaceMesh2FTSurf function, and make a seperate
//      function that reads the OFF file directly into FronTier.
void TriMeshOFF2Surf(INTERFACE* intfc, COMPONENT pos_comp,
        COMPONENT neg_comp, Mesh* mesh, SURFACE** surf)
{
    //save a copy of the current interface
    INTERFACE* saved_intfc = current_interface();

    //create a surface for the current interface
    set_current_interface(intfc);
    SURFACE* newsurf = make_surface(pos_comp,neg_comp,NULL,NULL);

    //allocate memory for the mesh vertices
    POINT** points;
    const unsigned int num_vtx = mesh->number_of_vertices();
    FT_VectorMemoryAlloc((POINTER*)&points, num_vtx, sizeof(POINT*));

    //read in the mesh vertices
    std::unordered_map<Vertex_index,int> vmap;
    Mesh::Vertex_range vrange = mesh->vertices();
    Mesh::Vertex_range::iterator vit = vrange.begin();
    for( int i = 0; i < num_vtx; i++ )
    {
        double vcoords[3];
        Point3 p = mesh->point(*vit);
        vcoords[0] = p.x();
        vcoords[1] = p.y();
        vcoords[2] = p.z();

        points[i] = Point(vcoords);
        points[i]->num_tris = 0;
        vmap[*vit] = i;
        vit++;
    }
    
    //allocate memory for the mesh triangles
    TRI** tris;
    const unsigned int num_faces = mesh->number_of_faces();
    FT_VectorMemoryAlloc((POINTER*)&tris, num_faces, sizeof(TRI*));

    //read in the mesh triangles
    Mesh::Face_range frange = mesh->faces();
    Mesh::Face_range::iterator fit = frange.begin();
    for( int i = 0; i < num_faces; i++ )
    { 
        CGAL::Vertex_around_face_circulator<Mesh> vb(mesh->halfedge(*fit), *mesh); 
        CGAL::Vertex_around_face_circulator<Mesh> ve(vb);
        std::vector<Vertex_index> vidx;

        do {            
            vidx.push_back(*vb);
            vb++;
        } while( vb != ve );

        int i0 = vmap[vidx[0]];
        int i1 = vmap[vidx[1]];
        int i2 = vmap[vidx[2]];

        tris[i] = make_tri(points[i0], points[i1], points[i2],
                NULL, NULL, NULL, NO);

        tris[i]->surf = newsurf;
        points[i0]->num_tris++;
        points[i1]->num_tris++;
        points[i2]->num_tris++;

        fit++;
    }
    
    //allocate storage in current interface table
    intfc->point_tri_store = (TRI**)store(3*num_faces*sizeof(TRI*));

    //distribute the storage to the points
    TRI** ptris = intfc->point_tri_store;
    for( int i = 0; i < num_vtx; i++ )
    {
        points[i]->tris = ptris;
        ptris += points[i]->num_tris;
        points[i]->num_tris = 0;
    }

    //build the DCEL corresponding to the surface mesh
    POINT* p;
    for( int j = 0; j < 3; j++ )
    {
        p = Point_of_tri(tris[0])[j];
        p->tris[p->num_tris++] = tris[0];
    }//Note: i = 0 in the for loop below

    for( int i = 1; i < num_faces; i++ )
    {
        tris[i]->prev = tris[i-1];
        tris[i-1]->next = tris[i];
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
                for( int l = 0; l < 3; l++ )
                {
                    if( Point_of_tri(tri1)[m] == Point_of_tri(tri2)[(l+1)%3] &&
                        Point_of_tri(tri2)[l] == Point_of_tri(tri1)[(m+1)%3] )
                    {
                        Tri_on_side(tri1,m) = tri2;
                        Tri_on_side(tri2,l) = tri1;
                    }
                }
            }
        }
    }

    //build the tri list for the surface
    newsurf->num_tri = num_faces;
    first_tri(newsurf) = tris[0];
    last_tri(newsurf) = tris[num_faces-1];
    last_tri(newsurf)->next = tail_of_tri_list(newsurf);
    first_tri(newsurf)->prev = head_of_tri_list(newsurf);
    reset_intfc_num_points(newsurf->interface);

    *surf = newsurf;
    set_current_interface(saved_intfc);
    FT_FreeThese(2,tris,points);
}


