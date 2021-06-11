#include "BVH_util.h"

//TODO: This function should also determine and return
//      the necessary number of grid blocks, gmax, in order
//      for the interface to propagate correctly. 
std::pair<std::vector<double>,std::vector<double> >
getInputMeshDimensionsWithPad(Mesh* mesh, double pad)
{
    std::vector<double> L(3,HUGE);
    std::vector<double> U(3,-HUGE);
    Mesh::Vertex_range vrange = mesh->vertices();
    Mesh::Vertex_range::iterator vit = vrange.begin();
    for( vit; vit != vrange.end(); ++vit )
    {
        double vcoords[3];
        cgalPoint3 p = mesh->point(*vit);
        vcoords[0] = p.x();
        vcoords[1] = p.y();
        vcoords[2] = p.z();
        for( int i = 0; i < 3; ++i )
        {
            if( vcoords[i] < L[i] )
                L[i] = vcoords[i];
            if( vcoords[i] > U[i] )
                U[i] = vcoords[i];            
        }
    }

    for( int i = 0; i < 3; ++i )
    {
        L[i] -= pad;
        U[i] += pad;
    }
    return std::make_pair(L,U);
}

void TriMeshOFF2MonoCompSurf(Front* front, Mesh* mesh)
{
    SURFACE* surf;
    INTERFACE* intfc = front->interf;
    COMPONENT amb_comp = intfc->default_comp;
    TriMeshOFF2Surf(intfc,amb_comp,amb_comp,mesh,&surf);
    wave_type(surf) = ELASTIC_BOUNDARY;
    FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);
    if (consistent_interface(front->interf) == NO)
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
    std::map<Vertex_index,int> vmap;
    Mesh::Vertex_range vrange = mesh->vertices();
    Mesh::Vertex_range::iterator vit = vrange.begin();
    for( int i = 0; i < num_vtx; i++ )
    {
        double vcoords[3];
        cgalPoint3 p = mesh->point(*vit);
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





