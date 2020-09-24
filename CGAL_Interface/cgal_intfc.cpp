#include "cgal_intfc.h"

extern void CGAL_MakeSphericalSurf(
        Front* front,
        double* center,
        double radius,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        int refinement_level,
        SURFACE** surf)
{
        double radii[3] = {radius,radius,radius};
        CGAL_MakeEllipsoidalSurf(front,center,radii,neg_comp,pos_comp,
                            w_type,refinement_level,surf);
}

//TODO: Max local feature size, max_lfs, can probably be tightened up
//      in the below cgal surface functions.
extern void CGAL_MakeEllipsoidalSurf(
        Front* front,
        double* center,
        double* radii,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        int refinement_level,
        SURFACE** surf)
{
        ellipsoid_function func(center,radii);
       
        double max_radius = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            if (radii[i] > max_radius)
                max_radius = radii[i];
        }

        FT bounding_sphere_squared_radius = sqr(1.0+max_radius);
        Point_3 bounding_sphere_center(center[0],center[1],center[2]);

        Sphere_3 bounding_sphere(bounding_sphere_center,
                                 bounding_sphere_squared_radius);

        double error_bound = 1.0e-06;//from the level surface
        CGAL::Implicit_surface_3<GT,ellipsoid_function>
            cgal_surface(func,bounding_sphere,error_bound);
        
        //TODO: Follow this pattern for notion of refinement_level
        //      in remain cgal meshing functions.
        double epsilon = 0.16;
        double max_lfs = max_radius;
        epsilon /= (int)refinement_level;

        CGAL::Surface_mesh_default_criteria_3<Tr>
            cgal_mesh_criteria = CGAL_GenerateMeshCriteria(epsilon,max_lfs);

        CGAL::Manifold_tag cgal_manifold_tag;

        CGAL_MakeLevelSurface<
            CGAL::Implicit_surface_3<GT,ellipsoid_function>,
            CGAL::Surface_mesh_default_criteria_3<Tr>,
            CGAL::Manifold_tag
                >(front,surf,neg_comp,pos_comp,w_type,
                        &cgal_surface,&cgal_mesh_criteria,&cgal_manifold_tag);
}

extern void CGAL_MakeCuboidSurf(
        Front* front,
        double* center,
        double* edges,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        int refinement_level,
        SURFACE **surf)
{
        cuboid_function func(center,edges);

        double max_edge_dist = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            if (edges[i] > max_edge_dist)
                max_edge_dist = edges[i];
        }
        
        FT bounding_sphere_squared_radius = sqr(1.0+max_edge_dist);
        Point_3 bounding_sphere_center(center[0],center[1],center[2]);

        Sphere_3 bounding_sphere(bounding_sphere_center,
                                 bounding_sphere_squared_radius);

        CGAL::Implicit_surface_3<GT,cuboid_function>
            cgal_surface(func,bounding_sphere,1.0e-06);

        //TODO: provide justification for this value of epsilon.
        double epsilon = 0.0425;
        double max_lfs = max_edge_dist;
        CGAL::Surface_mesh_default_criteria_3<Tr>
            cgal_mesh_criteria = CGAL_GenerateMeshCriteria(epsilon,max_lfs);

        CGAL::Manifold_tag cgal_manifold_tag;

        CGAL_MakeLevelSurface<
            CGAL::Implicit_surface_3<GT,cuboid_function>,
            CGAL::Surface_mesh_default_criteria_3<Tr>,
            CGAL::Manifold_tag
                >(front,surf,neg_comp,pos_comp,w_type,
                        &cgal_surface,&cgal_mesh_criteria,&cgal_manifold_tag);
}

extern void CGAL_MakeCylindricalSurf(
        Front* front,
        double* center,
        double radius,
        double height,
        int idir,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        int refinement_level,
        SURFACE** surf)
{
        cylinder_function func(center,radius,height,idir);
       
        double max_radius = std::max(radius,0.5*height);

        FT bounding_sphere_squared_radius = sqr(1.0+max_radius);
        Point_3 bounding_sphere_center(center[0],center[1],center[2]);

        Sphere_3 bounding_sphere(bounding_sphere_center,
                                 bounding_sphere_squared_radius);

        CGAL::Implicit_surface_3<GT,cylinder_function>
            cgal_surface(func,bounding_sphere,1.0e-06);
        
        //TODO: provide justification for this value of epsilon.
        double epsilon = 0.0425;
        double max_lfs = max_radius;
        CGAL::Surface_mesh_default_criteria_3<Tr>
            cgal_mesh_criteria = CGAL_GenerateMeshCriteria(epsilon,max_lfs);

        CGAL::Manifold_tag cgal_manifold_tag;

        CGAL_MakeLevelSurface<
            CGAL::Implicit_surface_3<GT,cylinder_function>,
            CGAL::Surface_mesh_default_criteria_3<Tr>,
            CGAL::Manifold_tag
                >(front,surf,neg_comp,pos_comp,w_type,
                        &cgal_surface,&cgal_mesh_criteria,&cgal_manifold_tag);
}

extern void CGAL_MakeConeSurf(
        Front* front,
        double* center,
        double slope,
        double height,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        int refinement_level,
        SURFACE** surf)
{
        cone_function func(center,slope,height);
        double radius = height/slope;
        double max_dist = std::max(radius,height);

        FT bounding_sphere_squared_radius = sqr(1.0+max_dist);
        Point_3 bounding_sphere_center(center[0],center[1],center[2]);

        Sphere_3 bounding_sphere(bounding_sphere_center,
                                 bounding_sphere_squared_radius);

        CGAL::Implicit_surface_3<GT,cone_function>
            cgal_surface(func,bounding_sphere,1.0e-06);
        
        //TODO: provide justification for this value of epsilon.
        double epsilon = 0.0425;
        double max_lfs = max_dist;
        CGAL::Surface_mesh_default_criteria_3<Tr>
            cgal_mesh_criteria = CGAL_GenerateMeshCriteria(epsilon,max_lfs);

        CGAL::Manifold_tag cgal_manifold_tag;

        CGAL_MakeLevelSurface<
            CGAL::Implicit_surface_3<GT,cone_function>,
            CGAL::Surface_mesh_default_criteria_3<Tr>,
            CGAL::Manifold_tag
                >(front,surf,neg_comp,pos_comp,w_type,
                        &cgal_surface,&cgal_mesh_criteria,&cgal_manifold_tag);

}

//TODO: Read about cgal meshing criteria and develop
//      a notion of refinement_level to adjust the values
//      of local feature size and epsilon.

//For details of mesh criteria see the following:
//
//  1. https://doc.cgal.org/latest/Surface_mesher/index.html
//
//  2. "Provably good sampling and meshing of surfaces"
//      Graphical Models, 67:405–451, 2005 
//      By Jean-Daniel Boissonnat and Steve Oudot.


CGAL::Surface_mesh_default_criteria_3<Tr>
CGAL_GenerateMeshCriteria(double epsilon, double max_lfs)
{
        //epsilon must be < 0.16 for theoretical gaurantees descibed
        //in the CGAL 3d surface mesh manual.

        //max_lfs is the max local feature size, i.e max distance
        //from points sampled on the level surface to the medial axis
        //of the surface, and should be approximated conservatively.
        
        //e.g. For a sphere, max_lfs should be set to the sphere radius.
        //     For an ellipsoid, max_lfs should be set to the largest of
        //     the 3 provided radii.
        
        //Lower bound on the minimum angle of surface mesh facets.
        //Must be <= 30.0 in order to guarantee convergence.
        double lb_ang = 30.0;

        //Upper bound on radius of surface Delauney balls.
        //Follows from the definition of (loose) epsilon-sampled surfaces.
        double ub_rad = epsilon*max_lfs;
        
        //Upper bound on the distance between centers of surface mesh
        //facets and centers of the Delauney Balls circumscribing them.
        double ub_dist = 4.5*sqr(epsilon)*max_lfs;

        CGAL::Surface_mesh_default_criteria_3<Tr> criteria(lb_ang,ub_rad,ub_dist);
        return criteria;
}       /* end CGAL_GenerateMeshCriteria */

template <typename CGAL_Surface,
          typename CGAL_MeshCritera,
          typename CGAL_ManifoldTag>

void CGAL_MakeLevelSurface(
        Front* front,
        SURFACE** surf,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        CGAL_Surface* cgal_surface,
        CGAL_MeshCritera* cgal_mesh_criteria,
        CGAL_ManifoldTag* cgal_manifold_tag)
{
        Tr tr;
        C2T3 c2t3(tr);

        CGAL::make_surface_mesh(c2t3,*cgal_surface,
                *cgal_mesh_criteria,*cgal_manifold_tag);

        CGAL_C2T3_to_FronTier(front,surf,neg_comp,pos_comp,&c2t3);
        wave_type(*surf) = w_type;
    
        if (consistent_interface(front->interf) == NO)
            clean_up(ERROR);

        front->interf->modified = YES;
        interface_reconstructed(front->interf) = NO;
}       /* end CGAL_MakeLevelSurface */

void CGAL_C2T3_to_FronTier(
        Front* front,
        SURFACE** surf,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        C2T3* c2t3)
{
        INTERFACE* intfc = front->interf;

        //save a copy of the current interface
        INTERFACE* saved_intfc = current_interface();

        //create a surface for the current interface
        set_current_interface(intfc);
        SURFACE* newsurf = make_surface(neg_comp,pos_comp,NULL,NULL);

        //output c2t3 facets to a CGAL Triangle Mesh
        Mesh mesh;
        CGAL::facets_in_complex_2_to_triangle_mesh<C2T3,Mesh>(*c2t3,mesh);

        //allocate memory for the mesh vertices
        POINT** points;
        const unsigned int num_vtx = mesh.number_of_vertices();
        FT_VectorMemoryAlloc((POINTER*)&points, num_vtx, sizeof(POINT*));

        //read in the mesh vertices
        std::map<Vertex_index,int> vmap;
        Mesh::Vertex_range vrange = mesh.vertices();
        Mesh::Vertex_range::iterator vit = vrange.begin();
        for( int i = 0; i < num_vtx; i++ )
        {
            double vcoords[3];
            cgalPoint3 p = mesh.point(*vit);
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
        const unsigned int num_faces = mesh.number_of_faces();
        FT_VectorMemoryAlloc((POINTER*)&tris, num_faces, sizeof(TRI*));

        //read in the mesh triangles
        Mesh::Face_range frange = mesh.faces();
        Mesh::Face_range::iterator fit = frange.begin();
        for( int i = 0; i < num_faces; i++ )
        { 
            CGAL::Vertex_around_face_circulator<Mesh> vb(mesh.halfedge(*fit),
                                    mesh); 
            CGAL::Vertex_around_face_circulator<Mesh> ve(vb);
            std::vector<Vertex_index> vidx;

            do {            
                vidx.push_back(*vb);
                vb++;
            } while( vb != ve );

            int i0 = vmap[vidx[0]];
            int i1 = vmap[vidx[1]];
            int i2 = vmap[vidx[2]];

            tris[i] = make_tri(points[i0],points[i1],points[i2],NULL,NULL,NULL,
                                NO);

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
        }//Note: This was the i = 0 iteration of the for loop below

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
                        if (Point_of_tri(tri1)[m] == Point_of_tri(tri2)[(l+1)%3]
                           &&
                           Point_of_tri(tri2)[l] == Point_of_tri(tri1)[(m+1)%3])
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
}       /* end CGAL_C2T3_to_FronTier */

//TODO: Save PosOrientPts() and GenerateCgalSurface_OLD() for now.
//      Can remove if CGAL_C2T3_to_FronTier() shown to work correctly.
static bool PosOrientPts(
        POINT *p1,
        POINT *p2,
        POINT *p3,
        double* bs_center)
{
        double v1[MAXD],v2[MAXD],v[MAXD];
        double norm[MAXD];
        for (int i = 0; i < 3; ++i)
        {
            v1[i] = Coords(p2)[i] - Coords(p1)[i];
            v2[i] = Coords(p3)[i] - Coords(p1)[i];
                
            double x0 = Coords(p1)[i] - bs_center[i];
            double x1 = Coords(p2)[i] - bs_center[i];
            double x2 = Coords(p3)[i] - bs_center[i];
            norm[i] = (x0 + x1 + x2)/3.0;

        }

        Cross3d(v1,v2,v);
        if (Dot3d(v,norm) > 0.0) return true;
        else return false;
}       /* end OrientOfPoints */

void GenerateCgalSurface_OLD(
        Front* front,
        SURFACE** surf,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        C2T3* c2t3,
        double* bs_center)
{
        INTERFACE* intfc = front->interf;
        INTERFACE* sav_intfc = current_interface();
        set_current_interface(intfc);

        SURFACE* newsurf = make_surface(neg_comp,pos_comp,NULL,NULL);

        std::map<Vertex_handle,int> V;

        int num_vtx = 0;
        Vertex_iterator vit;
        vit = c2t3->vertices_begin();

        //Create vertex map and find number of vertices
        while(vit != c2t3->vertices_end())
        {
            V[vit] = num_vtx++;
            vit++;
        }
        printf("NumMeshVertices = %d\n",num_vtx);

        POINT** points;
        uni_array(&points,num_vtx,sizeof(POINT*));
    
        int ii = 0;
        vit = c2t3->vertices_begin();

        while(vit != c2t3->vertices_end())
        {
            double pcoords[MAXD];

            pcoords[0] = (double) vit->point()[0];
            pcoords[1] = (double) vit->point()[1];
            pcoords[2] = (double) vit->point()[2];
        
            points[ii] = Point(pcoords);
            points[ii]->num_tris = 0;

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

        while (fit != c2t3->facets_end())
        {
            typename C2T3::Cell_handle cell = fit->first;
            int index = fit->second;

            int i1 = V[cell->vertex(tr1.vertex_triple_index(index,0))];
            int i2 = V[cell->vertex(tr1.vertex_triple_index(index,1))];
            int i3 = V[cell->vertex(tr1.vertex_triple_index(index,2))];

            if (PosOrientPts(points[i1],points[i2],points[i3],bs_center))
                tris[ii] = make_tri(points[i1],points[i2],points[i3],
                                               NULL,NULL,NULL,NO);
            else
                tris[ii] = make_tri(points[i1],points[i3],points[i2],
                                               NULL,NULL,NULL,NO);
        
            tris[ii]->surf = newsurf;

            points[i1]->num_tris++;
            points[i2]->num_tris++;
            points[i3]->num_tris++;

            ii++;
            fit++;
        }

        //From here this is the same as lines 1264-1319 of cgal.cpp
    
        int num_point_tris = 3*num_tris;
        intfc->point_tri_store = (TRI**) store(num_point_tris*sizeof(TRI*));
        TRI** ptris = intfc->point_tri_store;
        POINT* p;

        for (int i = 0; i < num_vtx; i++)
        {
            points[i]->tris = ptris;
            ptris += points[i]->num_tris;
            points[i]->num_tris = 0;
        }

        for (int i = 0; i < num_tris; i++)
        {
            if (i != 0)
            {
                tris[i]->prev = tris[i-1];
                tris[i-1]->next = tris[i];
            }

            for (int j = 0; j < 3; j++)
            {
                p = Point_of_tri(tris[i])[j];
                p->tris[p->num_tris++] = tris[i];
            }
        }

        for (int i = 0; i < num_vtx; i++)
        {
            ptris = points[i]->tris;
            int num_ptris = points[i]->num_tris;
            for (int j = 0; j < num_ptris; j++)
            {
                for (int k = 0; k < j; k++)
                {
                    TRI* tri1 = ptris[j];
                    TRI* tri2 = ptris[k];
                    for (int m = 0; m < 3; m++)
                    for (int l = 0; l < 3; l++)
                    {
                        if (Point_of_tri(tri1)[m] == Point_of_tri(tri2)[(l+1)%3] &&
                            Point_of_tri(tri1)[(m+1)%3] == Point_of_tri(tri2)[l])
                        {
                            Tri_on_side(tri1,m) = tri2;
                            Tri_on_side(tri2,l) = tri1;
                        }
                        else if (Point_of_tri(tri1)[m] == Point_of_tri(tri2)[l]
                                &&
                                Point_of_tri(tri1)[(m+1)%3] == Point_of_tri(tri2)[(l+1)%3])
                        {
                            printf("Orientation is opposite!\n");
                            clean_up(ERROR);
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
