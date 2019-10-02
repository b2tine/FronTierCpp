#ifndef CGAL_INTFC_H
#define CGAL_INTFC_H

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Origin.h>

#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <FronTier.h>

#include<cstdlib>

using Tr= CGAL::Surface_mesh_default_triangulation_3; 
using C2T3 = CGAL::Complex_2_in_triangulation_3<Tr>;
using GT = Tr::Geom_traits;
using Sphere_3 = GT::Sphere_3;
using Point_3 = GT::Point_3;
using FT = GT::FT;

using Vertex_handle =  C2T3::Vertex_handle;
using Vertex_iterator = C2T3::Vertex_iterator;
using Facet_iterator = C2T3::Facet_iterator;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using cgalPoint3 = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<cgalPoint3>;
using Vertex_index = Mesh::Vertex_index;
using Face_index = Mesh::Face_index;

extern void CGAL_MakeSphericalSurf(Front*,double*,double*,COMPONENT,COMPONENT,
                            int,int,SURFACE**);
extern void CGAL_MakeEllipsoidalSurf(Front*,double*,double*,COMPONENT,COMPONENT,
                            int,int,SURFACE**);
extern void CGAL_MakeCuboidSurf(Front*,double*,double*,COMPONENT,COMPONENT,int,
                            int,SURFACE**);

template <typename CGAL_Surface,
          typename CGAL_MeshCriteria,
          typename CGAL_ManifoldTag>
void CGAL_MakeLevelSurface(Front*,SURFACE**,COMPONENT,COMPONENT,int,
        CGAL_Surface*,CGAL_MeshCriteria*,CGAL_ManifoldTag*);

CGAL::Surface_mesh_default_criteria_3<Tr> CGAL_GenerateMeshCriteria(double,double);

void CGAL_C2T3_to_FronTier(Front*,SURFACE**,COMPONENT,COMPONENT,C2T3*);

//Remove this function if above CGAL_C2T3_to_FronTier() shown to work correctly
void GenerateCgalSurface_OLD(Front*,SURFACE**,COMPONENT,COMPONENT,C2T3*,double*);


//TODO: any downside to moving these to cgal_intfc.cpp?
struct ellipsoid_function
{
    double* center;
    double* radii;

    ellipsoid_function(double* cen, double* rad)
        : center{cen}, radii{rad}
    {}

    FT operator()(Point_3 p) const
    {
        const FT x = p.x() - center[0];
        const FT y = p.y() - center[1];
        const FT z = p.z() - center[2];
    
        return sqr(x)/sqr(radii[0]) +
               sqr(y)/sqr(radii[1]) +
               sqr(z)/sqr(radii[2]) - 1.0;
    }    
};

struct cuboid_function
{
    double* center;
    double* edges;

    cuboid_function(double* cen, double* edg)
        : center{cen}, edges{edg}
    {}

    FT operator()(Point_3 p) const
    {
        const FT x = p.x() - center[0];
        const FT y = p.y() - center[1];
        const FT z = p.z() - center[2];

        const FT X[3] = {x,y,z};

        FT dist = -HUGE;
        for (int i = 0; i < 3; ++i)
        {
            const FT i_dist = std::abs(X[i]) - edges[i];
            if (dist < i_dist) dist = i_dist;
        }
        
        return dist;
    }
};

#endif
