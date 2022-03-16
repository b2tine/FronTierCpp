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
#include<cmath>

using Tr= CGAL::Surface_mesh_default_triangulation_3; 
using C2T3 = CGAL::Complex_2_in_triangulation_3<Tr>;
using GT = Tr::Geom_traits;
using Sphere_3 = GT::Sphere_3;
using Point_3 = GT::Point_3;
using FT = GT::FT;

using Vertex_handle = C2T3::Vertex_handle;
using Vertex_iterator = C2T3::Vertex_iterator;
using Facet_iterator = C2T3::Facet_iterator;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using cgalPoint3 = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<cgalPoint3>;
using Vertex_index = Mesh::Vertex_index;
using Face_index = Mesh::Face_index;

//cgal_surf.cpp
extern void CGAL_MakeSphericalSurf(Front*,double*,double*,COMPONENT,COMPONENT,
                            int,int,SURFACE**);
extern void CGAL_MakeEllipsoidalSurf(Front*,double*,double*,COMPONENT,COMPONENT,
                            int,int,SURFACE**);
extern void CGAL_MakeCuboidSurf(Front*,double*,double*,COMPONENT,COMPONENT,int,
                            int,SURFACE**);

extern void CGAL_MakeConeSurf(Front*,double*,double,double,COMPONENT,COMPONENT,
                            int,int,SURFACE**);

extern void CGAL_MakeCapsuleSurf(Front*,double*,double,COMPONENT,COMPONENT,
                            int,int,SURFACE**);

extern void CGAL_MakeCylindricalSurf(Front*,double*,double,double,int,COMPONENT,
                            COMPONENT,int,int,SURFACE**);

extern void CGAL_MakeCylindricalShellSurf(Front*,double*,double,double,int,COMPONENT,
                            COMPONENT,int,int,SURFACE**);

//cgal_intfc.cpp
template <typename CGAL_Surface,
          typename CGAL_MeshCriteria,
          typename CGAL_ManifoldTag>
void CGAL_MakeLevelSurface(Front*,SURFACE**,COMPONENT,COMPONENT,int,
        CGAL_Surface*,CGAL_MeshCriteria*,CGAL_ManifoldTag*);

CGAL::Surface_mesh_default_criteria_3<Tr> CGAL_GenerateMeshCriteria(double,double);

void CGAL_C2T3_to_FronTier(Front*,SURFACE**,COMPONENT,COMPONENT,C2T3*);

//Remove this function if above CGAL_C2T3_to_FronTier() shown to work correctly
void GenerateCgalSurface_OLD(Front*,SURFACE**,COMPONENT,COMPONENT,C2T3*,double*);


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

struct cylinder_function
{
    double* center;
    double radius;
    double height;
    int idir;

    cylinder_function(double* cen, double rad, double h, int dir)
        : center{cen}, radius{rad}, height{h}, idir{dir}
    {}

    FT operator()(Point_3 p) const
    {
        const FT coords[3] = {p.x(), p.y(), p.z()};
   
        const FT x = coords[(idir+1)%3] - center[(idir+1)%3];
        const FT y = coords[(idir+2)%3] - center[(idir+2)%3];
        const FT z = coords[idir] - center[idir];

        if (z > -0.5*height && z < 0.5*height)
            return sqr(x) + sqr(y) - sqr(radius);
        else
            return 1.0;
    }    
};

struct cone_function
{
    double* center; //center of lower boundary disk
    double slope;
    double height;

    cone_function(double* cen, double sl, double h)
        : center{cen}, slope{sl}, height{h}
    {
        //for now assume positive slope
        assert(slope > 0.0);
    }

    FT operator()(Point_3 p) const
    {
        const FT x = p.x() - center[0];
        const FT y = p.y() - center[1];
        const FT z = p.z() - center[2];

        if (z > 0 && z < height)
            return sqr(slope)*(sqr(x) + sqr(y)) - sqr(z);
        else
            return 1.0;
    }
};

struct capsule_function
{
    double* nose;
    double radius;

    capsule_function(double* cen, double rad)
        : nose{cen}, radius{rad}
    {}

    double lower_cone_height = radius/std::tan(70.0*M_PI/180.0);
    double lower_slope = lower_cone_height/radius;
    
    //TODO: Add spherical nose section
    //
    //double nose_radius = 0.05;
    
    double mid_height = lower_cone_height*0.1374;

    double upper_cone_height = radius*std::tan(44.5*M_PI/180.0);
        //double upper_cone_height = radius*std::tan(47.0*M_PI/180.0);
    double upper_slope = upper_cone_height/radius;
    double upper_frustrum_height = lower_cone_height*1.672;

    double top_frustrum_height = lower_cone_height*0.25;
    double radius_top_frustrum = (upper_cone_height - upper_frustrum_height)*std::tan(44.5*M_PI/180.0);

    FT operator()(Point_3 p) const
    {
        const FT x = p.x() - nose[0];
        const FT y = p.y() - nose[1];
        const FT z = p.z() - nose[2];

        if (z > 0 && z <= lower_cone_height)
        {
            return std::sqrt(sqr(x) + sqr(y)) - (radius/lower_cone_height)*z;
        }
        else if (z > 0 && (z - lower_cone_height) <= mid_height)
        {
            return std::sqrt(sqr(x) + sqr(y)) - radius;
        }
        else if (z > 0 && (z - lower_cone_height - mid_height) <= upper_frustrum_height)
        {
            return std::sqrt(sqr(x) + sqr(y))
                - (radius - 1.0*(radius/upper_cone_height)*(z - lower_cone_height - mid_height));
        }
        else if (z > 0 && (z - lower_cone_height - mid_height - upper_frustrum_height) < top_frustrum_height)
        {
            return std::sqrt(sqr(x) + sqr(y))
                - (radius_top_frustrum - 
                        1.0*(radius_top_frustrum/top_frustrum_height)*(z - lower_cone_height - mid_height - upper_frustrum_height));
        }
        else
        {
            return 1.0;
        }
    }
};

#endif
