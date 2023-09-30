#ifndef CGAL_SURF_H
#define CGAL_SURF_H

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

        epsilon /= (double)refinement_level;
        max_lfs *= (double)refinement_level;

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

extern void CGAL_MakeCapsuleSurf(
        Front* front,
        double* nose,
        double radius,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        int refinement_level,
        SURFACE** surf)
{
        capsule_function func(nose,radius);
        
        double height = radius/std::tan(70*M_PI*180.0);
        double max_dist = std::max(radius,height);

        FT bounding_sphere_squared_radius = sqr(1.0+max_dist);
        Point_3 bounding_sphere_nose(nose[0],nose[1],nose[2]);

        Sphere_3 bounding_sphere(bounding_sphere_nose,
                                 bounding_sphere_squared_radius);

        CGAL::Implicit_surface_3<GT,capsule_function>
            cgal_surface(func,bounding_sphere,1.0e-06);
        
        //TODO: provide justification for this value of epsilon.
        double epsilon = 0.0425;
        double max_lfs = max_dist;
        
        epsilon /= (double)refinement_level;
            //max_lfs *= (double)refinement_level;

        CGAL::Surface_mesh_default_criteria_3<Tr>
            cgal_mesh_criteria = CGAL_GenerateMeshCriteria(epsilon,max_lfs);

        CGAL::Manifold_tag cgal_manifold_tag;

        CGAL_MakeLevelSurface<
            CGAL::Implicit_surface_3<GT,capsule_function>,
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
        
        epsilon /= (double)refinement_level;
            //max_lfs *= (double)refinement_level;

        
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

extern void CGAL_MakeCylindricalShellSurf(
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
    double height_cylinder = 2.0*height;
    CGAL_MakeCylindricalSurf(front,center,radius,height_cylinder,idir,
            neg_comp,pos_comp,w_type,refinement_level,surf);

    interface_reconstructed(front->interf) = YES;

    PLANE_PARAMS plane_params;
    plane_params.N[0] = plane_params.N[1] = plane_params.N[2] = 0.0;
    plane_params.P[0] = center[0];
    plane_params.P[1] = center[1];
    plane_params.P[2] = center[2];

    plane_params.N[idir] = -1.0;
    plane_params.P[idir] = center[idir] + 0.5*height;
    FT_CutSurfBdry(*surf,plane_constr_func,(POINTER)&plane_params,NULL,0,0);

    plane_params.N[idir] = 1.0;
    plane_params.P[idir] = center[idir] - 0.5*height;
    FT_CutSurfBdry(*surf,plane_constr_func,(POINTER)&plane_params,NULL,0,0);

    interface_reconstructed(front->interf) = NO;
}

#endif
