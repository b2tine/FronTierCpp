#include "cgal_intfc.h"

#include<stdio.h>
#include<iostream>
#include<fstream>


char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	char dname[100];

	f_basic.dim = 3;
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = -1.0;	f_basic.L[1] = -1.0; 	f_basic.L[2] = -1.0;
	f_basic.U[0] = 1.0; 	f_basic.U[1] = 1.0; 	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 32;	f_basic.gmax[1] = 32;   f_basic.gmax[2] = 32;
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


	level_func_pack.neg_component = 1;
	level_func_pack.pos_component = 2;

    int neg_comp = 1;
    int pos_comp = 2;
    int w_type = NEUMANN_BOUNDARY;

    //CGAL Ellipsoidal Level Surface
    
    FT_StartUp(&front,&f_basic);
	FT_InitIntfc(&front,&level_func_pack);
    
    SURFACE* ellipsoid_surf;
    double ellipsoid_center[3] = {0.25, -0.1, 0.1};
    double ellipsoid_radii[3] = {0.25, 0.5, 0.25};

    CGAL_MakeEllipsoidalSurf(&front,ellipsoid_center,ellipsoid_radii,
            neg_comp,pos_comp,w_type,1,&ellipsoid_surf);
    
	sprintf(dname,"%s/CGAL_ellipsoid_intfc",out_name);
	gview_plot_interface(dname,front.interf);

    delete_interface(front.interf);

    
    //CGAL Cuboid Level Surface

    FT_StartUp(&front,&f_basic);
	FT_InitIntfc(&front,&level_func_pack);
    
    SURFACE* cuboid_surf;
    double cuboid_center[3] = {0.0, 0.1, -0.1};
    double cuboid_edges[3] = {0.5, 0.5, 0.25};

    CGAL_MakeCuboidSurf(&front,cuboid_center,cuboid_edges,
            neg_comp,pos_comp,w_type,1,&cuboid_surf);

	sprintf(dname,"%s/CGAL_cuboid_intfc",out_name);
	gview_plot_interface(dname,front.interf);

    delete_interface(front.interf);

    
    //CGAL Cylindrical Level Surface
    
    FT_StartUp(&front,&f_basic);
	FT_InitIntfc(&front,&level_func_pack);
    
    SURFACE* cylindrical_surf;
    double cylinder_center[3] = {0.0, 0.0, 0.0};
    double cylinder_radius = 0.25;
    double cylinder_height = 0.5;
    int cylinder_dir = 2;

    CGAL_MakeCylindricalSurf(&front,cylinder_center,cylinder_radius,
            cylinder_height,cylinder_dir,neg_comp,pos_comp,w_type,
            1,&cylindrical_surf);
    
	sprintf(dname,"%s/CGAL_cylindrical_intfc",out_name);
	gview_plot_interface(dname,front.interf);

    delete_interface(front.interf);

    
    //FronTier Ellipsoidal Level Surface

    FT_StartUp(&front,&f_basic);
	FT_InitIntfc(&front,&level_func_pack);
    
    SURFACE* ft_ellipsoid_surf;
    double* ft_ellipsoid_center = ellipsoid_center;
    double* ft_ellipsoid_radii = ellipsoid_radii;
        //double ft_ellipsoid_center[3] = {0.25, -0.1, 0.1};
        //double ft_ellipsoid_radii[3] = {0.25, 0.5, 0.25};

    FT_MakeEllipticSurf(&front,ft_ellipsoid_center,ft_ellipsoid_radii,
            neg_comp,pos_comp,w_type,1,&ft_ellipsoid_surf);

    sprintf(dname,"%s/FT_ellipsoid_intfc",out_name);
	gview_plot_interface(dname,front.interf);

    delete_interface(front.interf);
   

    //FronTier Cuboid Level Surface

    FT_StartUp(&front,&f_basic);
	FT_InitIntfc(&front,&level_func_pack);
    
    SURFACE* ft_cuboid_surf;
    double* ft_cuboid_center = cuboid_center;
    double* ft_cuboid_edges = cuboid_edges; 
        //double ft_cuboid_center[3] = {0.0, 0.1, -0.1};
        //double ft_cuboid_edges[3] = {0.5, 0.5, 0.25};

    FT_MakeCuboidSurf(&front,ft_cuboid_center,ft_cuboid_edges,
            neg_comp,pos_comp,w_type,1,&ft_cuboid_surf);

    sprintf(dname,"%s/FT_cuboid_intfc",out_name);
	gview_plot_interface(dname,front.interf);

    delete_interface(front.interf);

    
    //FronTier Cylindrical Level Surface
   
    FT_StartUp(&front,&f_basic);
	FT_InitIntfc(&front,&level_func_pack);
    
    SURFACE* ft_cylindrical_surf;
    double* ft_cylinder_center = cylinder_center;
    double ft_cylinder_radius = cylinder_radius;
    double ft_cylinder_height = cylinder_height;
    int ft_cylinder_dir = cylinder_dir;
        //double ft_cylinder_center[3] = {0.0, 0.0, 0.0};
        //double ft_cylinder_radius = 0.25;
        //double ft_cylinder_height = 0.5;
        //int ft_cylinder_dir = 2;

    FT_MakeCylinderSurf(&front,ft_cylinder_center,ft_cylinder_radius,
            0.5*ft_cylinder_height,ft_cylinder_dir,neg_comp,pos_comp,
            w_type,&ft_cylindrical_surf);
    
	sprintf(dname,"%s/FT_cylindrical_intfc",out_name);
	gview_plot_interface(dname,front.interf);

    delete_interface(front.interf);
    
	clean_up(0);
}
