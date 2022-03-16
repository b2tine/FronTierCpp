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

	f_basic.L[0] = -2.0;	f_basic.L[1] = -2.0; 	f_basic.L[2] = -2.0;
	f_basic.U[0] = 2.0; 	f_basic.U[1] = 2.0; 	f_basic.U[2] = 2.0;
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

    //CGAL Cone Level Surface
    
    FT_StartUp(&front,&f_basic);
	FT_InitIntfc(&front,&level_func_pack);
    
    SURFACE* capsule_surf;
    double capsule_nose[3] = {0.0, 0.0, -0.5};
    double capsule_radius = 0.5;

    int refinement_level = 2;

    CGAL_MakeCapsuleSurf(&front,capsule_nose,capsule_radius,
            neg_comp,pos_comp,w_type,refinement_level,&capsule_surf);
    
	//sprintf(dname,"%s/CGAL_capsule_intfc",out_name);
	//gview_plot_interface(dname,front.interf);

    FT_Draw(&front);

    delete_interface(front.interf);

	clean_up(0);
}
