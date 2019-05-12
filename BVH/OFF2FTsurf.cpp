#include "BVH_util.h"
#include <ifluid_state.h>

static void initSTATEvelocity(Front*, double*);
static void dummySpringSolver(Front*);
static void elastic_point_propagate(Front*, POINTER, POINT*,
        POINT*, HYPER_SURF_ELEMENT*, HYPER_SURF*, double, double*);
static void mono_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
static void propagation_driver(Front*);


int main(int argc, char* argv[])
{
    static Front front;
    static RECT_GRID comp_grid;
    static F_BASIC_DATA f_basic;
    static LEVEL_FUNC_PACK level_func_pack;

    f_basic.dim = 3;
    FT_Init(argc,argv,&f_basic);

    char* in_name = f_basic.in_name;
    char* out_name = f_basic.out_name;

    Mesh inmesh;
    std::ifstream input(in_name);
    if (!input || !(input >> inmesh) || !CGAL::is_triangle_mesh(inmesh))
    {
        std::cerr << "Error: input file must be a triangular mesh OFF file\n";
        return 1;
    }

    std::ofstream outfile(std::string(out_name) + "/input-mesh.off");
    outfile << inmesh;
    outfile.close();

    auto Bounds = getInputMeshDimensionsWithPad(&inmesh,1.0);
    auto lb = Bounds.first;
    auto ub = Bounds.second;

    std::cout << "Domain Lower Bound: ";
    std::cout << lb[0] << " " << lb[1] << " " << lb[2] << "\n";
    std::cout << "Domain Upper Bound: ";
    std::cout << ub[0] << " " << ub[1] << " " << ub[2] << "\n";

    f_basic.L[0] = lb[0];   f_basic.L[1] = lb[1];    f_basic.L[2] = lb[2];
    f_basic.U[0] = ub[0];   f_basic.U[1] = ub[1];    f_basic.U[2] = ub[2];
    f_basic.gmax[0] = 50;   f_basic.gmax[1] = 50;   f_basic.gmax[2] = 40;

    f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
    f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
    f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
    
    f_basic.size_of_intfc_state = sizeof(STATE);

    FT_StartUp(&front,&f_basic);
    add_to_debug("trace");
    //add_to_debug("optimize_intfc");
    add_to_debug("BVH");

    level_func_pack.pos_component = 1;
    FT_InitIntfc(&front,&level_func_pack);
    
    TriMeshOFF2MonoCompSurf(&front,&inmesh);
    //optimizeElasticMesh(&front);
    static_mesh(front.interf) = YES;

    char dname[100];
    sprintf(dname,"%s/geomview-interface",out_name);
    gview_plot_interface(dname,front.interf);
    
    bool drawBVH = true;
    BVH bvh(&front,drawBVH);
    auto root_bv = bvh.getRoot()->getBV();
    root_bv.print();

    //velocity function parameters
    TRANS_PARAMS trans_params;
    trans_params.dim = 3;
    trans_params.vel[0] = 0.25;
    trans_params.vel[1] = 1.0/3.0;
    trans_params.vel[2] = -0.5;

    front.vparams = &trans_params;
    front.vfunc = translation_vel;

    front.curve_propagate = mono_curve_propagate;
    PointPropagationFunction(&front) = elastic_point_propagate;

    propagation_driver(&front);
    clean_up(0);
}

void propagation_driver(Front* front)
{
	front->max_time = 1.0; 
	front->max_step = 1000;
	front->print_time_interval = 0.01;
	front->movie_frame_interval = 0.01;

    double CFL = Time_step_factor(front) = 0.75;
    Tracking_algorithm(front) = STRUCTURE_TRACKING;

    FT_ResetTime(front);

    //Output the initial interface.
    FT_Save(front);
    FT_Draw(front);

    //Startup
    //TODO: This function computes the force and torque
    //      exerted on rigid bodies by the fluid.
    //FrontPreAdvance(front);
    
    FT_Propagate(front);
    FT_SetTimeStep(front);
    FT_SetOutputCounter(front);

	FT_TimeControlFilter(front);
	FT_PrintTimeStamp(front);

	for (;;)
    {
	    /* Propagating interface for time step dt */
	    FT_Propagate(front);
	    FT_AddTimeStepToCounter(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.
        FT_SetTimeStep(front);

	    /* Output section */
	    FT_TimeControlFilter(front);
	    FT_PrintTimeStamp(front);

        if (FT_IsSaveTime(front))
            FT_Save(front);

        if (FT_IsDrawTime(front))
            FT_Draw(front);

        if (FT_TimeLimitReached(front))
            break;
	}

}

void mono_curve_propagate(Front* front, POINTER wave,
        CURVE* oldc, CURVE* newc, double dt)
{
    POINT* oldp, *newp;
	BOND *oldb,*newb;
	double V[MAXD];
	BOND_TRI **btris;
	HYPER_SURF_ELEMENT *oldhse;
	HYPER_SURF         *oldhs;

	for (oldb = oldc->first, newb = newc->first;
            oldb != NULL; oldb = oldb->next, newb = newb->next)
	{
	    oldp = oldb->end;
	    newp = newb->end;
	    for (btris = Btris(oldb); btris && *btris; ++btris)
	    {
	    	oldp->hse = oldhse = Hyper_surf_element((*btris)->tri);
	    	oldp->hs = oldhs = Hyper_surf((*btris)->surface);
		    elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	    }
	}
}


void elastic_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
    fourth_order_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
    ft_assign(left_state(oldp),left_state(newp),front->sizest);
    ft_assign(right_state(oldp),right_state(newp),front->sizest);
    //ft_assign(left_state(newp),left_state(oldp),front->sizest);
    //ft_assign(right_state(newp),right_state(oldp),front->sizest);
}

//Forward Euler
void dummySpringSolver(Front* front)
{
	INTERFACE *intfc = front->interf;
	POINT *p;
	HYPER_SURF *hs;
    HYPER_SURF_ELEMENT *hse;
	
    (void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs)) 
    {
        if(Boundary_point(p))
            continue;
        
        STATE* sl = (STATE*)left_state(p);
        for (int i = 0; i < FT_Dimension(); ++i)
            Coords(p)[i] = sl->x_old[i] + front->dt*(sl->vel[i]);
    }

	CURVE**c;
	BOND* bond;
	intfc_curve_loop(intfc,c)
    {
        for ((bond) = (*c)->first; (bond) != (*c)->last; 
                (bond) = (bond)->next)
        {
            p = bond->end;
            STATE* sl = (STATE*)left_state(p);
            for (int i = 0; i < FT_Dimension(); ++i)	
                Coords(p)[i] = sl->x_old[i] + front->dt*(sl->vel[i]);
	    }
	}

	NODE** n;
	intfc_node_loop(intfc,n)
    {
        p = (*n)->posn;
        STATE* sl = (STATE*)left_state(p);
        for (int i = 0; i < FT_Dimension(); ++i)
            Coords(p)[i] = sl->x_old[i] + front->dt*(sl->vel[i]);
	}
}

void initSTATEvelocity(Front* front, double* vel)
{
    STATE* sl;
    STATE* sr;
    
    TRI* tri;
    POINT* p;

    SURFACE** s;
    INTERFACE* intfc = front->interf;
    intfc_surface_loop(intfc,s)
    {
        if( is_bdry(*s) ) continue;
        surf_tri_loop(*s,tri)
        {
            for( int i = 0; i < 3; ++i )
            {
                p = Point_of_tri(tri)[i];
                sl = (STATE*)left_state(p);
                sr = (STATE*)right_state(p);
                for( int j = 0; j < 3; ++j )
                {
                    sl->vel[j] = sr->vel[j] = vel[j];
                    sl->x_old[j] = Coords(p)[j];
                }
            }
        }
    }
}
