#include "BVH_util.h"
#include <ifluid_state.h>

extern void optimizeElasticMesh(Front*);

static void initSTATEvelocity(Front*, double*);
static void dummySpringSolver(Front*);
static void collision_point_propagate(Front*, POINTER, POINT*,
                                    POINT*, HYPER_SURF_ELEMENT*,
                                    HYPER_SURF*, double, double*);
static void collision_curve_propagate(Front*,POINTER,CURVE*,CURVE*,double);
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
    add_to_debug("optimize_intfc");
    add_to_debug("BVH");

    level_func_pack.pos_component = 1;
    FT_InitIntfc(&front,&level_func_pack);
    
    TriMeshOFF2MonoCompSurf(&front,&inmesh);

    //FT_RedistMesh(&front);
    optimizeElasticMesh(&front);
    static_mesh(front.interf) = YES;

    /*
    //velocity function parameters
    TRANS_PARAMS trans_params;
    trans_params.dim = 3;
    trans_params.vel[0] = 0.25;
    trans_params.vel[1] = 1.0/3.0;
    trans_params.vel[2] = -0.5;

    front.vparams = &trans_params;
    front.vfunc = translation_vel;
    */

    double initvel[3] = {0.0, 0.0, -2.5};
    initSTATEvelocity(&front, initvel);
    front.vfunc = NULL;

    front.curve_propagate = collision_curve_propagate;
    PointPropagationFunction(&front) = collision_point_propagate;

    propagation_driver(&front);
    clean_up(0);
}


void propagation_driver(Front *front)
{
	front->max_time = 1.0; 
	front->max_step = 2;
	front->print_time_interval = 0.01;
	front->movie_frame_interval = 0.01;
    double CFL = 0.75;

	//CollisionSolver *collision_solver = new CollisionSolver3d();

    Time_step_factor(front) = CFL;
	Tracking_algorithm(front) = STRUCTURE_TRACKING;

	Frequency_of_redistribution(front,GENERAL_WAVE) = 100000;
    printf("CFL = %f\n",CFL);
    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
            Frequency_of_redistribution(front,GENERAL_WAVE));

    //FT_RedistMesh(front);
    FT_ResetTime(front);

    // Always output the initial interface.
    FT_Save(front);
    FT_Draw(front);
    
    // This is a virtual propagation to get maximum front 
    // speed to determine the first time step.
    FT_Propagate(front);
    FT_SetTimeStep(front);
    FT_SetOutputCounter(front);

    FT_TimeControlFilter(front);
    FT_PrintTimeStamp(front);

    BVH bvh(front);
    //exit(0);

    for (;;)
    {
        /* Propagating interface for time step dt */
    
        if(debugging("CLOCK"))
            reset_clock();

        dummySpringSolver(front);
        FT_Propagate(front);
        
        //TODO: How can we update the hypersurf elements in the heirarchy
        //      when FT_Propagate calls this function?
        //          assign_interface_and_free_front(front,newfront);
        //      The Hse pointers in the leaf nodes become invalid after this call.
        

        //collision detect and handling
        /*
        collision_solver->assembleFromInterface(front->interf,front->dt);
        collision_solver->setFrictionConstant(0.0);
        collision_solver->resolveCollision();
        */

        FT_AddTimeStepToCounter(front);
        //bvh.drawUnlock();
        bvh.updateHeirarchy();
        exit(0);


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

	//delete collision_solver;

}       /* end propagation_driver */


void collision_curve_propagate(
	Front* front,
	POINTER wave,
	CURVE* oldc,
	CURVE* newc,
	double dt)
{
	int dim = 3;

    POINT* oldp = oldc->start->posn;
    POINT* newp = newc->start->posn;

    ft_assign(left_state(newp),left_state(oldp),front->sizest);
    ft_assign(right_state(newp),right_state(oldp),front->sizest);

    STATE* newsl = (STATE*)left_state(newp);
	STATE* oldsl = (STATE*)left_state(oldp);
	
    for (int i = 0; i < dim; ++i)
    {
        newsl->vel[i] = oldsl->vel[i];
        Coords(newp)[i] = Coords(oldp)[i] + dt*oldsl->vel[i];
        newsl->collsnImpulse[i] = 0.0;
        newsl->x_old[i] = Coords(oldp)[i];
    }

	oldp = oldc->end->posn;
    newp = newc->end->posn;
    newsl = (STATE*)left_state(newp);
	oldsl = (STATE*)left_state(oldp);
    ft_assign(left_state(newp),left_state(oldp),front->sizest);
    ft_assign(right_state(newp),right_state(oldp),front->sizest);

    for (int i = 0; i < dim; ++i)
    {
        newsl->vel[i] = oldsl->vel[i];
        Coords(newp)[i] = Coords(oldp)[i] + dt*oldsl->vel[i];
        newsl->collsnImpulse[i] = 0.0;
        newsl->x_old[i] = Coords(oldp)[i];
    }

	BOND* oldb;
    BOND* newb;

	for (oldb = oldc->first, newb = newc->first; oldb != oldc->last;
                oldb = oldb->next, newb = newb->next)
    {
        oldp = oldb->end;
        newp = newb->end;
        newsl = (STATE*)left_state(newp);
        oldsl = (STATE*)left_state(oldp);
        ft_assign(left_state(newp),left_state(oldp),front->sizest);
        ft_assign(right_state(newp),right_state(oldp),front->sizest);
    
        for (int i = 0; i < dim; ++i)
        {
            newsl->vel[i] = oldsl->vel[i];
            Coords(newp)[i] = Coords(oldp)[i] + dt*oldsl->vel[i];
            newsl->collsnImpulse[i] = 0.0;
            newsl->x_old[i] = Coords(oldp)[i];
        }
    }
}


void collision_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{

	STATE* sl = (STATE*)left_state(oldp);
    STATE* newsl = (STATE*)left_state(newp);
    ft_assign(left_state(newp),left_state(oldp),front->sizest);

    double vel[MAXD],s;
    int dim = front->rect_grid->dim;

	for (int i = 0; i < dim; ++i)
        vel[i] = sl->vel[i];

	for (int i = 0; i < dim; ++i)
	{
	    newsl->vel[i] = vel[i];
        Coords(newp)[i] = Coords(oldp)[i];
	    newsl->collsnImpulse[i] = 0.0;
	    newsl->x_old[i] = Coords(oldp)[i];
	}

	newsl->collsn_num = 0;
    s = mag_vector(V,dim);
    set_max_front_speed(dim,s,NULL,Coords(newp),front);
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
        int tid = 0;
        if( is_bdry(*s) ) continue;
        surf_tri_loop(*s,tri)
        {
            Tri_index(tri) = tid;
            //std::cout << "tri_index: " << tid << "\n";
            for( int i = 0; i < 3; ++i )
            {
                p = Point_of_tri(tri)[i];
                sl = (STATE*)left_state(p);
                sr = (STATE*)right_state(p);
                for( int j = 0; j < 3; ++j )
                {
                    sl->x_old[j] = Coords(p)[j];
                        
                    sl->vel[j] = 0.0;
                    if( tid >= 250 && tid <= 300 )
                        sl->vel[j] = vel[j];
                }
            }
            tid++;
        }
    }
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
