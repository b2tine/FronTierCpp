#include "bending.h"


static void calculateBendingForce3d2003(POINT* p1, TRI* tri, TRI* n_tri, double bends, double bendd);
static void calculateBendingForce3d2006(POINT* p1, TRI* tri, TRI* n_tri, double bends, double bendd);
static void calculateBendingForce3dparti(POINT* p1, TRI* tri, TRI* n_tr, double bends, double bendd);

static double calOriLeng(int index1, int index2, TRI* tri, TRI* n_tri);
static double divEx(double numerator, double denominator);
static void DebugShow(const double& sva);
    
static void computeCurvatureBinormal(BOND* b1, BOND* b2);
static void computeGradCurvatureBinormal(BOND* b1, BOND* b2);
static void computeStringPointBendingForce(BOND* b1, BOND* b2);

static std::vector<std::vector<double>> crossMat(double* u);
static std::vector<std::vector<double>> outerProduct(double* u, double* v, int dim);
static std::vector<std::vector<double>> scalarMultMatrix(
        double c, std::vector<std::vector<double>> M, int dimX, int dimY);
static std::vector<std::vector<double>> matAdd(
        std::vector<std::vector<double>> A, std::vector<std::vector<double>> B, int dimX, int dimY);
static std::vector<std::vector<double>> matMinus(
        std::vector<std::vector<double>> A, std::vector<std::vector<double>> B, int dimX, int dimY);
static std::vector<double> multMatVec(std::vector<std::vector<double>> A, std::vector<double> u, int m, int n);



//TODO: NEED TO CALL SOMEWHERE AFTER INTERFACE INITIALIZATION
void addStringBenders(Front* front)
{
    double string_bends = 0.0;

    FILE* infile = fopen(InName(front),"r");
    if (CursorAfterStringOpt(infile,"Enter yes for string bending force:"))
    {
        char string[100];
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] != 'y' || string[0] != 'Y')
        {
            CursorAfterString(infile,"Enter string bending stiffness constant:");
            fscanf(infile,"%lf",&string_bends);
        }
    }
    else
    {
        fclose(infile);
        return;
    }
    fclose(infile);
    

    INTERFACE* intfc = front->interf;
    CURVE** curve;
    BOND* bond;

    intfc_curve_loop(intfc,curve)
    {
        if (is_bdry(*curve)) continue;
        if (hsbdry_type(*curve) != STRING_HSBDRY) continue;

        //curve_bond_loop(*curve,bond)
        for (bond = (*curve)->first; bond != (*curve)->last; bond = bond->next)
        {
            BOND_BENDER* bond_bender = nullptr;
            FT_ScalarMemoryAlloc((POINTER*)&bond_bender,sizeof(BOND_BENDER));
            bond_bender->bends = string_bends;
            POINT* pt = bond->end;
            pt->extra = bond_bender;
        }
    }

    //TODO: Need to be able to restart with BOND_BENDER in pt->extra
}

void computeSurfBendingForce(INTERFACE* intfc, const double bends, const double bendd)
{
    //TODO: add switch to turn off (early retur)

    SURFACE **surf;
    TRI *tri;

    intfc_surface_loop(intfc, surf)
    {
        if (is_bdry(*surf)) continue;
        if (wave_type(*surf) != ELASTIC_BOUNDARY) continue;
        
        unsort_surf_point(*surf);//TODO: can remove?
        
        //NOTE: next surf_tri_loop() replaces the function
        //      clear_surf_point_force(*surf);
        surf_tri_loop(*surf, tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                POINT* p = Point_of_tri(tri)[i];
                STATE* sl = (STATE*)left_state(p);
                for (int j = 0; j < 3; ++j)
                    sl->bendforce[j] = 0.0;
            }
        }
        
        surf_tri_loop(*surf, tri)
        {
            if (is_bdry(*surf)) continue;
            if (wave_type(*surf) != ELASTIC_BOUNDARY) continue;
            
            //TODO: edge traversal of surface (requires half the computation)
            for (int i = 0; i < 3; ++i)
            {
                POINT *p = Point_of_tri(tri)[i];
                TRI* n_tri = Tri_on_side(tri, (i+1)%3);

                if (is_side_bdry(tri,(i+1)%3)) continue; 
                
                calculateBendingForce3d2006(p,tri,n_tri,bends,bendd);
                //calculateBendingForce3d2003(p,tri,n_tri,bends,bendd);//TODO: test this one
                    //calculateBendingForce3dparti(p,tri,n_tri,bends,bendd);
            }
        }
    }
}

void computeStringBendingForce(INTERFACE* intfc)
{
    CURVE** curve;
    BOND *b;

    intfc_curve_loop(intfc,curve)
    {
        if (is_bdry(*curve)) continue;
        if (hsbdry_type(*curve) != STRING_HSBDRY) continue;

        curve_bond_loop(*curve,b)
        {
            STATE* ss = (STATE*)left_state(b->start);
            STATE* se = (STATE*)left_state(b->end);
            for (int j = 0; j < 3; ++j)
            {
                ss->bendforce[j] = 0.0;
                se->bendforce[j] = 0.0;
            }
        }

        //TODO: Continue string bending stiffness implementation.
        //      Skeleton is below; need implementations for the
        //      following functions, and storage for intermediate
        //      computations.
        
        for (b = (*curve)->first; b != (*curve)->last; b = b->next)
        {
            //compute and store (kb)_i in x_i joining e_{i-1} and e_i
            computeCurvatureBinormal(b,b->next);
        }
        
        for (b = (*curve)->first; b != (*curve)->last; b = b->next)
        {
            //compute and store grad_{i}(kb)_j for j = i-1, i, i+1 at x_i
            computeGradCurvatureBinormal(b,b->next);
        }

        for (b = (*curve)->first; b != (*curve)->last; b = b->next)
        {
            //compute bending force on x_i
            computeStringPointBendingForce(b,b->next);
        }
    }

} /* setBendingForce3d */


//Appears to be from "Simple Linear Bending Stiffness in Particle Systems"
//Pascal Volino and Nadia Magnenat-Thalmann
void calculateBendingForce3d2006(
        POINT* p1,
        TRI* t1,
        TRI* t2,
        const double bends,
        const double bendd)
{
        if (Mag3d(Tri_normal(t1)) < MACH_EPS ||
            Mag3d(Tri_normal(t2)) < MACH_EPS) return;

        int index = Vertex_of_point(t1, p1);
        POINT *p3 = Point_of_tri(t1)[(index+1)%3];
        POINT *p4 = Point_of_tri(t1)[(index+2)%3];
        index = 3 - Vertex_of_point(t2, p3) - Vertex_of_point(t2, p4);
        POINT *p2 = Point_of_tri(t2)[index];

        double x13[3], x14[3], x23[3], x24[3], E[3];
        for (int i = 0; i < 3; ++i)
        {
            x13[i] = Coords(p1)[i] - Coords(p3)[i];
            x14[i] = Coords(p1)[i] - Coords(p4)[i];
            x23[i] = Coords(p2)[i] - Coords(p3)[i];
            x24[i] = Coords(p2)[i] - Coords(p4)[i];
            E[i] = Coords(p4)[i] - Coords(p3)[i];
        }
        double N1[3], N2[3], N3[3], N4[3];
        Cross3d(x13, x14, N1);
        Cross3d(x24, x23, N2);
        Cross3d(x23, x13, N3);
        Cross3d(x14, x24, N4);
        double N1_mag, N2_mag, N3_mag, N4_mag;
        N1_mag = Mag3d(N1);
        N2_mag = Mag3d(N2);
        N3_mag = Mag3d(N3);
        N4_mag = Mag3d(N4);
        double a1, a2, a3, a4;
        a1 = N2_mag / (N1_mag + N2_mag);
        a2 = 1 - a1;
        a3 = - N4_mag / (N3_mag + N4_mag);
        a4 = -1 - a3;
        double R[3];
        for (int i = 0; i < 3; ++i)
            R[i] = a1 * Coords(p1)[i] + a2 * Coords(p2)[i] +
                   a3 * Coords(p3)[i] + a4 * Coords(p4)[i];
        double E_mag, h1, h2;
        E_mag = Mag3d(E);
        h1 = fabs(Dot3d(x13, E) / E_mag);
        h1 = sqrt(Dot3d(x13, x13) - h1 * h1);
        h2 = fabs(Dot3d(x23, E) / E_mag);
        h2 = sqrt(Dot3d(x23, x23) - h2 * h2);

	    //double bend_stiff = getBendStiff(); 
        double bend_stiff = bends;
        double lambda = 2.0*(h1 + h2)*E_mag*bend_stiff / (3.0*h1*h2*h1*h2);

    STATE* state[4];
    POINT* pts[4] = {p1,p2,p3,p4};
    double a[4] = {a1,a2,a3,a4};

    /*
    for (int j = 0; j < 4; ++j)
    {
        if (state[j]->is_fixed || state[j]->is_registeredpt) continue;
        
        //TODO: Should force be redistributed to non fixed/registered point?
        //      See below for naive redistribution ...
        
        for (int i = 0; i < 3; ++i)
        {
            state[j]->bendforce[i] += -1.0*lambda*a[j]*R[i]*0.5;
        }
    }
    */

    for (int j = 0; j < 2; ++j)
    {
        int j0 = 2*j;
        int j1 = 2*j + 1;

        state[j0] = static_cast<STATE*>(left_state(pts[j0]));
        state[j1] =  static_cast<STATE*>(left_state(pts[j1]));

        if (!(state[j0]->is_fixed || state[j0]->is_registeredpt) &&
            !(state[j1]->is_fixed || state[j1]->is_registeredpt))
        {
            for (int i = 0; i < 3; ++i)
            {
                state[j0]->bendforce[i] += -1.0*lambda*a[j0]*R[i]*0.5;
                state[j1]->bendforce[i] += -1.0*lambda*a[j1]*R[i]*0.5;
            }
        }
        else if (!(state[j0]->is_fixed || state[j0]->is_registeredpt) &&
                 (state[j1]->is_fixed || state[j1]->is_registeredpt))
        {
            //TODO: correct???
            for (int i = 0; i < 3; ++i)
            {
                state[j0]->bendforce[i] += -1.0*lambda*a[j0]*R[i]*0.5;
                state[j0]->bendforce[i] += -1.0*lambda*a[j1]*R[i]*0.5;
            }
        }
        else if ((state[j0]->is_fixed || state[j0]->is_registeredpt) &&
                 !(state[j1]->is_fixed || state[j1]->is_registeredpt))
        {
            //TODO: correct???
            for (int i = 0; i < 3; ++i)
            {
                state[j1]->bendforce[i] += -1.0*lambda*a[j0]*R[i]*0.5;
                state[j1]->bendforce[i] += -1.0*lambda*a[j1]*R[i]*0.5;
            }
        }
        /*
        else
        {
            //Both points pinned; do nothing.
        }
        */    
    }
	
    /*
    for (int i = 0; i < 3; ++i)
    {
        //p1->force[i] += -lambda * a1 * R[i] * 0.5;
        //p2->force[i] += -lambda * a2 * R[i] * 0.5;
        //p3->force[i] += -lambda * a3 * R[i] * 0.5;
        //p4->force[i] += -lambda * a4 * R[i] * 0.5;
    }
	
    //TODO: use state->bendforce
    if (fabs(Coords(p1)[0] -0.75) < 0.1)
	{
	    double R_mag = Mag3d(R);
	    int count = 0;
	    if (Coords(p1)[2] > 0.8) count++;
	    if (Coords(p2)[2] > 0.8) count++;
	    if (Coords(p3)[2] > 0.8) count++;
	    if (Coords(p4)[2] > 0.8) count++;
	    if (R_mag > 0 && count == 3)
	    {
	    printf("R = [%f %f %f], R_mag = %e, lambda = %e\n", R[0],R[1],R[2],R_mag,lambda);
	    printf("nor R = [%f %f %f]\n", R[0]/R_mag, R[1]/R_mag, R[2]/R_mag);
	    printf("h1 = %f, h2 = %f, E_mag = %f\n", h1, h2, E_mag);
	    printf("p1 = [%f %f %f]\n", Coords(p1)[0], Coords(p1)[1], Coords(p1)[2]); 
	    printf("p2 = [%f %f %f]\n", Coords(p2)[0], Coords(p2)[1], Coords(p2)[2]); 
	    printf("p3 = [%f %f %f]\n", Coords(p3)[0], Coords(p3)[1], Coords(p3)[2]); 
	    printf("p4 = [%f %f %f]\n", Coords(p4)[0], Coords(p4)[1], Coords(p4)[2]); 
	    }
	}
        
    if (Mag3d(p1->force) == 0)
    {
        printf("f1 = [%f %f %f]\n", p1->force[0], p1->force[1], p1->force[2]);
        printf("f2 = [%f %f %f]\n", p2->force[0], p2->force[1], p2->force[2]);
        printf("f3 = [%f %f %f]\n", p3->force[0], p3->force[1], p3->force[2]);
        printf("f4 = [%f %f %f]\n", p4->force[0], p4->force[1], p4->force[2]);
        printf("a = [%f %f %f %f]\n\n", a1, a2, a3, a4);
        printf("h1 = %e, h2 = %e\n", h1, h2);
    }
    */
}       /* end calculateBendingForce3d */


// Appears to be from "Simulation of Clothing with Folds and Wrinkles"
// R. Bridson, S. Marino and R. Fedkiw
void calculateBendingForce3d2003(
        POINT* p1,
        TRI* t1,
        TRI* t2,
        const double bends,
        const double bendd)
{
        if (Mag3d(Tri_normal(t1)) < MACH_EPS ||
            Mag3d(Tri_normal(t2)) < MACH_EPS) return;

        int index = Vertex_of_point(t1, p1);
        POINT *p3 = Point_of_tri(t1)[(index+1)%3];
        POINT *p4 = Point_of_tri(t1)[(index+2)%3];
        index = 3 - Vertex_of_point(t2, p3) - Vertex_of_point(t2, p4);
        POINT *p2 = Point_of_tri(t2)[index];
        double x13[3], x14[3], x23[3], x24[3], E[3];

	std::transform(Coords(p1), Coords(p1) + 3, Coords(p3), 
		x13, std::minus<double>());
	std::transform(Coords(p1), Coords(p1) + 3, Coords(p4), 
		x14, std::minus<double>());
	std::transform(Coords(p2), Coords(p2) + 3, Coords(p3), 
		x23, std::minus<double>());
	std::transform(Coords(p2), Coords(p2) + 3, Coords(p4), 
		x24, std::minus<double>());
	std::transform(Coords(p4), Coords(p4) + 3, Coords(p3), 
		E, std::minus<double>());

        double N1[3], N2[3], E_mag, N1_mag, N2_mag;
        double n1[3], n2[3];

        Cross3d(x13, x14, N1);
        Cross3d(x24, x23, N2);
        E_mag = Mag3d(E);
        N1_mag = Mag3d(N1);
        N2_mag = Mag3d(N2);
        std::transform(E, E + 3, E, 
		std::bind2nd(std::divides<double>(), E_mag));
	std::transform(N1, N1 + 3, n1, 
		std::bind2nd(std::divides<double>(), N1_mag));
	std::transform(N2, N2 + 3, n2, 
		std::bind2nd(std::divides<double>(), N2_mag));
        std::transform(n1, n1 + 3, N1, 
		std::bind2nd(std::divides<double>(), N1_mag));
	std::transform(n2, n2 + 3, N2, 
		std::bind2nd(std::divides<double>(), N2_mag));

        double u1[3], u2[3], u3[3], u4[3];

	std::transform(N1, N1 + 3, u1,
                std::bind1st(std::multiplies<double>(), Dot3d(x14, E)));
	std::transform(N2, N2 + 3, u2,
                std::bind1st(std::multiplies<double>(), Dot3d(x24, E)));
	std::transform(u1, u1 + 3, u2, u3, std::plus<double>()); 
        std::transform(N1, N1 + 3, u1,
                std::bind1st(std::multiplies<double>(), -Dot3d(x13, E)));
        std::transform(N2, N2 + 3, u2,
                std::bind1st(std::multiplies<double>(), -Dot3d(x23, E)));
        std::transform(u1, u1 + 3, u2, u4, std::plus<double>());
	std::transform(N1, N1 + 3, u1, 
	        std::bind1st(std::multiplies<double>(), E_mag));
	std::transform(N2, N2 + 3, u2, 
	        std::bind1st(std::multiplies<double>(), E_mag));

	//double bend_stiff = getBendStiff();
    double bend_stiff = bends;
        double coeff = bend_stiff * sqr(E_mag) / (N1_mag + N2_mag);

        if (Dot3d(n1, n2) > 1.0 + 1.0e-10)
        {
	    std::cout << std::fixed << std::setprecision(14); 
	    std::cout << "t1 = "; 
	    std::for_each(Tri_normal(t1),Tri_normal(t1) + 3, DebugShow); 
	    std::cout << std::endl; 
	    std::cout << "n1 = ";
            std::for_each(n1,n1 + 3, DebugShow);
            std::cout << std::endl;
	    std::cout << "t2 = ";
            std::for_each(Tri_normal(t2),Tri_normal(t2) + 3, DebugShow);
            std::cout << std::endl;
	    std::cout << "n2 = ";
            std::for_each(n2,n2 + 3, DebugShow);
            std::cout << std::endl;
	    std::cout << "Dot3d(n1, n2) = "; 
            DebugShow(Dot3d(n1, n2)); 
            std::cout << std::endl; 
	    std::cout << "u1 = ";
            std::for_each(u1,u1 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "u2 = ";
            std::for_each(u2,u2 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "u3 = ";
            std::for_each(u3,u3 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "u4 = ";
            std::for_each(u4,u4 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "N1 = ";
            std::for_each(N1,N1 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "N2 = ";
            std::for_each(N2,N2 + 3, DebugShow);
            std::cout << std::endl;
            clean_up(0);
        }

        double sine_half_theta = sqrt(0.5 * std::max(0.0, 1.0 - Dot3d(n1, n2)));
        double tmp[3];

        Cross3d(n1, n2, tmp);
        if (Dot3d(tmp, E) < 0)
            sine_half_theta *= -1.0;
        coeff *= sine_half_theta;
	
	//double bend_damp = getBendDamp(); 
	double bend_damp = bendd;
        double dtheta = 0.0;

        dtheta = Dot3d(u1, p1->vel) + Dot3d(u2, p2->vel) +
                 Dot3d(u3, p3->vel) + Dot3d(u4, p4->vel);
        if (fabs(dtheta) < 1.0e-10) dtheta = 0.0;
// used for debugging        
//	double sum[3] = {0.0};

//	for (int i = 0; i < 3; i++) 
//	     sum[i] = u1[i] + u2[i] + u3[i] + u4[i]; 
//	int summ = Mag3d(sum); 
//	if (summ > 1.0e-6)
//	{
//	    std::cout << "u vector is too large\n"; 
//	    std::cout << "sum is: " << sum[0] << ' ' << sum[1] 
//		<< ' ' << sum[2];
//	    std::cout << "p1 is: " << Coords(p1)[0] << ' ' 
//		<< Coords(p1)[1] << ' ' << Coords(p1)[2] << std::endl;    
//	}
        

// used for debugging
//      if (dtheta > 0.0)
//        {
//            double sum[3];
//            for (int i = 0; i < 3; ++i)
//                sum[i] = u1[i] + u2[i] + u3[i] + u4[i];
//           printf("sum = %20.14f %20.14f %20.14f\n",sum[0],sum[1],sum[2]);
//            printf("p1->vel = %f %f %f\n",p1->vel[0],p1->vel[1],p1->vel[2]);
//            printf("dtheta = %20.14f\n", dtheta);
//            clean_up(0);
//        }     

    coeff += -1.0*bend_damp*E_mag*dtheta;

	STATE* state1 = static_cast<STATE*>(left_state(p1));
	STATE* state2 = static_cast<STATE*>(left_state(p2));
	STATE* state3 = static_cast<STATE*>(left_state(p3));
	STATE* state4 = static_cast<STATE*>(left_state(p4));

    // each tri_pair will be calculated twice
    
    
    if (!state1->is_fixed && !state1->is_registeredpt)
    {
        std::transform(u1, u1 + 3, tmp, std::bind1st(std::multiplies<double>(), coeff)); 
        std::transform(state1->bendforce, state1->bendforce + 3, tmp, state1->bendforce, std::plus<double>());
            //std::transform(p1->force, p1->force + 3, tmp, p1->force, std::plus<double>());
    }
	
    if (!state2->is_fixed && !state2->is_registeredpt)
    {
        std::transform(u2, u2 + 3, tmp, std::bind1st(std::multiplies<double>(), coeff));
        std::transform(state2->bendforce, state2->bendforce + 3, tmp, state2->bendforce, std::plus<double>());
            //std::transform(p2->force, p2->force + 3, tmp, p2->force, std::plus<double>());
    }
    
    if (!state3->is_fixed && !state3->is_registeredpt)
    {
        std::transform(u3, u3 + 3, tmp, std::bind1st(std::multiplies<double>(), coeff));
        std::transform(state3->bendforce, state3->bendforce + 3, tmp, state3->bendforce, std::plus<double>());
            //std::transform(p3->force, p3->force + 3, tmp, p3->force, std::plus<double>());
    }
	
    if (!state4->is_fixed && !state4->is_registeredpt)
    {
        std::transform(u4, u4 + 3, tmp, std::bind1st(std::multiplies<double>(), coeff));
        std::transform(state4->bendforce, state4->bendforce + 3, tmp, state4->bendforce, std::plus<double>());
            //std::transform(p4->force, p4->force + 3, tmp, p4->force, std::plus<double>());
    }

}       /* calculateBendingForce3d */


//Appears to be a simple naive bending model
void calculateBendingForce3dparti(
        POINT* p1,
        TRI* tri,
        TRI* n_tri,
        const double bends,
        const double bendd)
{
    int index1, index2; 
    std::set<POINT*> pointset; 

    for (int i = 0; i < 3; i++)
    {
        if (Point_of_tri(tri)[i] == p1)
            index1 = i;
        pointset.insert(Point_of_tri(tri)[i]); 
    }
 
	POINT* p2; 

    //find opp point in second tri
    for (int i = 0; i < 3; i++) 
    {
	     if (pointset.find(Point_of_tri(n_tri)[i]) == pointset.end())
         {
             p2 = Point_of_tri(n_tri)[i]; 
		    index2 = i; 
		    break; 
	     }
    }

    double length0 = calOriLeng(index1, index2, tri, n_tri); 
	double length = separation(p1, p2, 3);
	STATE* state = static_cast<STATE*>(left_state(p1));
	double *vel1 = state->vel; 
	state = static_cast<STATE*>(left_state(p2)); 
	double *vel2 = state->vel;  
/*
	std::cout << "p1: " << p1 << std::endl; 
	std::cout << "p2: " << p2 << std::endl; 
	std::cout << std::setprecision(16) << "length: " << length << ' ' << "length0: " << std::setprecision(16) << length0 << std::endl; 
*/
	double velr[3] = {0.0}; 
	double dir[3] = {0.0}; 
    double dot = 0;

	for (int i = 0; i < 3; i++) {
	     velr[i] = vel2[i] - vel1[i]; 
	     dir[i] = Coords(p2)[i] - Coords(p1)[i]; 
             dot += velr[i]*dir[i];
	} 
        
    if (length > 1.0e-10)
    {
	    for (int i = 0; i < 3; i++) 
	         dir[i] /= length; 
    }


	STATE* statep1 = static_cast<STATE*>(left_state(p1));
    for (int i = 0; i < 3; i++)
    {
         statep1->bendforce[i] += bends * (length - length0) * dir[i]; 
            //p1->force[i] += bends * (length - length0) * dir[i]; 
         
         // bending force excerts along virtual edge
         statep1->bendforce[i] += bendd * velr[i]; 
                //p1->force[i] -= bendd * dot * dir[i]; 
         //p1->force[i] += bendd * velr[i]; 
    }
}

double calOriLeng(int index1, int index2, TRI* tri, TRI* n_tri)
{
    double c = tri->side_length0[(index1+1)%3]; 
    double b1 = tri->side_length0[index1]; 
    double a1 = tri->side_length0[3-index1-(index1+1)%3]; 
    double a2, b2; 
    POINT* p1 = Point_of_tri(tri)[(index1+1)%3]; 
    int pInN; 

    for (int i = 0; i < 3; i++) 
	 if (Point_of_tri(n_tri)[i] == p1) {
	     pInN = i; 
	     break; 
	 }
    if (pInN < index2 || pInN > index2 && pInN == 2) {
        b2 = n_tri->side_length0[pInN]; 
        a2 = n_tri->side_length0[index2]; 
    }
    else {
	b2 = n_tri->side_length0[index2]; 
	a2 = n_tri->side_length0[3-index2-pInN]; 
    }

    double cangle1, cangle2;

    try {
        cangle1 = divEx(sqr(a1) + sqr(c) - sqr(b1), (2 * a1 * c));
    } catch (std::overflow_error e) {
        std::cout << e.what() << " -> ";
        return a2;
    }
    try {
        cangle2 = divEx(sqr(a2) + sqr(c) - sqr(b2), (2 * a2 * c));
    } catch (std::overflow_error e) {
        std::cout << e.what() << " -> ";
        return a1;
    }

    double sangle1 = sqrt(1 - std::min(sqr(cangle1), 1.0)); 
    double sangle2 = sqrt(1 - std::min(sqr(cangle2), 1.0));
    double cangle = cangle1 * cangle2 - sangle1 * sangle2; 
    
    return sqrt(sqr(a1) + sqr(a2) - 2 * cangle * a1 * a2); 
}

double divEx(double numerator, double denominator)
{
    if (fabs(denominator) < 1.0e-10)
        throw std::overflow_error("Divide by zero exception!\n");
    return numerator / denominator;
}

void DebugShow(const double & sva)
{
    std::cout << std::setw(20) << sva << " ";
}

void computeCurvatureBinormal(BOND* b1, BOND* b2)
{
    POINT* pt = b1->end;
    if (!pt->extra) return;
    BOND_BENDER* bond_bender = (BOND_BENDER*)pt->extra;

    //Compute the curvature binormal vector:
    //  (kb)_i = (2*e_{i-1} ^ e_i)/(|e_{i-1}||e_i| + < e_{i-1}, e_i >
    
    double cross12[MAXD];
    vector_product_on_bonds(b1,b2,3,cross12);//cross product

    double len1 = bond_length(b1);
    double len2 = bond_length(b2);
    double dot12 = scalar_product_on_bonds(b1,b2,3);
    double denom = len1*len2 + dot12;

    for (int i = 0; i < 3; ++i)
    {
        bond_bender->kb[i] = 2.0*cross12[i]/denom;
    }
}

void computeGradCurvatureBinormal(BOND* b1, BOND* b2)
{
    POINT* p2 = b1->end;
    if (!p2->extra) return;
    BOND_BENDER* bbender2 = (BOND_BENDER*)p2->extra;

    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
        bbender2->gradprev_kb[i][j] = 0.0;
        bbender2->grad_kb[i][j] = 0.0;
        bbender2->gradnext_kb[i][j] = 0.0;
    }

    //TODO: take fixed points into consideration ...

    //compute grad_{i}(kb[i-1])
    if (b1->prev != nullptr)
    {
        BOND* b0 = b1->prev;
        double len0 = bond_length(b0);
        double len1 = bond_length(b1);
        double dot01 = scalar_product_on_bonds(b0,b1,3);
        double denom01 = len0*len1 + dot01;

        double e0[3];
        for (int i = 0; i < 3; ++i)
        {
            e0[i] = Coords(b0->end)[i] - Coords(b0->start)[i];
        }
        auto e0cross_mat = crossMat(e0);
        
        POINT* p1 = b0->end;
        BOND_BENDER* bbender1 = (BOND_BENDER*)p1->extra;
        double* kb1 = bbender1->kb;
        auto op10_mat = outerProduct(kb1,e0,3);
        
        auto grad2_kb1_mat = scalarMultMatrix(2.0,e0cross_mat,3,3);
        grad2_kb1_mat = matMinus(grad2_kb1_mat,op10_mat,3,3);
        grad2_kb1_mat = scalarMultMatrix(1.0/denom01,grad2_kb1_mat,3,3);

        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            bbender2->gradprev_kb[i][j] = grad2_kb1_mat[i][j];
        }
    }

    //compute grad_{i}(kb[i+1])
    if (b2->next != nullptr)
    {
        BOND* b3 = b2->next;
        double len2 = bond_length(b2);
        double len3 = bond_length(b3);
        double dot23 = scalar_product_on_bonds(b2,b3,3);
        double denom23 = len2*len3 + dot23;

        double e3[3];
        for (int i = 0; i < 3; ++i)
        {
            e3[i] = Coords(b3->end)[i] - Coords(b3->start)[i];
        }
        auto e3cross_mat = crossMat(e3);
        
        POINT* p3 = b3->start;
        BOND_BENDER* bbender3 = (BOND_BENDER*)p3->extra;
        double* kb3 = bbender3->kb;
        auto op33_mat = outerProduct(kb3,e3,3);
        
        auto grad2_kb3_mat = scalarMultMatrix(2.0,e3cross_mat,3,3);
        grad2_kb3_mat = matAdd(grad2_kb3_mat,op33_mat,3,3);
        grad2_kb3_mat = scalarMultMatrix(1.0/denom23,grad2_kb3_mat,3,3);

        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            bbender2->gradnext_kb[i][j] = grad2_kb3_mat[i][j];
        }
    }

    //compute grad_{i}(kb[i])
    double len1 = bond_length(b1);
    double len2 = bond_length(b2);
    double dot12 = scalar_product_on_bonds(b1,b2,3);
    double denom12 = len1*len2 + dot12;

    double e1[3], e2[3];
    double diff21[3];
    for (int i = 0; i < 3; ++i)
    {
        e1[i] = Coords(b1->end)[i] - Coords(b1->start)[i];
        e2[i] = Coords(b2->end)[i] - Coords(b2->start)[i];
        diff21[i] = e2[i] - e1[i];
    }
    auto e1cross_mat = crossMat(e1);
    auto e2cross_mat = crossMat(e2);
        
    double* kb2 = bbender2->kb;
        //auto op21_mat = outerProduct(kb2,e1,3);
        //auto op22_mat = outerProduct(kb2,e2,3);
    auto op2diff21_mat = outerProduct(kb2,diff21,3);

    auto grad2_kb2_mat = matAdd(e1cross_mat,e2cross_mat,3,3);
    grad2_kb2_mat = scalarMultMatrix(2.0,grad2_kb2_mat,3,3);
    grad2_kb2_mat = matAdd(grad2_kb2_mat,op2diff21_mat,3,3);
    grad2_kb2_mat = scalarMultMatrix(-1.0/denom12,grad2_kb2_mat,3,3);

    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
        bbender2->grad_kb[i][j] = grad2_kb2_mat[i][j];
    }
}

void computeStringPointBendingForce(BOND* b1, BOND* b2)
{
    POINT* p2 = b1->end;
    if (!p2->extra) return;
    BOND_BENDER* bbender2 = (BOND_BENDER*)p2->extra;

    double bendforce[3] = {0.0};
    double bends = bbender2->bends;

    //contribution from i-1
    if (b1->prev != nullptr)
    {
        BOND* b0 = b1->prev;
        double restlen0 = bond_length0(b0);
        double restlen1 = bond_length0(b1);
        double l1_bar = restlen0 + restlen1;

        std::vector<std::vector<double>> gradprev_kb_T(3,std::vector<double>(3,0.0));
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            gradprev_kb_T[i][j] = bbender2->gradprev_kb[j][i];
        }

        POINT* p1 = b0->end;
        BOND_BENDER* bbender1 = (BOND_BENDER*)p1->extra;
        std::vector<double> kb1 = std::vector<double>(bbender1->kb, bbender1->kb + 3);
        
        auto bendforce1 = multMatVec(gradprev_kb_T,kb1,3,3);
        for (int i = 0; i < 3; ++i)
        {
            bendforce1[i] *= -2.0*bends/l1_bar;
            bendforce[i] += bendforce1[i];
        }
    }

    //contribution from i+1
    if (b2->next != nullptr)
    {
        BOND* b3 = b2->next;
        double restlen2 = bond_length0(b2);
        double restlen3 = bond_length0(b3);
        double l3_bar = restlen2 + restlen3;

        std::vector<std::vector<double>> gradnext_kb_T(3,std::vector<double>(3,0.0));
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            gradnext_kb_T[i][j] = bbender2->gradnext_kb[j][i];
        }

        POINT* p3 = b3->start;
        BOND_BENDER* bbender3 = (BOND_BENDER*)p3->extra;
        std::vector<double> kb3 = std::vector<double>(bbender3->kb, bbender3->kb + 3);
        
        auto bendforce3 = multMatVec(gradnext_kb_T,kb3,3,3);
        for (int i = 0; i < 3; ++i)
        {
            bendforce3[i] *= -2.0*bends/l3_bar;
            bendforce[i] += bendforce3[i];
        }
    }

    //contribution from i
    double restlen1 = bond_length0(b1);
    double restlen2 = bond_length0(b2);
    double l2_bar = restlen1 + restlen2;

    std::vector<std::vector<double>> grad_kb_T(3,std::vector<double>(3,0.0));
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
        grad_kb_T[i][j] = bbender2->grad_kb[j][i];
    }

    std::vector<double> kb2 = std::vector<double>(bbender2->kb, bbender2->kb + 3);
    
    auto bendforce2 = multMatVec(grad_kb_T,kb2,3,3);
    for (int i = 0; i < 3; ++i)
    {
        bendforce2[i] *= -2.0*bends/l2_bar;
        bendforce[i] += bendforce2[i];
    }

    STATE* sl2 = static_cast<STATE*>(left_state(p2));
    for (int i = 0; i < 3; ++i)
    {
        sl2->bendforce[i] = bendforce[i];
    }
}

std::vector<std::vector<double>> crossMat(double* u)
{
    std::vector<std::vector<double>> cross_mat(3,std::vector<double>(3,0.0));

    cross_mat[0][0] = 0.0;   cross_mat[0][1] = -u[2];   cross_mat[0][2] = u[1];
    cross_mat[1][0] = u[2];  cross_mat[1][1] = 0.0;     cross_mat[1][2] = -u[0];
    cross_mat[2][0] = -u[1]; cross_mat[2][1] = u[0];    cross_mat[2][2] = 0.0;

    return cross_mat;
}

std::vector<std::vector<double>> outerProduct(double* u, double* v, int dim)
{
    std::vector<std::vector<double>> op_mat(dim,std::vector<double>(dim,0.0));
    
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        op_mat[i][j] = u[i]*v[j];
    }
    
    return op_mat;
}

std::vector<std::vector<double>> scalarMultMatrix(
        double c,
        std::vector<std::vector<double>> M,
        int dimX,
        int dimY)
{
    std::vector<std::vector<double>> cM(dimX,std::vector<double>(dimY,0.0));

    for (int i = 0; i < dimX; ++i)
    for (int j = 0; j < dimY; ++j)
    {
        cM[i][j] = c*M[i][j];
    }

    return cM;
}

std::vector<std::vector<double>> matAdd(
        std::vector<std::vector<double>> A,
        std::vector<std::vector<double>> B,
        int dimX,
        int dimY)
{
    std::vector<std::vector<double>> sum_mat(dimX,std::vector<double>(dimY,0.0));
    
    for (int i = 0; i < dimX; ++i)
    for (int j = 0; j < dimY; ++j)
    {
        sum_mat[i][j] = A[i][j] + B[i][j];
    }

    return sum_mat;
}

std::vector<std::vector<double>> matMinus(
        std::vector<std::vector<double>> A,
        std::vector<std::vector<double>> B,
        int dimX,
        int dimY)
{
    std::vector<std::vector<double>> diff_mat(dimX,std::vector<double>(dimY,0.0));
    
    for (int i = 0; i < dimX; ++i)
    for (int j = 0; j < dimY; ++j)
    {
        diff_mat[i][j] = A[i][j] - B[i][j];
    }

    return diff_mat;
}

std::vector<double> multMatVec(
        std::vector<std::vector<double>> A,
        std::vector<double> u,
        int m,
        int n)
{
    std::vector<double> vec(m,0.0);
    
    for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
    {
        vec[i] = A[i][j]*u[j];
    }

    return vec;
}

