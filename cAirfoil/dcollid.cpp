#include <FronTier.h>
#include "collid.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cfenv>

#include <omp.h>

//union find functions
inline POINT*& root(POINT*);
inline POINT*& tail(POINT*);

//default parameters for collision detection
double CollisionSolver3d::s_dt = DT;
double CollisionSolver3d::s_cr = 1.0;
bool CollisionSolver3d::s_detImpZone = false;

double CollisionSolver3d::s_eps = EPS;
double CollisionSolver3d::s_thickness = 0.001;
double CollisionSolver3d::s_k = 5000;
double CollisionSolver3d::s_m = 0.001;
double CollisionSolver3d::s_mu = 0.5;

double CollisionSolver3d::l_eps = 10.0*EPS;
double CollisionSolver3d::l_thickness = 0.005;
double CollisionSolver3d::l_k = 50000;
double CollisionSolver3d::l_m = 0.002;
double CollisionSolver3d::l_mu = 0.5;

//debugging variables
int CollisionSolver3d::moving_edg_to_edg = 0;
int CollisionSolver3d::moving_pt_to_tri = 0;
int CollisionSolver3d::is_coplanar = 0;
int CollisionSolver3d::edg_to_edg = 0;
int CollisionSolver3d::pt_to_tri = 0;

std::vector<double> CollisionSolver3d::CollisionTimes;


void CollisionSolver3d::setTimeStepSize(double new_dt){s_dt = new_dt;}
double CollisionSolver3d::getTimeStepSize(){return s_dt;}

void CollisionSolver3d::setFabricRoundingTolerance(double neweps){s_eps = neweps;}
double CollisionSolver3d::getFabricRoundingTolerance(){return s_eps;}

void CollisionSolver3d::setStringRoundingTolerance(double neweps){l_eps = neweps;}
double CollisionSolver3d::getStringRoundingTolerance(){return l_eps;}


//set restitution coefficient between rigid bodies
void   CollisionSolver3d::setRestitutionCoef(double new_cr){s_cr = new_cr;}
double CollisionSolver3d::getRestitutionCoef(){return s_cr;}

//fabric points
void CollisionSolver3d::setFabricThickness(double h){s_thickness = h;}
double CollisionSolver3d::getFabricThickness(){return s_thickness;}

void   CollisionSolver3d::setFabricSpringConstant(double new_k){s_k = new_k;}
double CollisionSolver3d::getFabricSpringConstant(){return s_k;}

void   CollisionSolver3d::setFabricFrictionConstant(double new_mu){s_mu = new_mu;}
double CollisionSolver3d::getFabricFrictionConstant(){return s_mu;}

void   CollisionSolver3d::setFabricPointMass(double new_m){s_m = new_m;}
double CollisionSolver3d::getFabricPointMass(){return s_m;}

//string points
void CollisionSolver3d::setStringThickness(double h){l_thickness = h;}
double CollisionSolver3d::getStringThickness(){return l_thickness;}

void   CollisionSolver3d::setStringSpringConstant(double new_k){l_k = new_k;}
double CollisionSolver3d::getStringSpringConstant(){return l_k;}

void   CollisionSolver3d::setStringFrictionConstant(double new_mu){l_mu = new_mu;}
double CollisionSolver3d::getStringFrictionConstant(){return l_mu;}

void   CollisionSolver3d::setStringPointMass(double new_m){l_m = new_m;}
double CollisionSolver3d::getStringPointMass(){return l_m;}

void CollisionSolver3d::setStrainLimit(double slim) {strain_limit = slim;}
void CollisionSolver3d::setStrainRateLimit(double srlim) {strainrate_limit = srlim;}

double CollisionSolver3d::setVolumeDiff(double vd){vol_diff = vd;}


CollisionSolver3d::~CollisionSolver3d()
{
    abt_proximity.reset();
    abt_collision.reset();
    clearHseList();
}

void CollisionSolver3d::clearHseList()
{
	for (unsigned i = 0; i < hseList.size(); ++i)
    {
		delete hseList[i];
	}
	hseList.clear();
}

//NOTE: Must be called before calling the spring solver
void CollisionSolver3d::assembleFromInterface(
	const INTERFACE* intfc, const double dt)
{
	setTimeStepSize(dt);
	clearHseList();

	SURFACE** s;
	CURVE** c;
	TRI *tri;
	BOND *b;

	int n_tri = 0;
    int n_bond = 0;
	
    //TODO: Collect each CD_HSE_TYPE in seperate hseLists?
    
    intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    unsort_surface_point(*s);
	    
        surf_tri_loop(*s,tri)
	    {
            CD_HSE_TYPE tag;
            if (wave_type(*s) == MOVABLE_BODY_BOUNDARY || 
                wave_type(*s) == NEUMANN_BOUNDARY)
            {
                tag = CD_HSE_TYPE::RIGID_TRI;
            }
            else 
            {
                tag = CD_HSE_TYPE::FABRIC_TRI;
            }
            
            hseList.push_back(new CD_TRI(tri,tag));
		    n_tri++;
	    }
	}

	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY)
            continue; 

        CD_HSE_TYPE tag = CD_HSE_TYPE::STRING_BOND;
	    curve_bond_loop(*c,b)
	    {
            hseList.push_back(new CD_BOND(b,tag));
		    n_bond++;
	    }
	}

    setSizeCollisionTimes(hseList.size());

	makeSet(hseList);
	createImpZoneForRG(intfc);
	setDomainBoundary(intfc->table->rect_grid.L, intfc->table->rect_grid.U);

	if (debugging("intfc_assembly")){
	    printf("%d num of tris, %d num of bonds\n",n_tri,n_bond);
	    printf("%lu number of elements is assembled\n",hseList.size());
	}
}

//NOTE: Must be called after assembleFromInterface()
//      and before calling the spring solver.
void CollisionSolver3d::recordOriginalPosition()
{
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            STATE* sl = (STATE*)left_state(pt); 
            
            sl->has_collsn = false;
            sl->has_strainlim_prox = false;
            sl->collsn_dt = -1.0;

            if (isMovableRigidBody(pt))
                continue;

            for (int j = 0; j < 3; ++j)
            {
                sl->x_old[j] = Coords(pt)[j];
            
                if (std::isnan(sl->x_old[j]))
                {
                    std::cout << "nan_x_old" << std::endl;
                    clean_up(ERROR);
                }
            }
	    }
	}
}

void CollisionSolver3d::setDomainBoundary(double* L, double* U)
{
	for (int i = 0; i < m_dim; ++i)
    {
	    Boundary[i][0] = L[i];
	    Boundary[i][1] = U[i];
	}
}

//TODO: investigate
void CollisionSolver3d::detectDomainBoundaryCollision() {
	double dt = getTimeStepSize();
	double mu = getFabricFrictionConstant();
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
                it < hseList.end(); ++it) {
	    for (int i = 0; i < (*it)->num_pts(); ++i) {
		POINT* pt = (*it)->Point_of_hse(i);
		//if (isMovableRigidBody(pt)) continue;
                STATE* sl = (STATE*)left_state(pt);

        double cand_coords[3]; //candidate position
		//try to modify the average velocity 
		//according to the new candidate position
		double dv = 0;
		for (int j = 0; j < m_dim; ++j) {
		    cand_coords[j] = sl->x_old[j] + dt*sl->avgVel[j];
		    double L, U;
		    L = getDomainBoundary(j,0);
		    U = getDomainBoundary(j,1);
		    if (cand_coords[j] <= L) 
		    {
			sl->has_collsn = true;
			cand_coords[j] = L + s_thickness;
			dv = fabs(sl->avgVel[j]);
                        // ytb
                        for (int k = 0; k < 3; k++)
                             sl->avgVel[k] = 0.0;
		    	//sl->avgVel[j] = 0.0;
		        Coords(pt)[j] = cand_coords[j];
		    }
		    else if (cand_coords[j] >= U)
		    {
			sl->has_collsn = true;
			cand_coords[j] = U - s_thickness;
			dv = fabs(sl->avgVel[j]);
		    	sl->avgVel[j] = 0.0;
		        Coords(pt)[j] = cand_coords[j];
		    }
		}
		//reduce tangential velocity with friction
		double preVt = Mag3d(sl->avgVel);
		if (preVt > MACH_EPS)
		for (int j = 0; j < m_dim; ++j) 
		    sl->avgVel[j] *= std::max(1.0-mu*dv/preVt,0.0); 
	    }
	}
}

void CollisionSolver3d::computeAverageVelocity()
{
    double max_speed = 0.0;
    double* max_vel = nullptr;
    POINT* max_pt = nullptr;

    double dt = getTimeStepSize();

    std::vector<CD_HSE*>::iterator it;
    for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            STATE* sl = (STATE*)left_state(pt); 

            for (int j = 0; j < 3; ++j)
    	    {
                if (dt > ROUND_EPS)
                {
                    sl->avgVel[j] = (Coords(pt)[j] - sl->x_old[j])/dt;
                    sl->avgVel_old[j] = sl->avgVel[j];
                }
                else
                {
                    sl->avgVel[j] = 0.0;
                    sl->avgVel_old[j] = 0.0;
                }
                
                if (std::isnan(sl->avgVel[j]) || std::isinf(sl->avgVel[j]))
                {
                    std::cout << "nan avgVel" << std::endl;
                    printf("dt = %e, x_old = %f, x_new = %f\n",
                    dt,sl->x_old[j],Coords(pt)[j]);
                    clean_up(ERROR);
                }
		
            }

            if (debugging("average_velocity"))
            {
                if (Mag3d(sl->avgVel) >= max_speed)
                {
                    max_speed = Mag3d(sl->avgVel);
                    max_vel = sl->avgVel;
                    max_pt = pt;
                }
            }
        
        }
        
    }

    if (debugging("average_velocity"))
    {
        if (max_vel)
        {
            std::cout << "Maximum average velocity is "
                << max_vel[0] << " "
                << max_vel[1] << " "
                << max_vel[2] << std::endl; 
        }
        
        if (max_pt)
        {
            STATE* sl = (STATE*)left_state(max_pt);
            printf("x_old = [%f %f %f]\n",
                    sl->x_old[0],sl->x_old[1],sl->x_old[2]);
	        printf("x_new = [%f %f %f]\n",
                    Coords(max_pt)[0],Coords(max_pt)[1],Coords(max_pt)[2]);
	        printf("dt = %f\n",dt);
            printf("Gindex(max_pt) = %d\n",Gindex(max_pt));
        }
    }

    resetPositionCoordinates();
}

void CollisionSolver3d::resetPositionCoordinates()
{
    std::vector<CD_HSE*>::iterator it;
    for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            STATE* sl = (STATE*)left_state(pt);

            for (int j = 0; j < m_dim; ++j)
                Coords(pt)[j] =  sl->x_old[j];
        }
    }
}

void CollisionSolver3d::turnOffImpZone(){s_detImpZone = false;}
void CollisionSolver3d::turnOnImpZone(){s_detImpZone = true;}
bool CollisionSolver3d::getImpZoneStatus(){return s_detImpZone;}

void CollisionSolver3d::computeImpactZone()
{
    std::cout<<"Starting compute Impact Zone: "<<std::endl;

	int niter = 0;
	bool is_collision = true;

	turnOnImpZone();
    
    while(is_collision)
    {
        is_collision = false;

        start_clock("dynamic_AABB_collision");
        aabbCollision();
        abt_collision->query();
        stop_clock("dynamic_AABB_collision");

        is_collision = abt_collision->getCollsnState();

        if (debugging("collision"))
        {
            infoImpactZones();

            std::cout << "    #"<<niter++ << ": "
                      << abt_collision->getCount() 
                      << " collision pairs";
            
            if (is_collision)
            {
                std::cout << ",  avg_collsn_dt = "
                    << getAverageCollisionTime();
            }
            std::cout << std::endl;

            std::cout << "     " << numImpactZones
                      << " impact zones" << std::endl;
            std::cout << "     " << numImpactZonePoints
                      << " total impact zone points" << std::endl;
        }
    }
	
    turnOffImpZone();
}

void CollisionSolver3d::infoImpactZones()
{
	numImpactZones = 0;
	numImpactZonePoints = 0;

	unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
		    POINT* pt = (*it)->Point_of_hse(i);
		    POINT* head = findSet(pt);
		
            //skip traversed or isolated pts
            if (sorted(pt) || weight(head) == 1)
                continue;
            else
            {
                markImpactZonePoints(head);
                numImpactZonePoints += weight(head);
                numImpactZones++;
            }
	    }
	}
}

void CollisionSolver3d::markImpactZonePoints(POINT* head)
{
    POINT* p = head;
    while (p)
    {
        sorted(p) = YES;
        p = next_pt(p);
    }
}

void CollisionSolver3d::updateImpactZoneVelocityForRG()
{
	unsortHseList(hseList);

    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            
            //skip traversed or isolated pts
            if (sorted(pt) || weight(findSet(pt)) == 1) continue;
            else if (!isMovableRigidBody(pt))
            {
                sorted(pt) = YES;
                continue;
            }
            else
                updateImpactListVelocity(findSet(pt));
	    }
	}
}

//For jacobi style update of impact zones
void CollisionSolver3d::updateImpactZoneVelocity()
{
	numImpactZones = 0;
	numImpactZonePoints = 0;

	unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
		    POINT* pt = (*it)->Point_of_hse(i);
		    POINT* head = findSet(pt);
		
            //skip traversed or isolated pts
            if (sorted(pt) || weight(head) == 1)
                continue;
            else
            {
                updateImpactListVelocity(head);
                numImpactZonePoints += weight(head);
                numImpactZones++;
            }
	    }
	}
}

void CollisionSolver3d::resolveCollision()
{
	//catch floating point exception: nan/inf
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	computeAverageVelocity();

    //TODO: Fix limitStrain() and limitStrainRate() edge traversal,
    //      and implement the zero relative velocity in direction of edge
    //      constraint (currently only enforces the position constraint).
    /*
    if (!debugging("strainlim_off"))
    {
        limitStrain();
        //limitStrainRate();
    }
    */

    //static proximity handling
	detectProximity();

    /*
    if (!debugging("strainlim_off"))
    {
        limitStrainRate();
    }
    */

	//check linear trajectories for collisions
	detectCollision();

    //TODO: fix this function
	//detectDomainBoundaryCollision();

	//update position using final midstep velocity
	updateFinalPosition();

    //TODO: Verify detectProximity() should not be called here again.
    //detectProximity();
	
	updateFinalVelocity();
}

// function to perform AABB tree building, updating structure
// and query for proximity detection process
void CollisionSolver3d::aabbProximity()
{
    if (!abt_proximity)
    {
        abt_proximity =
            std::unique_ptr<AABBTree>(new AABBTree(MotionState::STATIC));

        for (auto it = hseList.begin(); it != hseList.end(); it++)
        {
            double tol = CollisionSolver3d::getFabricThickness();
            if ((*it)->type == CD_HSE_TYPE::STRING_BOND)
                tol = CollisionSolver3d::getStringThickness();

            AABB* ab = new AABB(tol,*it);
            abt_proximity->addAABB(ab);
        }
        abt_proximity->updatePointMap(hseList);
        volume = abt_proximity->getVolume();
    }
    else
    {
        abt_proximity->isProximity = false;
        abt_proximity->updateAABBTree(hseList);
        if (fabs(abt_proximity->getVolume()-volume) > vol_diff*volume)
        {
            abt_proximity->updateTreeStructure();
            volume = abt_proximity->getVolume();
            build_count_pre++;
        }
    }
}

void CollisionSolver3d::detectProximity()
{
    start_clock("dynamic_AABB_proximity");
    aabbProximity();
    abt_proximity->query();
    stop_clock("dynamic_AABB_proximity");

	if (debugging("proximity"))
    {
        std::cout << abt_proximity->getCount()
            << " proximity pairs" << std::endl;
    }

    if (abt_proximity->isProximity)
        updateAverageVelocity();

    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects
}

// AABB tree for collision detection process
void CollisionSolver3d::aabbCollision()
{
    clearCollisionTimes();

    if (!abt_collision)
    {
        abt_collision =
            std::unique_ptr<AABBTree>(new AABBTree(MotionState::MOVING));

        for (auto it = hseList.begin(); it != hseList.end(); it++)
        {
            double tol = CollisionSolver3d::getFabricRoundingTolerance();
            if ((*it)->type == CD_HSE_TYPE::STRING_BOND)
                tol = CollisionSolver3d::getStringRoundingTolerance();

            AABB* ab = new AABB(tol,*it,s_dt);
            abt_collision->addAABB(ab);
        }
        abt_collision->updatePointMap(hseList);
        volume = abt_collision->getVolume();
    }
    else
    {
        abt_collision->isCollsn = false;
        abt_collision->setTimeStep(s_dt);
        abt_collision->updateAABBTree(hseList);

        if (fabs(abt_collision->getVolume() - volume) > vol_diff * volume)
        {
            build_count_col++;
            abt_collision->updateTreeStructure();
            volume = abt_collision->getVolume();
        }
    }
}

void CollisionSolver3d::detectCollision()
{
	std::cout << "Starting collision handling: " << std::endl;
	
	const int MAX_ITER = 12;
	
    bool is_collision = true; 
	setHasCollision(false);
	
    int niter = 1;
	int cd_count = 0;
   
    while(is_collision)
    {
	    is_collision = false;
	    
        start_clock("dynamic_AABB_collision");
        aabbCollision();
        abt_collision->query();
        stop_clock("dynamic_AABB_collision");

        is_collision = abt_collision->getCollsnState();

	    if (debugging("collision"))
        {
            std::cout << "    #" << niter << ": "
                << abt_collision->getCount() 
                << " collision pairs";
            
            if (is_collision)
            {
                std::cout << ",  avg_collsn_dt = "
                    << getAverageCollisionTime();
            }
            std::cout << std::endl;
        }

        if (++niter > MAX_ITER)
            break;
	}

    start_clock("computeImpactZone");
	if (is_collision) 
	    computeImpactZone();
    stop_clock("computeImpactZone");
}

//Note: num has default value of 4,
//and first has default value of false
extern void createImpZone(POINT* pts[], int num, bool first)
{
	for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
	        if ((!first) && (isMovableRigidBody(pts[i]) || 
				 isMovableRigidBody(pts[j])))
            {
                continue;
            }
            mergePoint(pts[i],pts[j]); 
	    }
	}
}

void CollisionSolver3d::limitStrainRate()
{
	const int MAX_ITER = 5;
    for (int iter = 0; iter <= MAX_ITER; ++iter)
    {
        modifyStrainRate();
        
        if (debugging("strain_limiting"))
        {
            printf("   %d Strain Rate Edges\n",numStrainRateEdges);
        }

	    bool excess_strain = false;
        if (numStrainRateEdges > 0)
            excess_strain = true;

        if (!excess_strain)
            break;
	}
}

//TODO: FIX TRAVERSAL, and enforce relative velocity constraint
//gauss-seidel iteration
void CollisionSolver3d::modifyStrainRate()
{
    numStrainRateEdges = 0;
	double dt = getTimeStepSize();
    double TOL = strainrate_limit;
    
    unsortHseList(hseList);

    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        if (isRigidBody(*it)) continue;
	    
        POINT* p[2];
        STATE* sl[2];

        int np = (*it)->num_pts();
        for (int i = 0; i < ((np == 2) ? 1 : np); ++i)
        {
            p[0] = (*it)->Point_of_hse(i%np);	
            p[1] = (*it)->Point_of_hse((i+1)%np);

            //TODO: THIS DOES NOT WORK! MISSES EDGES.
            if (sorted(p[0]) && sorted(p[1])) continue;

            sl[0] = (STATE*)left_state(p[0]);
            sl[1] = (STATE*)left_state(p[1]);

            double x_cand0[3], x_cand1[3];
            for (int j = 0; j < 3; ++j)
            {
                x_cand0[j] = sl[0]->x_old[j] + sl[0]->avgVel[j]*dt;
                x_cand1[j] = sl[1]->x_old[j] + sl[1]->avgVel[j]*dt;
            }

		    double lnew = distBetweenCoords(x_cand0,x_cand1);
		    double lold = distBetweenCoords(sl[0]->x_old,sl[1]->x_old);

            if (fabs(lnew - lold) > TOL*lold)
            {
                double I = (fabs(lnew - lold) - TOL*lold)/dt;
                double I0, I1;

                if (lnew > lold) //Tension
                {
                    I0 = 0.5*I;
                    I1 = -0.5*I;
                }
                else             //Compression
                {
                    I0 = -0.5*I;
                    I1 = 0.5*I;
                }

                double e01[3];
                Pts2Vec(p[0],p[1],e01);
                scalarMult(1.0/lold,e01,e01);

                for (int j = 0; j < 3; ++j)
                {
                    sl[0]->avgVel[j] += I0*e01[j];
                    sl[1]->avgVel[j] += I1*e01[j];
                }
                sl[0]->has_strainlim_collsn = true;
                sl[1]->has_strainlim_collsn = true;
                numStrainRateEdges++;
            }

            sorted(p[0]) = YES;
            sorted(p[1]) = YES;
        }
    }
}

void CollisionSolver3d::limitStrain()
{
	const int MAX_ITER = 2;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        numStrainEdges = 0;

        modifyStrain();
        
        if (debugging("strain_limiting"))
        {
            printf("   %d Strain Edges\n",numStrainEdges);
        }

	    bool excess_strain = false;
        if (numStrainEdges > 0)
            excess_strain = true;

        if (!excess_strain)
            break;
	}
}

//TODO: FIX TRAVERSAL, and enforce relative velocity constraint
//jacobi iteration
void CollisionSolver3d::modifyStrain()
{
    double dt = getTimeStepSize();
    double TOL = strain_limit;
    double CTOL = 0.01;//TODO: add input parameter for this

	unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        if (isRigidBody(*it)) continue;
	    
        POINT* p[2];
        STATE* sl[2];

        int np = (*it)->num_pts();
        for (int i = 0; i < ((np == 2) ? 1 : np); ++i)
        {
            p[0] = (*it)->Point_of_hse(i%np);	
            p[1] = (*it)->Point_of_hse((i+1)%np);

            //TODO: THIS DOES NOT WORK! MISSES EDGES.
                //if (sorted(p[0]) && sorted(p[1])) continue;

            sl[0] = (STATE*)left_state(p[0]);
            sl[1] = (STATE*)left_state(p[1]);

            double x_cand0[3], x_cand1[3];
            for (int j = 0; j < 3; ++j)
            {
                x_cand0[j] = sl[0]->x_old[j] + sl[0]->avgVel[j]*dt;
                x_cand1[j] = sl[1]->x_old[j] + sl[1]->avgVel[j]*dt;
            }

		    double lnew = distBetweenCoords(x_cand0,x_cand1);
            double len0;

            if ((*it)->type == CD_HSE_TYPE::STRING_BOND)
            {
                CD_BOND* cd_bond = dynamic_cast<CD_BOND*>(*it);
                len0 = cd_bond->m_bond->length0;
            }
            else
            {
                CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(*it);
                len0 = cd_tri->m_tri->side_length0[i];
            }
            
            double strain = lnew - len0;
            if (strain > TOL*len0 || strain < -CTOL*len0)
            {
                double I0, I1;
                if (strain > TOL*len0) //Tension
                {
                    double I = (strain - TOL*len0)/dt;
                    I0 = 0.5*I;
                    I1 = -0.5*I;
                }
                else                   //Compression
                {
                    double I = -(strain + CTOL*len0)/dt;
                    I0 = -0.5*I;
                    I1 = 0.5*I;
                } 
                
                double e01[3];
                Pts2Vec(p[0],p[1],e01);
                scalarMult(1.0/lnew,e01,e01);
		    
                for (int j = 0; j < 3; ++j)
                {
                    sl[0]->strainImpulse[j] += I0*e01[j];
                    sl[1]->strainImpulse[j] += I1*e01[j];
                }

                sl[0]->strain_num++;
                sl[1]->strain_num++;
                numStrainEdges++;
            }
        
            //sorted(p[0]) = YES;
            //sorted(p[1]) = YES;
        }
    }
    
	unsortHseList(hseList);
    
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        if (isRigidBody(*it)) continue;
	    
        int np = (*it)->num_pts();
        for (int i = 0; i < np; ++i)
	    {
            POINT* p = (*it)->Point_of_hse(i);
            if (sorted(p)) continue;

		    STATE* sl = (STATE*)left_state(p);
            if (sl->strain_num > 0)
            {
                sl->has_strainlim_prox = true;
                for (int j = 0; j > 3; ++j)
                {
                    sl->avgVel[j] += sl->strainImpulse[j];
                    sl->strainImpulse[j] = 0.0;
                }
                sl->strain_num = 0;
            }

            sorted(p) = YES;
        }
    }
}

void CollisionSolver3d::updateFinalPosition()
{
    unsortHseList(hseList);
	double dt = getTimeStepSize();

    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            if (sorted(pt) || isStaticRigidBody(pt))
                continue;

            STATE* sl = (STATE*)left_state(pt);
            for (int j = 0; j < 3; ++j)
            {
                double ncoord = sl->x_old[j] + sl->avgVel[j]*dt;
                if (std::isnan(ncoord))
                {
                    printf("nan coords, x_old = %f, avgVel = %f\n",
                            sl->x_old[j],sl->avgVel[j]);
                    clean_up(ERROR);
                }
                
                Coords(pt)[j] = ncoord;
            }

            sorted(pt) = YES;
        }
	}
}

void CollisionSolver3d::updateFinalVelocity()
{
    unsortHseList(hseList);
	
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            if (sorted(pt) || isStaticRigidBody(pt))
                continue;

            STATE* sl = (STATE*)left_state(pt);
            if (!sl->has_collsn && 
                !sl->has_strainlim_prox && !sl->has_strainlim_collsn) continue;

            for (int j = 0; j < 3; ++j)
            {
                pt->vel[j] = sl->avgVel[j];
                sl->vel[j] = sl->avgVel[j];
                
                if (std::isnan(pt->vel[j]))
                {
                    printf("nan vel and avgVel\n");
                    clean_up(ERROR);
                }
            }

            sorted(pt) = YES;
        }
        
    }
    
    updateFinalForRG();
}

void CollisionSolver3d::updateFinalForRG()
{
	POINT* pt;
        STATE* sl;
        double dt = getTimeStepSize();
	std::vector<int> mrg;
	std::map<int, bool> visited;

        for (auto it = hseList.begin(); it < hseList.end(); it++) {
             for (int i = 0; i < (*it)->num_pts(); i++) {
                  pt = (*it)->Point_of_hse(i);
                  sl = (STATE*)left_state(pt);
                  if (!isMovableRigidBody(pt)) continue;

                  int rg_index = body_index(pt->hs);

                  if ((visited.count(rg_index) == 0) || (!visited[rg_index]))
                  {
                       double* com = center_of_mass(pt->hs);

                       mrg_com[rg_index] = std::vector<double>(com, com+3);
                       visited[rg_index] = true;
                  }
             }
        }

	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
             it < hseList.end(); ++it)
        {
            for (int i = 0; i < (*it)->num_pts(); ++i){
                pt = (*it)->Point_of_hse(i);
                sl = (STATE*)left_state(pt);
                if (!isMovableRigidBody(pt)) continue;
		int rg_index = body_index(pt->hs);
                if (sl->has_collsn && 
		    std::find(mrg.begin(), mrg.end(), rg_index) == mrg.end())
                {
                    mrg.push_back(rg_index); 
                    for (int j = 0; j < 3; ++j)
                    {
                        center_of_mass_velo(pt->hs)[j] = sl->avgVel[j];
                        center_of_mass(pt->hs)[j] = sl->avgVel[j] * dt + 
                                                (mrg_com[rg_index])[j];
                    }
		    visited[rg_index] = false;
		    if (debugging("rigid_body"))
		    {
			printf("After collision handling: \n");
			printf("Body Index: %d\n", rg_index);
			printf("center_of_mass = %f %f %f\n", 
				center_of_mass(pt->hs)[0], 
				center_of_mass(pt->hs)[1], 
				center_of_mass(pt->hs)[2]);
			printf("center_of_mass_velo = %f %f %f\n", 
				center_of_mass_velo(pt->hs)[0], 
				center_of_mass_velo(pt->hs)[1], 
				center_of_mass_velo(pt->hs)[2]);
		    }
                }
		if ((visited.count(rg_index) == 0) || (!visited[rg_index]))
		{
		    double* com = center_of_mass(pt->hs);
		    mrg_com[rg_index] = std::vector<double>(com, com+3);
		    visited[rg_index] = true;
		}
            }
        }
}

//For Jacobi velocity update
void CollisionSolver3d::updateAverageVelocity()
{
	POINT *p;
	STATE *sl;
	double maxSpeed = 0;
	double* maxVel = nullptr;

    unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        if (isRigidBody(*it)) continue;
	    
        int np = (*it)->num_pts(); 
	    for (int j = 0; j < np; ++j)
	    {
            p = (*it)->Point_of_hse(j);
            
            if (sorted(p) || isStaticRigidBody(p))
                continue;

            sl = (STATE*)left_state(p);
            if (sl->collsn_num > 0)
            {
                sl->has_collsn = true;

                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] +=
                        (sl->collsnImpulse[k] + sl->friction[k])/((double)sl->collsn_num);
                
                    if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
                    {
                        printf("inf/nan vel[%d]: impulse = %f, friction = %f, collsn_num = %d\n",
                        k,sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
                        clean_up(ERROR);
                    }
                
                    sl->collsnImpulse[k] = 0.0;
                    sl->friction[k] = 0.0;
                }

                sl->collsn_num = 0;
            }

            if (sl->collsn_num_RG > 0)
            {
                sl->has_collsn = true;
                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] += sl->collsnImpulse_RG[k]/((double)sl->collsn_num_RG);
                }
                sl->collsn_num_RG = 0;
            }

            if (debugging("average_velocity"))
            {
                double speed = Mag3d(sl->avgVel);
                if (speed > maxSpeed)
                {
                    maxVel = sl->avgVel;
                    maxSpeed = speed;
                }
            }
		
            sorted(p) = YES;
	    }
	}
	
    /*
    //TODO: Moved into detect_proximity(), if successful can remove.
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects
    */

	if (debugging("average_velocity"))
    {
	    if (maxVel != nullptr)
        {
	        printf("    max velocity = [%f %f %f]\n",
                    maxVel[0],maxVel[1],maxVel[2]);
        }
    }
}

bool getCollision(const CD_HSE* a, const CD_HSE* b)
{
	const CD_BOND *cd_b1, *cd_b2;
	const CD_TRI  *cd_t1, *cd_t2;

	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
            return false;
	    return MovingTriToTri(t1,t2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return MovingBondToBond(b1,b2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
	return false;
}

//This is checking the geometric primitive for intersection
bool getProximity(const CD_HSE* a, const CD_HSE* b)
{
	const CD_BOND *cd_b1, *cd_b2;
	const CD_TRI  *cd_t1, *cd_t2;

	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
            return false;
	    return TriToTri(t1,t2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return BondToBond(b1,b2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
	return false;
}

void CollisionSolver3d::printDebugVariable()
{
	std::cout << "Enter EdgeToEdge " << edg_to_edg 
		  << " times"<< std::endl;
	std::cout << "Enter PointToTri " << pt_to_tri 
		  << " times"<< std::endl;
	std::cout << "Enter isCoplanar " << is_coplanar
		  << " times"<< std::endl;
	
    moving_edg_to_edg = moving_pt_to_tri = is_coplanar = 0;
	edg_to_edg = pt_to_tri = 0;
}


/*******************************
* utility functions start here *
*******************************/

/* The followings are helper functions for vector operations. */
void Pts2Vec(const POINT* p1, const POINT* p2, double* v){
	for (int i = 0; i < 3; ++i)	
	    v[i] = Coords(p2)[i] - Coords(p1)[i];
}

double distBetweenCoords(double* x1, double* x2)
{
	double dist = 0.0;
	for (int i = 0; i < 3; ++i){
		dist += sqr(x1[i]-x2[i]);
	}
	return std::sqrt(dist);
}

void addVec(double* v1, double* v2, double* ans)
{
	for (int i = 0; i < 3; ++i)
	    ans[i] = v1[i]+v2[i];
}

void minusVec(double* v1, double* v2, double* ans)
{
	for (int i = 0; i < 3; ++i)
	    ans[i] = v1[i]-v2[i];
}

void scalarMult(double a, double* v, double* ans)
{
	for (int i = 0; i < 3; ++i)
            ans[i] = a*v[i];	
}

extern double myDet3d(double a[3][3]){
    return  a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) 
	  - a[0][1]*(a[1][0]*a[2][2] - a[2][0]*a[1][2]) 
	  + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]);
}

void unsortHseList(std::vector<CD_HSE*>& hseList)
{
    std::vector<CD_HSE*>::iterator it; 
    for (it = hseList.begin(); it < hseList.end(); ++it)
	{
	    int np = (*it)->num_pts();
	    for (int i = 0; i < np; ++i)
            sorted((*it)->Point_of_hse(i)) = NO;
            
        if ((*it)->type == CD_HSE_TYPE::FABRIC_TRI)
        {
            CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(*it);
            sorted(cd_tri->m_tri) = NO;
        }
	}
}

void unsort_surface_point(SURFACE *surf)
{
    TRI* tri = first_tri(surf);
    for (tri; !at_end_of_tri_list(tri,surf); tri = tri->next)
    {
        for (int i = 0; i < 3; ++i)
        {
            POINT* p = Point_of_tri(tri)[i];
            sorted(p) = NO;
        }
    }
}       /* end unsort_surface_point */



//functions for UF alogrithm
int& weight(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.num_pts;
}

inline POINT*& root(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.root;
}

POINT*& next_pt(POINT* p){
	STATE* sl = (STATE*)left_state(p);
        return sl->impZone.next_pt;
}

inline POINT*& tail(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.tail;
}

extern void makeSet(std::vector<CD_HSE*>& hseList)
{
	STATE* sl;
	POINT* pt;
    for (std::vector<CD_HSE*>::iterator it = hseList.begin();
            it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sorted(pt) = NO;
            sl = (STATE*)left_state(pt);
            sl->impZone.next_pt = NULL;
            sl->impZone.tail = pt;
            sl->impZone.root = pt;
            sl->impZone.num_pts = 1;
        }
    }
}

POINT* findSet(POINT* p){
	if (root(p) != p)
		root(p) = findSet(root(p));
	return root(p);
}

void mergePoint(POINT* X, POINT* Y){
	POINT* PX = findSet(X);
	POINT* PY = findSet(Y);
	if (PX == PY) return;
	if (weight(PX) > weight(PY)){
	    //update root after merge
	    weight(PX) += weight(PY);
	    root(PY) = PX;
	    //link two list, update tail
	    next_pt(tail(PX)) = PY;
	    tail(PX) = tail(PY); 
	}
	else{
	    //update root after merge
	    weight(PY) += weight(PX);
	    root(PX) = PY;
	    //link two list, update tail
	    next_pt(tail(PY)) = PX;
	    tail(PY) = tail(PX); 
	}
}
//end of UF functions

void printPointList(POINT** plist,const int n){
	for (int i = 0; i < n; ++i){
	    printf("pt[%d] = [%f %f %f]\n",i,Coords(plist[i])[0],
		Coords(plist[i])[1],Coords(plist[i])[2]);
	}
}

bool isStaticRigidBody(const POINT* p){
    STATE* sl = (STATE*)left_state(p);
    return sl->is_fixed;
}

bool isStaticRigidBody(const CD_HSE* hse){
    for (int i = 0; i < hse->num_pts(); ++i)
   	if (isStaticRigidBody(hse->Point_of_hse(i)))
	    return true;
    return false;
}

bool isMovableRigidBody(const POINT* p){
    STATE* sl = (STATE*)left_state(p);
    return sl->is_movableRG;
}

bool isMovableRigidBody(const CD_HSE* hse){
    for (int i = 0; i < hse->num_pts(); ++i)
        if (isMovableRigidBody(hse->Point_of_hse(i)))
            return true;
    return false;
}

bool isRigidBody(const POINT* p){
    return isStaticRigidBody(p) || isMovableRigidBody(p);
}

bool isRigidBody(const CD_HSE* hse){
    return isStaticRigidBody(hse) || isMovableRigidBody(hse);
}

extern void SpreadImpactZoneImpulse(
        POINT* p,
        double impulse,
        double* nor)
{
        POINT* root = findSet(p);
        while (root)
        {
            STATE *sl = (STATE*)left_state(root);
            for (int i = 0; i < 3; ++i)
                sl->collsnImpulse_RG[i] += impulse * nor[i];
            sl->collsn_num_RG += 1;
            root = next_pt(root);
        }
}
