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

int CollisionSolver3d::tstep;
std::string CollisionSolver3d::outdir;

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

void CollisionSolver3d::clearCollisionTimes()
{
    CollisionTimes.clear();
}

void CollisionSolver3d::setSizeCollisionTimes(unsigned int size)
{
    CollisionTimes.reserve(size);
}

void CollisionSolver3d::addCollisionTime(double collsn_dt)
{
    CollisionTimes.push_back(collsn_dt);
}

double CollisionSolver3d::getAverageCollisionTime()
{
    double avg_dt =
        std::accumulate(CollisionTimes.begin(),CollisionTimes.end(),0.0);
    avg_dt /= CollisionTimes.size();
    return avg_dt;
}


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

const std::vector<CD_HSE*>& CollisionSolver3d::getHseList() const
{
    return hseList;
}

void CollisionSolver3d::initializeSystem(Front* front)
{
    ft = front;
    setStep(front->step);
    setTimeStepSize(front->dt);
    setOutputDirectory(OutName(front));
    assembleFromInterface(front->interf);
    recordOriginalPosition();
}

//NOTE: Must be called before calling the spring solver
void CollisionSolver3d::assembleFromInterface(INTERFACE* intfc)
{
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
	    unsort_surf_point(*s);
	    
        surf_tri_loop(*s,tri)
	    {
            CD_HSE_TYPE tag;

            if (wave_type(*s) == ELASTIC_BOUNDARY)
            {
                tag = CD_HSE_TYPE::FABRIC_TRI;
            }
            else if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
            {
                tag = CD_HSE_TYPE::MOVABLE_RIGID_TRI;
            }
            else if (wave_type(*s) == NEUMANN_BOUNDARY)
            {
                tag = CD_HSE_TYPE::STATIC_RIGID_TRI;
            }
            else 
            {
                printf("assembleFromInterface() ERROR: "
                        "unknown surface type\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
            
            hseList.push_back(new CD_TRI(tri,tag));
		    n_tri++;
	    }
	}

	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY)
            continue; 

        unsort_curve_point(*c);

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
            sl->has_strainlim = false;
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

void CollisionSolver3d::turnOnImpZone(){s_detImpZone = true;}
void CollisionSolver3d::turnOffImpZone(){s_detImpZone = false;}
bool CollisionSolver3d::getImpZoneStatus(){return s_detImpZone;}

void CollisionSolver3d::computeImpactZone()
{
    std::cout<<"Starting compute Impact Zone: "<<std::endl;

	int niter = 0;
    const int MAXITER = 100;

	turnOnImpZone();
    
	bool is_collision = true;
    while(is_collision)
    {
        niter++;
        is_collision = false;

        start_clock("dynamic_AABB_collision");
        aabbCollision();
        abt_collision->query();
        stop_clock("dynamic_AABB_collision");

        is_collision = abt_collision->getCollsnState();

        if (debugging("collision"))
        {
            infoImpactZones();

            std::cout << "    #" << niter << ": "
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
        
        if (niter >= MAXITER)
        {
            printf("computeImpactZone(): ERROR\n\t\
                    maxiters %d without convergence!\n",
                    MAXITER);
            
            debugImpactZones();
            clean_up(EXIT_FAILURE);
        }
    }
	
    turnOffImpZone();
}

void CollisionSolver3d::debugImpactZones()
{
    std::string outdir = CollisionSolver3d::getOutputDirectory();
    
    unsortHseList(hseList);
	int numImpactZone = 0;

	for (auto it = hseList.begin(); it < hseList.end(); ++it)
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
                //markImpactZonePoints(head);
                std::vector<POINT*> impactzone_pts;
                std::string fname = outdir + "/impzone-" +
                    std::to_string(numImpactZone);
                
                printf("Impact Zone #%d -- %d points",
                        numImpactZone,weight(head));
                
                POINT* p = head;
                while (p)
                {
                    double* coords = Coords(p);
                    printf("\t\tGindex = %ld coords = %g %g %g\n",
                            Gindex(p),coords[0],coords[1],coords[2]);
                    
                    sorted(p) = YES;
                    impactzone_pts.push_back(p);
                    p = next_pt(p);
                }
                printf("\n\n");
                vtk_write_pointset(impactzone_pts,fname,numImpactZone);
                numImpactZone++;
            }
        }
    }
    
    FT_Save(ft);
    FT_Draw(ft);
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
	//unsortHseList(hseList);
    
    /*
    auto rgbList = getHseTypeList(CD_HSE_TYPE::MOVABLE_RIGID_TRI);
    unsortHseList(rgbList);
    
    //TODO: Try this. Sub in rgbList below and remove the if blocks
    //      we no longer need.
    */

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
        //feenableexcept(FE_INVALID | FE_OVERFLOW);

    start_clock("computeAverageVelocity");
	computeAverageVelocity();
    stop_clock("computeAverageVelocity");

    // Apply impulses to enforce strain and strain rate limiting
    // distance/positional constraints on adjacent mesh vertices. 
    if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
    {
        //limitStrainRatePosn();
        limitStrainPosn();
    }

    // Static proximity handling
    start_clock("detectProximity");
	detectProximity();
    stop_clock("detectProximity");

	// Check linear trajectories for collisions
    start_clock("detectCollision");
	detectCollision();
    stop_clock("detectCollision");


    /*
    //TODO: function needs fixing -- is this even worth correcting/using?
    detectDomainBoundaryCollision();
    */


	//update position using final midstep velocity
	updateFinalPosition();

    // Zero out the relative velocity between adjacent mesh vertices
    // in the direction of their connecting edge.
    if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
    {
        limitStrainVel();
    }

	updateFinalVelocity();

    /*
    //Consolidates updateFinalPosition() and updateFinalVelocity()
    //in order to avoid a second traversal of the points. However,
    //when strain limiting is eventually implemented, the separate
    //traversals may be necessary to adjust the velocities based on
    //the final positions in order to limit the strain/strain rate.
    
    start_clock("updateFinalStates");
    updateFinalStates(); 
    stop_clock("updateFinalStates");
    */

    updateFinalForRG();
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
	setHasCollision(false);//TODO: can remove?
	
    int niter = 0;
	int cd_count = 0;
   
    while(is_collision)
    {
        niter++;
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

        if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
        {
            //TODO: Write a gauss-seidel version of
            //      limitStrainRatePosn() for use here.
            limitStrainRatePosn();
        }

        if (niter >= MAX_ITER) break;
	}

    start_clock("computeImpactZone");
	if (is_collision) 
	    computeImpactZone();
    stop_clock("computeImpactZone");

	std::cout << "End collision handling. " << std::endl;
}

//TODO: The bool first and associated default value
//      is only use for creating an impact zone for rigid
//      bodies, another function should be made specifically
//      for this initialization.
//
//Note: num has default value of 4,
//and first has default value of false
extern void createImpZone(POINT* pts[], int num, bool first)
{
    //TODO: What is the minimum number of mergePoint() calls?
	for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
            //TODO: Should it check isRigidBody() instead?
            //      Otherwise, static rigid bodies can become
            //      part of the impact zone.
            //      This might be okay though ....
            if (!first)
            {
                if (isMovableRigidBody(pts[i]) ||
                    isMovableRigidBody(pts[j])) continue;
                //TODO: In this scenario a fabric element collides
                //      with a movable rigid body element, and only
                //      the fabric points are merged into an impact
                //      zone. This doesn't make sense because the
                //      fabric element can't collide with itself.
            }

            mergePoint(pts[i],pts[j]); 
	    }
	}
}

void createImpactZone(POINT* pts[], int num)
{
    //TODO: What is the minimum number of mergePoint() calls?
    //
    //      For TRI and POINT should only require 3:
    //
    //          mergePoint(pts[0],pts[1]);
    //          mergePoint(pts[0],pts[2]);
    //          mergePoint(pts[0],pts[3]);
    //
    //      For BOND and BOND should only require 3:
    //
    //          mergePoint(pts[0],pts[1]);
    //          mergePoint(pts[2],pts[3]);
    //          mergePoint(pts[0],pts[2]);
    //
    //      The current looping structure makes 6 calls to mergePoint().
    //      Each findSet() call within mergePoint() is O(n) in the worst
    //      case, i.e. when impact zones have grown large.
	
    for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
            if (isMovableRigidBody(pts[i]) ||
                isMovableRigidBody(pts[j])) continue;
            
            mergePoint(pts[i],pts[j]);
	    }
	}
}

void createImpactZoneForRigidBody(POINT* pts[], int num)
{
	for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
            if (!isMovableRigidBody(pts[i]) &&
                !isMovableRigidBody(pts[j])) continue;

            mergePoint(pts[i],pts[j]); 
	    }
	}
}

std::vector<CD_HSE*> CollisionSolver3d::getHseTypeList(CD_HSE_TYPE type)
{
    std::vector<CD_HSE*> hseTypeList(hseList.size());
    
    auto is_hsetype = [=](CD_HSE* hse)
    {
        return hse->type == type;
    };
    
    auto it_last = std::copy_if(hseList.begin(),hseList.end(),
            hseTypeList.begin(),is_hsetype);

    hseTypeList.resize(std::distance(hseTypeList.begin(),it_last));
    
    return hseTypeList;
}

//jacobi iteration
void CollisionSolver3d::limitStrainPosn()
{
    auto bondList = getHseTypeList(CD_HSE_TYPE::STRING_BOND);
    auto triList = getHseTypeList(CD_HSE_TYPE::FABRIC_TRI);

	const int MAX_ITER = 3;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        //bool bondStrain = computeStrainImpulsesPosn(bondList);
        //bool triStrain = computeStrainImpulsesPosn(triList);
        //if (!bondStrain && !triStrain) break;
        
        int numBondStrain = computeStrainImpulsesPosn(bondList);
        int numTriStrain = computeStrainImpulsesPosn(triList);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Edges\n",numBondStrain);
            printf("%d TRI Strain Edges\n",numTriStrain);
        }

        if (numBondStrain == 0 && numTriStrain == 0) break;

        applyStrainImpulses();
	}
}

//bool CollisionSolver3d::computeStrainImpulsesPosn(std::vector<CD_HSE*>& list)
int CollisionSolver3d::computeStrainImpulsesPosn(std::vector<CD_HSE*>& list)
{
    double TOL = strain_limit;
    double dt = getTimeStepSize();

    int numStrainEdges = 0;
	for (auto it = list.begin(); it < list.end(); ++it)
    {
        POINT* p[2];
        STATE* sl[2];

        int np = (*it)->num_pts();
        int ne = ((np == 2) ? 1 : np);

        for (int i = 0; i < ne; ++i)
        {
            double len0;
            
            if ((*it)->type == CD_HSE_TYPE::STRING_BOND)
            {
                CD_BOND* cd_bond = dynamic_cast<CD_BOND*>(*it);
                len0 = cd_bond->m_bond->length0;
            }
            else if ((*it)->type == CD_HSE_TYPE::FABRIC_TRI)
            {
                //      Check the three edges of each triangle.
                //      If there isn't another triangle adjacent
                //      to the edge, operate on the edge.
                //      If there is another adjacent triangle,
                //      operate on the edge only if the current
                //      triangle has a smaller pointer than its
                //      neighbor (use global_index instead).
                //      This way, each edge is considered only once.
                //
                //          From triangle.c : writeedges()

                CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(*it);
                TRI* tri = cd_tri->m_tri;
                TRI* tri_nb = Tri_on_side(tri,i);
                
                //TODO: unsure about is_side_bdry(tri,i)
                    //if (!is_side_bdry(tri,i) && Tri_on_side(tri,i) != nullptr)
                
                if (tri_nb != nullptr)
                {
                    if (Gindex(tri_nb) < Gindex(tri)) continue;
                }
                
                len0 = cd_tri->m_tri->side_length0[i];
            }
            else
            {
                printf("computeStrainImpulsesPosn() ERROR: "
                        "unknown CD_HSE_TYPE\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
            

            p[0] = (*it)->Point_of_hse(i%np);
            p[1] = (*it)->Point_of_hse((i+1)%np);

            sl[0] = (STATE*)left_state(p[0]);
            sl[1] = (STATE*)left_state(p[1]);

            double x_cand0[3], x_cand1[3];
            for (int j = 0; j < 3; ++j)
            {
                x_cand0[j] = sl[0]->x_old[j] + sl[0]->avgVel[j]*dt;
                x_cand1[j] = sl[1]->x_old[j] + sl[1]->avgVel[j]*dt;
            }

            double lnew = distBetweenCoords(x_cand0,x_cand1);

            double delta_len0 = lnew - len0;
            
            //if (delta_len0 > TOL*len0)
            if (delta_len0 > TOL*len0 || delta_len0 < -0.25*TOL*len0)
            {
                double I;
                if (delta_len0 > TOL*len0) //Tension
                { 
                    I = 0.5*(delta_len0 - TOL*len0)/dt;
                }
                else                       //Compression
                {
                    I = 0.5*(delta_len0 + 0.25*TOL*len0)/dt;
                }

                double vec01[MAXD];
                Pts2Vec(p[0],p[1],vec01);
                scalarMult(1.0/lnew,vec01,vec01);
                
                for (int j = 0; j < 3; ++j)
                {
                    sl[0]->strainImpulse[j] += I*vec01[j];
                    sl[1]->strainImpulse[j] -= I*vec01[j];
                }
                
                sl[0]->strain_num++;
                sl[1]->strain_num++;
                numStrainEdges++;
            }
        }
    }
    
    /*
    if (debugging("strain_limiting"))
    {
        printf("%d Strain Edges\n",numStrainEdges);
    }

    return numStrainEdges > 0;
    */

    return numStrainEdges;
}

//jacobi iteration
void CollisionSolver3d::limitStrainRatePosn()
{
    auto bondList = getHseTypeList(CD_HSE_TYPE::STRING_BOND);
    auto triList = getHseTypeList(CD_HSE_TYPE::FABRIC_TRI);

	const int MAX_ITER = 3;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        //bool bondStrainRate = computeStrainRateImpulsesPosn(bondList);
        //bool triStrainRate = computeStrainRateImpulsesPosn(triList);
        //if (!bondStrainRate && !triStrainRate) break;
        
        int numBondStrainRate = computeStrainRateImpulsesPosn(bondList);
        int numTriStrainRate = computeStrainRateImpulsesPosn(triList);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Rate Edges\n",numBondStrainRate);
            printf("%d TRI Strain Rate Edges\n",numTriStrainRate);
        }

        if (numBondStrainRate == 0 && numTriStrainRate == 0) break;

        applyStrainImpulses();
	}
}

//bool CollisionSolver3d::computeStrainRateImpulsesPosn(std::vector<CD_HSE*>& list)
int CollisionSolver3d::computeStrainRateImpulsesPosn(std::vector<CD_HSE*>& list)
{
    double TOL = strain_limit;
    double dt = getTimeStepSize();

    int numStrainRateEdges = 0;
	for (auto it = list.begin(); it < list.end(); ++it)
    {
        POINT* p[2];
        STATE* sl[2];

        int np = (*it)->num_pts();
        int ne = ((np == 2) ? 1 : np);

        for (int i = 0; i < ne; ++i)
        {
            if ((*it)->type == CD_HSE_TYPE::FABRIC_TRI)
            {
                CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(*it);
                TRI* tri = cd_tri->m_tri;
                TRI* tri_nb = Tri_on_side(tri,i);
                
                if (tri_nb != nullptr)
                {
                    if (Gindex(tri_nb) < Gindex(tri))
                        continue;
                }
            }

            p[0] = (*it)->Point_of_hse(i%np);
            p[1] = (*it)->Point_of_hse((i+1)%np);

            sl[0] = (STATE*)left_state(p[0]);
            sl[1] = (STATE*)left_state(p[1]);

            double x_cand0[3], x_cand1[3];
            for (int j = 0; j < 3; ++j)
            {
                x_cand0[j] = sl[0]->x_old[j] + sl[0]->avgVel[j]*dt;
                x_cand1[j] = sl[1]->x_old[j] + sl[1]->avgVel[j]*dt;
            }

            double lnew = distBetweenCoords(x_cand0,x_cand1);
            double lold = distance_between_positions(sl[0]->x_old,sl[1]->x_old,3);
            double delta_lold = lnew - lold;

            //if (delta_lold > TOL*lold)
            if (fabs(delta_lold) > TOL*lold)
            {
                double I;
                if (delta_lold > TOL*lold) //Tension
                { 
                    I = 0.5*(delta_lold - TOL*lold)/dt;
                }
                else                       //Compression
                {
                    I = 0.5*(delta_lold + TOL*lold)/dt;
                }
                
                //double I = 0.5*(delta_lold - TOL*lold)/dt;

                double vec01[MAXD];
                Pts2Vec(p[0],p[1],vec01);
                scalarMult(1.0/lnew,vec01,vec01);
                
                for (int j = 0; j < 3; ++j)
                {
                    sl[0]->strainImpulse[j] += I*vec01[j];
                    sl[1]->strainImpulse[j] -= I*vec01[j];
                }
                
                sl[0]->strain_num++;
                sl[1]->strain_num++;
                numStrainRateEdges++;
            }
        }
    }
    
    /*
    if (debugging("strain_limiting"))
    {
        printf("%d Strain Rate Edges\n",numStrainRateEdges);
    }

    return numStrainRateEdges > 0;
    */

    return numStrainRateEdges;
}

//jacobi iteration
void CollisionSolver3d::limitStrainVel()
{
    auto bondList = getHseTypeList(CD_HSE_TYPE::STRING_BOND);
    auto triList = getHseTypeList(CD_HSE_TYPE::FABRIC_TRI);
	
    const int MAX_ITER = 3;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        //bool bondStrain = computeStrainImpulsesVel(bondList);
        //bool triStrain = computeStrainImpulsesVel(triList);
        //if (!bondStrain && !triStrain) break;
        
        int numBondStrainVel = computeStrainImpulsesVel(bondList);
        int numTriStrainVel = computeStrainImpulsesVel(triList);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Velocity Edges\n",numBondStrainVel);
            printf("%d TRI Strain Velocity Edges\n",numTriStrainVel);
        }

        if (numBondStrainVel == 0 && numTriStrainVel == 0) break;

        applyStrainImpulses();
	}
}

//bool CollisionSolver3d::computeStrainImpulsesVel(std::vector<CD_HSE*>& list)
int CollisionSolver3d::computeStrainImpulsesVel(std::vector<CD_HSE*>& list)
{
    double TOL = strain_limit;
    double dt = getTimeStepSize();

    int numRelVelStrainEdges = 0;
	for (auto it = list.begin(); it < list.end(); ++it)
    {
        POINT* p[2];
        STATE* sl[2];

        int np = (*it)->num_pts();
        int ne = ((np == 2) ? 1 : np);

        for (int i = 0; i < ne; ++i)
        {
            if ((*it)->type == CD_HSE_TYPE::FABRIC_TRI)
            {
                CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(*it);
                TRI* tri = cd_tri->m_tri;
                TRI* tri_nb = Tri_on_side(tri,i);
                
                if (tri_nb != nullptr)
                {
                    if (Gindex(tri_nb) < Gindex(tri))
                        continue;
                }
            }

            p[0] = (*it)->Point_of_hse(i%np);
            p[1] = (*it)->Point_of_hse((i+1)%np);

            sl[0] = (STATE*)left_state(p[0]);
            sl[1] = (STATE*)left_state(p[1]);
            
            //skip edges that did not get strain limiting impulses
            //in limitStrainPos() or limitStrainRatePos()
            if (!(sl[0]->has_strainlim && sl[1]->has_strainlim)) continue;

            double vel_rel[MAXD];
            for (int j = 0; j < 3; ++j)
                vel_rel[j] = sl[1]->avgVel[j] - sl[0]->avgVel[j];
            
            double vec01[MAXD];
            Pts2Vec(p[0],p[1],vec01);//p1 - p0
            double len = Mag3d(vec01);
            scalarMult(1.0/len,vec01,vec01);

            //component of the relative velocity in the direction
            //of the edge joining points a and b (a-->b)
            double vcomp01 = Dot3d(vel_rel,vec01);

            //TODO: MACH_EPS may be too strict of a tolerance here
            if (fabs(vcomp01) > MACH_EPS)
            {
                double I = 0.5*vcomp01;
                
                for (int j = 0; j < 3; ++j)
                {
                    sl[0]->strainImpulse[j] += I*vec01[j];
                    sl[1]->strainImpulse[j] -= I*vec01[j];
                }
                
                sl[0]->strain_num++;
                sl[1]->strain_num++;
                numRelVelStrainEdges++;
            }
        }
    }
    
    /*
    if (debugging("strain_limiting"))
    {
        printf("%d Relative Velocity Strain Edges\n",numRelVelStrainEdges);
    }

    return numRelVelStrainEdges > 0;
    */

    return numRelVelStrainEdges;
}

//TODO: Could append bondList and triList into a
//      single list, and pass as an argument to
//      applyStrainImpulses(). Then would not need
//      to check for rigid body points.
void CollisionSolver3d::applyStrainImpulses()
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
            if (sorted(p)) continue;

            sl = (STATE*)left_state(p);
            if (sl->strain_num > 0)
            {
                sl->has_strainlim = true;

                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] += sl->strainImpulse[k]/sl->strain_num;
                
                    if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
                    {
                        printf("inf/nan vel[%d]: strain_impulse = %f, strain_num = %d\n",
                                k,sl->strainImpulse[k],sl->strain_num);
                        LOC(); clean_up(ERROR);
                    }
                
                    sl->strainImpulse[k] = 0.0;
                }

                sl->strain_num = 0;
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
	
	if (debugging("average_velocity"))
    {
	    if (maxVel != nullptr)
        {
	        printf("    max velocity = [%f %f %f]\n",
                    maxVel[0],maxVel[1],maxVel[2]);
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
            if (sorted(pt) || isStaticRigidBody(pt)) continue;

            STATE* sl = (STATE*)left_state(pt);
            for (int j = 0; j < 3; ++j)
            {
                double ncoord = sl->x_old[j] + sl->avgVel[j]*dt;
                if (std::isnan(ncoord))
                {
                    printf("nan coords, x_old = %f, avgVel = %f\n",
                            sl->x_old[j],sl->avgVel[j]);
                    LOC(); clean_up(ERROR);
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
            if (sorted(pt) || isStaticRigidBody(pt)) continue;
            sorted(pt) = YES;
            
            STATE* sl = (STATE*)left_state(pt);
            if (!sl->has_collsn && !sl->has_strainlim) continue;

            for (int j = 0; j < 3; ++j)
            {
                pt->vel[j] = sl->avgVel[j];
                sl->vel[j] = sl->avgVel[j];
                
                if (std::isnan(pt->vel[j]))
                {
                    printf("nan vel and avgVel\n");
                    LOC(); clean_up(ERROR);
                }
            }
        }
        
    }
    
}

//Consolidate updateFinalPosition() and updateFinalVelocity()
//since avgVel is not further modified by updateFinalVelocity()
//by solving the implicit equation needed for central differencing
//as in the Bridson and Fedkiw paper.
void CollisionSolver3d::updateFinalStates()
{
    unsortHseList(hseList);
	double dt = getTimeStepSize();

    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            if (sorted(pt) || isStaticRigidBody(pt)) continue;

            STATE* sl = (STATE*)left_state(pt);
            for (int j = 0; j < 3; ++j)
            {
                double ncoord = sl->x_old[j] + sl->avgVel[j]*dt;
                if (std::isnan(sl->avgVel[j]) || std::isnan(ncoord))
                {
                    printf("CollisionSolver3d::updateFinalState() ERROR:\n");
                    printf("\tx_old = %f\n \tavgVel = %f\n",
                            sl->x_old[j],sl->avgVel[j]);
                    LOC(); clean_up(ERROR);
                }
                Coords(pt)[j] = ncoord;
            }
            
            sorted(pt) = YES;
            if (!sl->has_collsn && !sl->has_strainlim) continue;
            
            for (int j = 0; j < 3; ++j)
            {
                pt->vel[j] = sl->avgVel[j];
                sl->vel[j] = sl->avgVel[j];
            }
        }
	}
    
}

//TODO: Figure out what's going on here.
void CollisionSolver3d::updateFinalForRG()
{
	POINT* pt;
    STATE* sl;
    double dt = getTimeStepSize();
	std::map<int,bool> visited;
	std::vector<int> mrg;

    for (auto it = hseList.begin(); it < hseList.end(); ++it)
    {
         for (int i = 0; i < (*it)->num_pts(); i++)
         {
              pt = (*it)->Point_of_hse(i);
              if (!isMovableRigidBody(pt)) continue;

              int rg_index = body_index(pt->hs);
              if ((visited.count(rg_index) == 0) || (!visited[rg_index]))
              {
                   double* com = center_of_mass(pt->hs);
                   mrg_com[rg_index] = std::vector<double>(com,com+3);
                   visited[rg_index] = true;
              }
         }
    }

    for (auto it = hseList.begin(); it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            if (!isMovableRigidBody(pt)) continue;

            sl = (STATE*)left_state(pt);
            int rg_index = body_index(pt->hs);
    
            if (sl->has_collsn &&
                std::find(mrg.begin(),mrg.end(),rg_index) == mrg.end())
            {
                mrg.push_back(rg_index); 
                for (int j = 0; j < 3; ++j)
                {
                    center_of_mass_velo(pt->hs)[j] = sl->avgVel[j];
                    center_of_mass(pt->hs)[j] = 
                        sl->avgVel[j]*dt + (mrg_com[rg_index])[j];
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
                mrg_com[rg_index] = std::vector<double>(com,com+3);
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
                        (sl->collsnImpulse[k] + sl->friction[k])/sl->collsn_num;
                
                    if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
                    {
                        printf("inf/nan avgVel[%d]: impulse = %f,"
                               " friction = %f, collsn_num = %d\n",k,
                               sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
                        LOC(); clean_up(ERROR);
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
                    sl->avgVel[k] += sl->collsnImpulse_RG[k]/sl->collsn_num_RG;
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

//TODO: If impact zone handling is enabled, should all
//      points of the hypersurface elements be added to
//      an impact zone all at once, rather than only adding
//      4 interfering points at a time in the point to tri
//      and edge to edge tests?
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
	}
}

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

//TODO: verify anything related to this and impact zones
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

POINT* findSet(POINT* p)
{
	if (root(p) != p)
        return root(p) = findSet(root(p));
    return p;
}

//TODO: Don't think next_pt pointers being properly updated.
void mergePoint(POINT* X, POINT* Y)
{
    POINT* PX = findSet(X);
	POINT* PY = findSet(Y);
	if (PX == PY) return;
	
    if (weight(PX) > weight(PY))
    {
	    //update root after merge
	    weight(PX) += weight(PY);
	    root(PY) = PX;
	    //link two list, update tail
	    next_pt(tail(PX)) = PY;
	    tail(PX) = tail(PY); 
	}
	else
    {
	    //update root after merge
	    weight(PY) += weight(PX);
	    root(PX) = PY;
	    //link two list, update tail
	    next_pt(tail(PY)) = PX;
	    tail(PY) = tail(PX); 
	}
}
//end of UF functions

void printPointList(POINT** plist,const int n)
{
	for (int i = 0; i < n; ++i)
    {
	    printf("pt[%d] = [%f %f %f]\n",i,Coords(plist[i])[0],
		Coords(plist[i])[1],Coords(plist[i])[2]);
	}
}

bool isStaticRigidBody(const POINT* p)
{
    STATE* sl = (STATE*)left_state(p);
    return sl->is_fixed;
}

bool isStaticRigidBody(const CD_HSE* hse)
{
    for (int i = 0; i < hse->num_pts(); ++i)
    {
        if (isStaticRigidBody(hse->Point_of_hse(i)))
            return true;
    }
    return false;
}

bool isMovableRigidBody(const POINT* p)
{
    STATE* sl = (STATE*)left_state(p);
    return sl->is_movableRG;
}

bool isMovableRigidBody(const CD_HSE* hse)
{
    for (int i = 0; i < hse->num_pts(); ++i)
    {
        if (isMovableRigidBody(hse->Point_of_hse(i)))
            return true;
    }
    return false;
}

bool isRigidBody(const POINT* p)
{
    return isStaticRigidBody(p) || isMovableRigidBody(p);
}

bool isRigidBody(const CD_HSE* hse)
{
    return isStaticRigidBody(hse) || isMovableRigidBody(hse);
}

//HACK, should never be used
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
