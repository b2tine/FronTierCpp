#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <FronTier.h>
#include "AABB.h"
#include "collid.h"

#include <omp.h>

/*****declaration of static functions starts here********/
//static void makeSet(std::vector<CD_HSE*>&);
static POINT* findSet(POINT*);
static void mergePoint(POINT*,POINT*);
inline POINT*& root(POINT*);
inline POINT*& tail(POINT*);
/*******************end of declaration*******************/

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;
typedef Kernel::Triangle_3                                    Triangle_3;

//define default parameters for collision detection
bool   CollisionSolver::s_detImpZone = false;
double CollisionSolver::s_eps = EPS;
double CollisionSolver::s_thickness = 0.0001;
double CollisionSolver::s_dt = DT;
double CollisionSolver::s_k = 1000;
double CollisionSolver::s_m = 0.01;
double CollisionSolver::s_lambda = 0.02;
double CollisionSolver::s_cr = 0.0;

//debugging variables
int CollisionSolver::moving_edg_to_edg = 0;
int CollisionSolver::moving_pt_to_tri = 0;
int CollisionSolver::is_coplanar = 0;
int CollisionSolver::edg_to_edg = 0;
int CollisionSolver::pt_to_tri = 0;

// To use Pimpl idiom by unique_ptr, special member function 
// should be explicit declared in header file and defined in 
// implementation file
CollisionSolver::CollisionSolver(int dim)
    : m_dim(dim), abt_proximity(nullptr)
{}

CollisionSolver::CollisionSolver() = default;
CollisionSolver::CollisionSolver(CollisionSolver&& rhs) = default;
CollisionSolver& CollisionSolver::operator=(CollisionSolver&& rhs) = default;
CollisionSolver::~CollisionSolver() = default;

void CollisionSolver::clearHseList(){
	for (unsigned i = 0; i < hseList.size(); ++i){
		delete hseList[i];
	}
	hseList.clear();
}
//set rounding tolerance
void CollisionSolver::setRoundingTolerance(double neweps)
{
	s_eps = neweps;
}

double CollisionSolver::getRoundingTolerance(){return s_eps;}

//set fabric thickness
void CollisionSolver::setFabricThickness(double h){s_thickness = h;}
double CollisionSolver::getFabricThickness(){return s_thickness;}
double CollisionSolver::setVolumeDiff(double vd) { vol_diff = vd; }

//this function should be called at every time step
void CollisionSolver::setTimeStepSize(double new_dt)
{
    s_dt = new_dt;
}

double CollisionSolver::getTimeStepSize(){return s_dt;}

//set spring constant
void   CollisionSolver::setSpringConstant(double new_k){s_k = new_k;}
double CollisionSolver::getSpringConstant(){return s_k;}

//set spring friction 
void   CollisionSolver::setFrictionConstant(double new_la){s_lambda = new_la;}
double CollisionSolver::getFrictionConstant(){return s_lambda;}

//set mass of fabric point
void   CollisionSolver::setPointMass(double new_m){s_m = new_m;}
double CollisionSolver::getPointMass(){return s_m;}

//set restitution coefficient between rigid bodies
void   CollisionSolver::setRestitutionCoef(double new_cr){s_cr = new_cr;}
double CollisionSolver::getRestitutionCoef(){return s_cr;}

void CollisionSolver::recordOriginalPosition(){
	POINT* pt;
	STATE* sl;
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it){
	    for (int i = 0; i < (*it)->num_pts(); ++i){
		pt = (*it)->Point_of_hse(i);
		sl = (STATE*)left_state(pt); 
		sl->has_collsn = false;
		if (isMovableRigidBody(pt)) continue;
		for (int j = 0; j < 3; ++j)
		    sl->x_old[j] = Coords(pt)[j];
		if (std::isnan(sl->x_old[0])) std::cout<<"nan_x_old"<<std::endl;
	    }
	}
}

void CollisionSolver::setDomainBoundary(double* L, double* U) {
	for (int i = 0; i < m_dim; ++i) {
	    Boundary[i][0] = L[i];
	    Boundary[i][1] = U[i];
	}
}

void CollisionSolver::detectDomainBoundaryCollision() {
	double dt = getTimeStepSize();
	double mu = getFrictionConstant();
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

void CollisionSolver::computeAverageVelocity()
{
    POINT* pt;
    STATE* sl; 
    double dt = getTimeStepSize();
    double max_speed = 0.0;
    double* max_vel = nullptr;
    POINT* max_pt=nullptr;

    /*
    printf("Test dt = %f  ROUND_EPS = %f\n",dt,ROUND_EPS);
    pt = gpoints[11622];
    sl = (STATE*)left_state(pt); 
    printf("The point coords: %f %f %f\n",Coords(pt)[0],Coords(pt)[1],
                                Coords(pt)[2]);
    printf("avgVel = %f %f %f\n",sl->avgVel[0],sl->avgVel[1],sl->avgVel[2]);
    printf("sl->x_old = %f %f %f\n",sl->x_old[0],sl->x_old[1],sl->x_old[2]);
    */

    /*
    printf("Test computeAverageVelocity() hseList.size() = %d\n",
                    hseList.size());
    */

    for (std::vector<CD_HSE*>::iterator it = hseList.begin();
            it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sl = (STATE*)left_state(pt); 
            for (int j = 0; j < 3; ++j)
    	    {
                //TODO: us dt > DT = 0.001 --- ROUND_EPS = 1.0e-10
                //if (dt > ROUND_EPS)
                if (dt > DT)
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
                    std::cout<<"nan avgVel" << std::endl;
                    printf("dt = %e, x_old = %f, x_new = %f\n",
                    dt,sl->x_old[j],Coords(pt)[j]);
                    clean_up(ERROR);
                }
		
            }

            if (debugging("collision"))
            {
                if (Mag3d(sl->avgVel) >= max_speed)
                {
                    max_speed = Mag3d(sl->avgVel);
                    max_vel = sl->avgVel;
                    max_pt = pt;
                    /*
                    if (Gindex(pt) == 11622)
                    {
                        printf("\n\nmax_speed = %f  max_vel = %f %f %f\n",
                            max_speed,max_vel[0],max_vel[1],max_vel[2]);
                        printf("Gindex(p) = %d\n",Gindex(pt));
                        printf("Coords(p) = %f %f %f\n\n\n",Coords(pt)[0],
                                    Coords(pt)[1],Coords(pt)[2]);
                    }
                    */
                }
            }
        
        }
        
    }

    if (debugging("collision"))
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
            sl = (STATE*)left_state(max_pt);
            printf("x_old = [%f %f %f]\n",
                    sl->x_old[0],sl->x_old[1],sl->x_old[2]);
	        printf("x_new = [%f %f %f]\n",
                    Coords(max_pt)[0],Coords(max_pt)[1],Coords(max_pt)[2]);
	        printf("dt = %f\n",dt);
            printf("Gindex(max_pt) = %d\n",Gindex(max_pt));
        }
    }
	
    //restore coords of points to old coords !!!
    //x_old is the only valid coords for each point 
    //Coords(point) is for temporary judgement
	
    for (std::vector<CD_HSE*>::iterator it = hseList.begin();
            it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sl = (STATE*)left_state(pt);
            for (int j = 0; j < m_dim; ++j)
                Coords(pt)[j] =  sl->x_old[j];
        }
        
    }
}

void CollisionSolver::turnOffImpZone(){s_detImpZone = false;}
void CollisionSolver::turnOnImpZone(){s_detImpZone = true;}
bool CollisionSolver::getImpZoneStatus(){ return s_detImpZone;}

//this function is needed if collision still happens
//after several iterations;
void CollisionSolver::computeImpactZone()
{
	bool is_collision = true;
    int numZones = 0;
	int niter = 0;

    std::cout<<"Starting compute Impact Zone: "<<std::endl;

	turnOnImpZone();
	//makeSet(hseList); //this is done in AABBTree::updatePointMap()
    
    //int impzone_counter = 0
    //TODO: This can enter infinite loop
    while(is_collision)
    {
        //if (++impzone_counter == 20) clean_up(ERROR);
        is_collision = false;

        //start UF alogrithm
        //merge four pts if collision happens

        start_clock("dynamic_AABB_collision");
            aabbCollision();
            is_collision = abt_collision->getCollsnState();
        stop_clock("dynamic_AABB_collision");

            updateAverageVelocity();

        updateImpactZoneVelocity(numZones);
            std::cout <<"    #"<<niter++ << ": " << abt_collision->getCount() 
                      << " pair of collision tris" << std::endl;
        std::cout <<"     "<< numZones
              <<" zones of impact" << std::endl;
    }
	turnOffImpZone();
	return;
}

void CollisionSolver::updateImpactZoneVelocityForRG()
{
	POINT* pt;
	unsortHseList(hseList);

	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            
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

void CollisionSolver::updateImpactZoneVelocity(int &nZones)
{
	POINT* pt;
	int numZones = 0;

	unsortHseList(hseList);
	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it){
	    for (int i = 0; i < (*it)->num_pts(); ++i){
		pt = (*it)->Point_of_hse(i);
		//skip traversed or isolated pts
		if (sorted(pt) ||
		    weight(findSet(pt)) == 1) continue;
		else{
		    updateImpactListVelocity(findSet(pt));
		    numZones++;
		}
	    }
	}	
	nZones = numZones;
}

//resolve collision in the input tris list
void CollisionSolver::resolveCollision()
{
	//catch floating point exception: nan/inf
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	start_clock("computeAverageVelocity");
	computeAverageVelocity();
	stop_clock("computeAverageVelocity");

	start_clock("detectProximity");
	detectProximity();
	stop_clock("detectProximity");

	if (debugging("collision"))
	    printDebugVariable();

	//check linear trajectories for collisions
	start_clock("detectCollision");
	detectCollision();
	stop_clock("detectCollision");

	if (debugging("collision"))
	    printDebugVariable();
	
	start_clock("detectDomainBoundaryCollision");
	detectDomainBoundaryCollision();
	stop_clock("detectDomainBoundaryCollision");

	//update position using final midstep velocity
	start_clock("updateFinalPosition");
	updateFinalPosition();
	stop_clock("updateFinalPosition");

    //TODO: implement this function correctly
	//start_clock("reduceSuperelast");
	    //reduceSuperelast();
	//stop_clock("reduceSuperelast");
	
	start_clock("updateFinalVelocity");
    //detectProximity();
    //TODO: implement this function correctly
	updateFinalVelocity();
	stop_clock("updateFinalVelocity");
}

// function to perform AABB tree building, updating structure
// and query for proximity detection process
void CollisionSolver::aabbProximity() {
    if (!abt_proximity) {
        double pre_tol = CollisionSolver3d::getFabricThickness();
        abt_proximity = std::unique_ptr<AABBTree>(new AABBTree(STATIC));
        for (auto it = hseList.begin(); it != hseList.end(); it++) {
             AABB* ab = new AABB(pre_tol, *it, abt_proximity->getType());
             abt_proximity->addAABB(ab);
        }
        abt_proximity->updatePointMap(hseList);
        volume = abt_proximity->getVolume();
    }
    /*
    else {
        delete abt_proximity.release();
        double pre_tol = CollisionSolver3d::getFabricThickness();

        abt_proximity = std::move(std::make_unique<AABBTree>(STATIC));
        for (auto it = hseList.begin(); it != hseList.end(); it++) {
             AABB* ab = new AABB (pre_tol, *it, abt_proximity->getType());
             abt_proximity->addAABB(ab);
        }
        abt_proximity->updatePointMap(hseList);
        volume = abt_proximity->getVolume();
    }
    */
    else {
        abt_proximity->updateAABBTree(hseList);
        // if current tree structure doesn't fit for the current 
        // surface, update structure of the tree
        if (fabs(abt_proximity->getVolume() - volume) > vol_diff * volume) {
            abt_proximity->updateTreeStructure();
            volume = abt_proximity->getVolume();
            build_count_pre++;
            std::cout << "build_count_pre is " << build_count_pre << std::endl; 
        }
    }
    // query for collision detection of AABB elements
    abt_proximity->query(this);
}

void CollisionSolver::detectProximity()
{
    start_clock("dynamic_AABB_proximity");
    aabbProximity();
    stop_clock("dynamic_AABB_proximity");

	updateAverageVelocity();
	if (debugging("collision"))
        std::cout << abt_proximity->getCount()
            << " pair of proximity" << std::endl;
}

// AABB tree for collision detection process
void CollisionSolver::aabbCollision() {
    if (!abt_collision) {
        abt_collision = std::unique_ptr<AABBTree>(new AABBTree(MOVING));
        for (auto it = hseList.begin(); it != hseList.end(); it++) {
             AABB* ab = new AABB(*it, abt_collision->getType(), s_dt);
             abt_collision->addAABB(ab);
        }
        abt_collision->updatePointMap(hseList);
        volume = abt_collision->getVolume();
    }
    /*
    else {
        delete abt_collision.release();
        abt_collision = std::move(std::make_unique<AABBTree>(MOVING));
        for (auto it = hseList.begin(); it != hseList.end(); it++) {
             AABB* ab = new AABB (*it, abt_collision->getType(), s_dt);
             abt_collision->addAABB(ab);
        }
        abt_collision->updatePointMap(hseList);
        volume = abt_collision->getVolume();
    }
    */
    else {
        abt_collision->setTimeStep(s_dt);
        abt_collision->updateAABBTree(hseList);
        if (fabs(abt_collision->getVolume() - volume) > vol_diff * volume) {
            build_count_col++;
            abt_collision->updateTreeStructure();
            volume = abt_collision->getVolume();
            std::cout << "build_count_col is " << build_count_col << std::endl; 
        }
    }
    abt_collision->query(this);
}

void CollisionSolver::detectCollision()
{
	bool is_collision = true; 
	const int MAX_ITER = 8;
	int niter = 1;

	std::cout<<"Starting collision handling: "<<std::endl;
	//record if has an actual collision
	//this is useful for adpative dt
	int cd_count = 0;
	setHasCollision(false);
	//
	while(is_collision){
	    is_collision = false;
	    
        start_clock("dynamic_AABB_collision");
	    aabbCollision();
            is_collision = abt_collision->getCollsnState();
	    stop_clock("dynamic_AABB_collision");

	    if (cd_count++ == 0 && is_collision)
                setHasCollision(true);

	    updateAverageVelocity();
	    std::cout<<"    #"<<niter << ": " << abt_collision->getCount() 
		     << " pair of collision tris" << std::endl;
	    
        if (++niter > MAX_ITER)
            break;
	}

	start_clock("computeImpactZone");
	if (is_collision) 
	    computeImpactZone();
	stop_clock("computeImpactZone");
}

/*
void CollisionSolver::detectCollision()
{
	bool is_collision = true; 
	int MAX_ITER = 8;
	//const int MAX_ITER = 8;
	int niter = 1;

    int npairs = static_cast<int>(HUGE);
    int prev_npairs = npairs;

	std::cout<<"Starting collision handling: "<<std::endl;
	//record if has an actual collision
	//this is useful for adpative dt
	int cd_count = 0;
	setHasCollision(false);
	//
	while(is_collision) {
	    
        is_collision = false;
	    start_clock("dynamic_AABB_collision");
	    
        aabbCollision();
        is_collision = abt_collision->getCollsnState();

        stop_clock("dynamic_AABB_collision");

        //TODO: This doesn't do anything.
        //      boolean set in setHasCollision() is retrieved using
        //      hasCollision(), but is never called anywhere.
	    if (cd_count++ == 0 && is_collision)
            setHasCollision(true);

	    updateAverageVelocity();
        
        prev_npairs = npairs;
        npairs = abt_collision->getCount(); 
	    
        std::cout<<"    #"<<niter << ": " << npairs
		     << " pair of collision tris" << std::endl;

	    //if (++niter > MAX_ITER) break;
	    if (++niter > MAX_ITER)
        {
            if (npairs >= prev_npairs/2 + 1)
            {
                break;
            }
            MAX_ITER++;
        }
	}
    
    start_clock("computeImpactZone");
	if (is_collision) 
	    computeImpactZone();
	stop_clock("computeImpactZone");
}
*/

//helper function to detect a collision between 
//a moving point and a moving triangle
//or between two moving edges 
extern void createImpZone(POINT* pts[], int num, bool first){
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

//TODO: Implement this correctly.
//      jacobi iteration style for strain and
//      gauss-seidel iteration style for strain rate.
//      Should be called after collisions have been handled.
bool CollisionSolver::reduceSuperelastOnce(int& num_edges)
{
	double dt = getTimeStepSize();
	const double superelasTol = 0.10;
	bool has_superelas = false;
	num_edges = 0;
	
    for (unsigned i = 0; i < hseList.size(); ++i)
    {
	    CD_HSE* hse = hseList[i];
	    if (isRigidBody(hse)) continue;
	    
	    int np = hse->num_pts();
        for (int j = 0; j < ((np == 2) ? 1 : np); ++j)
        {
            POINT* p[2];
            STATE* sl[2];
            p[0] = hse->Point_of_hse(j%np);	
            p[1] = hse->Point_of_hse((j+1)%np);
            sl[0]= (STATE*)left_state(p[0]);
            sl[1]= (STATE*)left_state(p[1]);

            double x_cand[2][3];    
            for (int k = 0; k < 2; ++k)
            {
                double tmp[3];
                scalarMult(dt,sl[k]->avgVel,tmp);
                addVec(sl[k]->x_old,tmp,x_cand[k]);
            }
            
            double len_new = distance_between_positions(x_cand[0],x_cand[1],3);
		    double len_old = distance_between_positions(sl[0]->x_old,sl[1]->x_old,3);
		    double len0;
		
            if (CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(hse))
                len0 = cd_tri->m_tri->side_length0[j];
            else if (CD_BOND* cd_bond = dynamic_cast<CD_BOND*>(hse))
                len0 = cd_bond->m_bond->length0;
            else
            {
                std::cout<<"Unknown type"<<std::endl;
                clean_up(ERROR);
            }
            
            double vec[3], v_rel[3];
            minusVec(sl[0]->x_old,sl[1]->x_old,vec);
	        minusVec(sl[0]->avgVel,sl[1]->avgVel,v_rel);
            
            if (len_old > ROUND_EPS && len_new > ROUND_EPS)
            {
                scalarMult(1/len_old,vec,vec); //normalize
                double strain_rate = (len_new-len_old)/len_old;
                double strain = (len_new-len0)/len0;
                
                if (fabs(strain) > superelasTol || fabs(strain_rate) > superelasTol)
                {
                    double v_tmp[3];
                    addVec(sl[0]->avgVel,sl[1]->avgVel,v_tmp);
                    scalarMult(0.5,v_tmp,v_tmp);
                    memcpy((void*)sl[0]->avgVel,(void*)v_tmp,3*sizeof(double)); 
                    memcpy((void*)sl[1]->avgVel,(void*)v_tmp,3*sizeof(double));
                    num_edges++;
                    has_superelas = true;
                }
            }
            else
            {
                printf("Warning: len0 = %e, len_new = %e, len_old = %e\n",
                        len0,len_new,len_old);
                printf("p0 = %p, p1 = %p\n",(void*)p[0],(void*)p[1]);
                printf("x_old[0] = [%f %f %f]\n",sl[0]->x_old[0],sl[0]->x_old[1],sl[0]->x_old[2]);
                printf("avgVel[0] = [%f %f %f]\n",sl[0]->avgVel[0],sl[0]->avgVel[1],sl[0]->avgVel[2]);
                printf("x_old[1] = [%f %f %f]\n",sl[1]->x_old[0],sl[1]->x_old[1],sl[1]->x_old[2]);
                printf("avgVel[1] = [%f %f %f]\n",sl[1]->avgVel[0],sl[1]->avgVel[1],sl[1]->avgVel[2]);
                
                double v_tmp[3];
                addVec(sl[0]->avgVel,sl[1]->avgVel,v_tmp);
                scalarMult(0.5,v_tmp,v_tmp);
                memcpy((void*)sl[0]->avgVel,(void*)v_tmp,3*sizeof(double)); 
                memcpy((void*)sl[1]->avgVel,(void*)v_tmp,3*sizeof(double));
                num_edges++;
                has_superelas = true;
            }	
	    
        }
	
    }

	return has_superelas;
}

void CollisionSolver::updateFinalPosition()
{
	POINT* pt;
	STATE* sl;
	double dt = getTimeStepSize();

	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
	     it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sl = (STATE*)left_state(pt);
            for (int j = 0; j < 3; ++j)
            {
                Coords(pt)[j] = sl->x_old[j]+sl->avgVel[j]*dt;
                if (std::isnan(Coords(pt)[j]))
                    printf("nan coords, x_old = %f, avgVel = %f\n",
                            sl->x_old[j],sl->avgVel[j]);
            }
        }
	}
}

void CollisionSolver::reduceSuperelast()
{
	bool has_superelas = true;
	int niter = 0;
    int num_edges;
	const int max_iter = 3;
	while(has_superelas && niter++ < max_iter)
    {
	    has_superelas = reduceSuperelastOnce(num_edges);
	}

	if (debugging("collision"))
        printf("    %d edges are over strain limit after %d iterations\n",num_edges,niter);
}

//TODO: This is not the correct update.
void CollisionSolver::updateFinalVelocity()
{
    //avgVel is actually the velocity at t(n+1/2)
    //need to call spring solver to get velocity at t(n+1)
    //for simplicity now set v(n+1) = v(n+1/2)


    //detectProximity();
	//detectCollision(); 
	
    POINT* pt;
	STATE* sl;
	double dt = getTimeStepSize();

	for (std::vector<CD_HSE*>::iterator it = hseList.begin();
             it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sl = (STATE*)left_state(pt);
            if (!sl->has_collsn) 
                continue;
            
            for (int j = 0; j < 3; ++j)
            {
                pt->vel[j] = sl->avgVel[j];
                sl->vel[j] = sl->avgVel[j];
                
                if (std::isnan(pt->vel[j]))
                    printf("nan vel and avgVel\n");
            }
        }
        
    }
    
    updateFinalForRG();
}

void CollisionSolver::updateFinalForRG()
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

void CollisionSolver::updateAverageVelocity()
{
	POINT *p;
	STATE *sl;
	double maxSpeed = 0;
	double* maxVel = nullptr;

#ifdef HAVE_VTK
	if (debugging("CollisionImpulse"))
    {
        char fname[200] = "vtk_test";
        static int count = 0;
        updateFinalPosition();
        if (create_directory(fname,NO))
        {
            sprintf(fname,"%s/surf-%03d.vtp",fname,count++);
            vtkplotVectorSurface(hseList,fname);
        }
	}
#endif

    unsortHseList(hseList);
	for (unsigned i = 0; i < hseList.size(); ++i)
	{
	    CD_HSE* hse = hseList[i];
	    int np = hse->num_pts(); 

	    for (int j = 0; j < np; ++j)
	    {
		p = hse->Point_of_hse(j);
		
        if (isStaticRigidBody(p)) continue;
        if (sorted(p)) continue;

		sl = (STATE*)left_state(p);
		if (sl->collsn_num > 0)
		{
		    sl->has_collsn = true;
		    for (int k = 0; k < 3; ++k)
		    {
			sl->avgVel[k] += (sl->collsnImpulse[k] + sl->friction[k])/sl->collsn_num;
			if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
			{
			    printf("inf/nan vel[%d]: impulse = %f, friction = %f, collsn_num = %d\n",
				k,sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
			    clean_up(ERROR);
			}
			//reset impulse and fricition to 0
			//collision handling will recalculate the impulse
			sl->collsnImpulse[k] = sl->friction[k] = 0.0;
		    }
		    sl->collsn_num = 0;
		}

		/* test for RG */
		if (sl->collsn_num_RG > 0)
		{
		    sl->has_collsn = true;
		    for (int k = 0; k < 3; ++k)
			sl->avgVel[k] += sl->collsnImpulse_RG[k]/sl->collsn_num_RG;
		    sl->collsn_num_RG = 0;
		}

		if (debugging("collision"))
        {
		    //debugging: print largest speed
		    double speed = Mag3d(sl->avgVel);
		    if (speed > maxSpeed) {
			maxVel = sl->avgVel;
                        maxSpeed = speed;
                    }
		}

		sorted(p) = YES;
	    }
	}
	
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects

	if (debugging("collision"))
	if (maxVel != nullptr)
	    printf("    max velocity = [%f %f %f]\n",maxVel[0],maxVel[1],maxVel[2]);
	if (debugging("collision"))
	    printDebugVariable();
}

bool CollisionSolver::isCollision(const CD_HSE* a, const CD_HSE* b){
	const CD_BOND *cd_b1, *cd_b2;
	const CD_TRI  *cd_t1, *cd_t2;
        //Commented code turns off collision detection involving the lines/strings
        /*
        if (a->name == "lines" && b->name == "lines") {
            return false;
        }
        if (a->name == "lines" && b->name == "tris_rigid" || 
            a->name == "tris_rigid" && b->name == "lines")
            return false;
       */
	double h = CollisionSolver3d::getRoundingTolerance();
	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
            return false;
	    return MovingTriToTri(t1,t2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
        //return false;
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return MovingBondToBond(b1,b2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
        //return false;
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1,h);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
        //return false;
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1,h);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
	return false;
}

//This is checking the geometric primitive for intersection
bool CollisionSolver::isProximity(const CD_HSE* a, const CD_HSE* b){
	const CD_BOND *cd_b1, *cd_b2;
	const CD_TRI  *cd_t1, *cd_t2;
	double h = CollisionSolver3d::getFabricThickness();

	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
            return false;
        /*
            if (Gindex(Point_of_tri(t1)[0]) == 18414 &&
                Gindex(Point_of_tri(t1)[1]) == 18415 &&
                Gindex(Point_of_tri(t1)[2]) == 12191 &&
                Gindex(Point_of_tri(t2)[0]) == 17741 &&
                Gindex(Point_of_tri(t2)[1]) == 10421 &&
                Gindex(Point_of_tri(t2)[2]) == 19335)
            {
                printf("t1:\n");
                print_tri_coords(t1);
                printf("t2:\n");
                print_tri_coords(t2);
            }
        */
	    return TriToTri(t1,t2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return BondToBond(b1,b2,h);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1,h);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1,h);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
	return false;
}

void CollisionSolver::printDebugVariable(){
	std::cout << "Enter EdgeToEdge " << edg_to_edg 
		  << " times"<< std::endl;
	std::cout << "Enter PointToTri " << pt_to_tri 
		  << " times"<< std::endl;
	std::cout << "Enter isCoplanar " << is_coplanar
		  << " times"<< std::endl;
	moving_edg_to_edg = moving_pt_to_tri = is_coplanar = 0;
	edg_to_edg = pt_to_tri = 0;
}

/********************************
* implementation for CD_HSE     *
*********************************/
double CD_BOND::max_static_coord(int dim){
    return std::max(Coords(m_bond->start)[dim],
		    Coords(m_bond->end)[dim]);
}

double CD_BOND::min_static_coord(int dim){
    return std::min(Coords(m_bond->start)[dim],
		    Coords(m_bond->end)[dim]);
}

double CD_BOND::max_moving_coord(int dim,double dt){
    double ans = -HUGE;
    for (int i = 0; i < 2; ++i){
	POINT* pt = (i == 0)? m_bond->start : m_bond->end;
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(ans,sl->x_old[dim]);
	ans = std::max(ans,sl->x_old[dim]+sl->avgVel[dim]*dt); 
    }    
    return ans;
}

double CD_BOND::min_moving_coord(int dim,double dt){
    double ans = HUGE;
    for (int i = 0; i < 2; ++i){
	POINT* pt = (i == 0)? m_bond->start : m_bond->end;
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(ans,sl->x_old[dim]);
	ans = std::min(ans,sl->x_old[dim]+sl->avgVel[dim]*dt); 
    }    
    return ans;
}

POINT* CD_BOND::Point_of_hse(int i) const{
    if (i >= num_pts())
	return NULL;
    else
        return (i == 0) ? m_bond->start : 
			  m_bond->end;
}

double CD_TRI::max_static_coord(int dim){
    double ans = -HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(sl->x_old[dim],ans);
    }
    return ans;
}

double CD_TRI::min_static_coord(int dim){
    double ans = HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(sl->x_old[dim],ans);
    }
    return ans;
}

double CD_TRI::max_moving_coord(int dim,double dt){
    double ans = -HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(ans,sl->x_old[dim]);
	ans = std::max(ans,sl->x_old[dim]+sl->avgVel[dim]*dt);
    }
    return ans;
}

double CD_TRI::min_moving_coord(int dim,double dt){
    double ans = HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(ans,sl->x_old[dim]);
	ans = std::min(ans,sl->x_old[dim]+sl->avgVel[dim]*dt);
    }
    return ans;
}

POINT* CD_TRI::Point_of_hse(int i) const{
    if (i >= num_pts())
	return NULL;
    else
        return Point_of_tri(m_tri)[i];
}

/*******************************
* utility functions start here *
*******************************/
/* The followings are helper functions for vector operations. */
void Pts2Vec(const POINT* p1, const POINT* p2, double* v){
	for (int i = 0; i < 3; ++i)	
	    v[i] = Coords(p1)[i] - Coords(p2)[i];
}

double distBetweenCoords(double* v1, double* v2)
{
	double dist = 0.0;
	for (int i = 0; i < 3; ++i){
		dist += sqr(v1[i]-v2[i]);
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

void scalarMult(double a,double* v, double* ans)
{
	for (int i = 0; i < 3; ++i)
            ans[i] = a*v[i];	
}

extern double myDet3d(double a[3][3]){
    return  a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) 
	  - a[0][1]*(a[1][0]*a[2][2] - a[2][0]*a[1][2]) 
	  + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]);
}

void unsortHseList(std::vector<CD_HSE*>& hseList){
	for (unsigned j = 0; j < hseList.size(); ++j)
	{
	    CD_HSE* hse = hseList[j];
	    int np = hse->num_pts();
	    for (int i = 0; i < np; ++i){
		sorted(hse->Point_of_hse(i)) = NO;
	    }
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
