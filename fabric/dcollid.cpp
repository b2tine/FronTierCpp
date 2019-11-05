#include <FronTier.h>
#include "collid.h"
#include "Proximity.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cfenv>

#include <omp.h>

//TODO: Make this a function after putting VTK
//      into configure.ac and build.

/*
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
*/

/*****declaration of static functions starts here********/
//static void makeSet(std::vector<CD_HSE*>&);
//static POINT* findSet(POINT*);
//static void mergePoint(POINT*,POINT*);
//inline POINT*& root(POINT*);
//inline POINT*& tail(POINT*);
/*******************end of declaration*******************/

//define default parameters for collision detection
double CollisionSolver3d::s_eps = EPS;
double CollisionSolver3d::s_thickness = 0.001;
double CollisionSolver3d::s_dt = DT;
double CollisionSolver3d::s_k = 1000;
double CollisionSolver3d::s_m = 0.01;
double CollisionSolver3d::s_mu = 0.0;
double CollisionSolver3d::s_cr = 1.0;
bool   CollisionSolver3d::s_detImpZone = false;

//debugging variables
int CollisionSolver3d::moving_edg_to_edg = 0;
int CollisionSolver3d::moving_pt_to_tri = 0;
int CollisionSolver3d::is_coplanar = 0;
int CollisionSolver3d::edg_to_edg = 0;
int CollisionSolver3d::pt_to_tri = 0;

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

void CollisionSolver3d::setRoundingTolerance(double neweps) {s_eps = neweps;}
double CollisionSolver3d::getRoundingTolerance() {return s_eps;}

void CollisionSolver3d::setFabricThickness(double h){s_thickness = h;}
double CollisionSolver3d::getFabricThickness() {return s_thickness;}

double CollisionSolver3d::setVolumeDiff(double vd) {vol_diff = vd;}

//this function should be called at every time step
void CollisionSolver3d::setTimeStepSize(double new_dt) {s_dt = new_dt;}
double CollisionSolver3d::getTimeStepSize() {return s_dt;}

void   CollisionSolver3d::setSpringConstant(double new_k) {s_k = new_k;}
double CollisionSolver3d::getSpringConstant() {return s_k;}

//the spring model static friction coefficent
void   CollisionSolver3d::setFrictionConstant(double new_mu) {s_mu = new_mu;}
double CollisionSolver3d::getFrictionConstant() {return s_mu;}

void   CollisionSolver3d::setPointMass(double new_m) {s_m = new_m;}
double CollisionSolver3d::getPointMass() {return s_m;}

//set restitution coefficient between rigid bodies
void   CollisionSolver3d::setRestitutionCoef(double new_cr) {s_cr = new_cr;}
double CollisionSolver3d::getRestitutionCoef() {return s_cr;}


void CollisionSolver3d::assembleFromInterface(const INTERFACE* intfc, double dt)
{
	//assemble tris list from input intfc
	//this function should be called before
	//spring interior dynamics computed
	SURFACE** s;
	CURVE** c;
	TRI *tri;
	BOND *b;
	int n_tri = 0, n_bond = 0;
	setTimeStepSize(dt);
	clearHseList();
	
    intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    unsort_surface_point(*s);
	    surf_tri_loop(*s,tri)
	    {
                if (wave_type(*s) == MOVABLE_BODY_BOUNDARY || 
                    wave_type(*s) == NEUMANN_BOUNDARY)
	            hseList.push_back(new CD_TRI(tri, "tris_rigid"));
                else 
                    hseList.push_back(new CD_TRI(tri, "tris"));
		n_tri++;
	    }
	}

	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue; 
	    curve_bond_loop(*c,b)
	    {
            hseList.push_back(new CD_BOND(b,m_dim, "lines"));
		    n_bond++;
	    }
	}

	makeSet(hseList);
	createImpZoneForRG(intfc);
	setDomainBoundary(intfc->table->rect_grid.L, intfc->table->rect_grid.U);

	if (debugging("assembleFromInterface")){
	    printf("%d num of tris, %d num of bonds\n",n_tri,n_bond);
	    printf("%lu number of elements is assembled\n",hseList.size());
	}
}

void unsort_surface_point(SURFACE *surf)
{
    TRI *tri;
    POINT *p;
    int i;

    for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
                    tri = tri->next)
    {
        for (i = 0; i < 3; ++i)
        {
            p = Point_of_tri(tri)[i];
            sorted(p) = NO;
        }
    }
}       /* end unsort_surface_point */


void CollisionSolver3d::recordOriginalPosition()
{
	POINT* pt;
	STATE* sl;

    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
		    pt = (*it)->Point_of_hse(i);
		    sl = (STATE*)left_state(pt); 
		    
            sl->has_collsn = false;
		
            //TODO: why?
            if (isMovableRigidBody(pt))
                continue;
		
            for (int j = 0; j < 3; ++j)
                sl->x_old[j] = Coords(pt)[j];
	
            if (std::isnan(sl->x_old[0]))
                std::cout<<"nan_x_old"<<std::endl;
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

void CollisionSolver3d::computeAverageVelocity()
{
    POINT* pt;
    STATE* sl; 
    double dt = getTimeStepSize();
    double max_speed = 0.0;
    double* max_vel = nullptr;
    POINT* max_pt=nullptr;

    for (std::vector<CD_HSE*>::iterator it = hseList.begin();
            it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sl = (STATE*)left_state(pt); 
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
                    std::cout<<"nan avgVel" << std::endl;
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
            sl = (STATE*)left_state(max_pt);
            printf("x_old = [%f %f %f]\n",
                    sl->x_old[0],sl->x_old[1],sl->x_old[2]);
	        printf("x_new = [%f %f %f]\n",
                    Coords(max_pt)[0],Coords(max_pt)[1],Coords(max_pt)[2]);
	        printf("dt = %f\n",dt);
            printf("Gindex(max_pt) = %d\n",Gindex(max_pt));
        }
    }

    //x_old is the only valid coords for each point 
    //Coords(point) is for temporary judgement
    resetPositionCoordinates();
}

void CollisionSolver3d::resetPositionCoordinates()
{
    POINT* pt;
    STATE* sl; 

    for (std::vector<CD_HSE*>::iterator it = hseList.begin();
            it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sl = (STATE*)left_state(pt);
            for (int j = 0; j < m_dim; ++j)
                Coords(pt)[j] = sl->x_old[j];
        }
    }
}

void CollisionSolver3d::turnOffImpZone()
{
    s_detImpZone = false;
}

void CollisionSolver3d::turnOnImpZone(){
    s_detImpZone = true;
}

bool CollisionSolver3d::getImpZoneStatus()
{
    return s_detImpZone;
}

void CollisionSolver3d::resolveCollision()
{
	//catch floating point exception: nan/inf
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	start_clock("computeAverageVelocity");
	computeAverageVelocity();
	stop_clock("computeAverageVelocity");
	
    start_clock("detectProximity");
	detectProximity();
	stop_clock("detectProximity");

	//if (debugging("printDebugVariable"))
	  //  printDebugVariable();

	//check linear trajectories for collisions
	start_clock("detectCollision");
	detectCollision();
	stop_clock("detectCollision");

	//if (debugging("printDebugVariable"))
	  //  printDebugVariable();
	
	//start_clock("detectDomainBoundaryCollision");
	//detectDomainBoundaryCollision();
	//stop_clock("detectDomainBoundaryCollision");

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

void CollisionSolver3d::detectProximity()
{
    aabbProximity();
    proximityCandidates.clear();
    proximityCandidates = abt_collision->getCandidates();

	if (debugging("proximity"))
        std::cout << candidates.size()
            << " pair of proximity candidates" << std::endl;

    processProximityCandidates();
}

//TODO: Only build proximity tree once at startup
//      using hilbert curves for bottom up construction.
//      Will need to use global triangle and bond indices,
//      and preserve restart functionality.
//
//Build AABB tree for proximity detection process
void CollisionSolver3d::aabbProximity()
{
    if (!abt_proximity)
    {
        double pre_tol = CollisionSolver3d::getFabricThickness();
        abt_proximity =
            std::unique_ptr<AABBTree>(new AABBTree(MotionState::STATIC));

        for (auto it = hseList.begin(); it != hseList.end(); it++)
        {
             AABB* ab = new AABB(pre_tol,*it,abt_proximity->getType());
             abt_proximity->addAABB(ab);
        }
        abt_proximity->updatePointMap(hseList);
        proximity_vol = abt_proximity->getVolume();
    }
}

void CollisionSolver3d::processProximityCandidates()
{
    Proximities.clear();

    std::vector<NodePair>::iterator nit;
    for (nit = proximityCandidates.begin(); nit < proximityCandidates.end(); ++nit)
    {
        Node* A = nit->first;
        Node* B = nit->second;
        
        CD_HSE* a = A->data->hse;
        CD_HSE* b = B->data->hse;

        std::unique_ptr<Proximity> proximity = checkProximity(a,b,s_thickness);
        if (proximity)
        {
            proximity->computeImpulse();
            Proximities.push_back(std::move(proximity));
        }
    }

    std::vector<unique_ptr<Proximity>>::iterator pit;
    for (pit = Proximities.begin(); pit < Proximities.end(); ++pit)
    {
        (*pit)->updateAverageVelocity();
    }
}

//TODO: use gauss-seidel updates
//TODO: finish updating for new data structures
void CollisionSolver3d::detectCollision()
{
    if (debugging("collision"))
        std::cout<<"Starting collision handling: "<<std::endl;
	
	const int MAX_ITER = 12;
    const double h = CollisionSolver3d::getRoundingTolerance();
	
    int niter = 1;
	int cd_count = 0;
   
    //TODO: keep track of total elapsed time?
    
    bool is_collision = true; 

    while(is_collision)
    {
	    is_collision = false;
	    
        //TODO: process one collision pair, then recompute
        //      tree and get new candidates?
        aabbCollision();
        collisionCandidates.clear();
        collisionCandidates = abt_collision->getCandidates();

	    if (debugging("collision"))
            std::cout<<"    #"<<niter << ": " << collisionCandidates.getCount() 
                << " pair of collision candidates" << std::endl;

        processCollisionCandidates();

        if (++niter > MAX_ITER)
            break;
	}

    //TODO: implement computeImpactZone() using new data structures
    //      SEE PROVOT PAPER
	if (is_collision) 
    {
        start_clock("computeImpactZone");
	    computeImpactZone();
        stop_clock("computeImpactZone");
    }
}

//Build/Update AABB tree for collision detection process
void CollisionSolver3d::aabbCollision()
{
    if (!abt_collision)
    {
        //double pre_tol = CollisionSolver3d::getRoundingTolerance();
        double pre_tol = CollisionSolver3d::getFabricThickness();
        abt_collision =
            std::unique_ptr<AABBTree>(new AABBTree(MotionState::MOVING));

        for (auto it = hseList.begin(); it != hseList.end(); it++)
        {
             AABB* ab = new AABB(pre_tol,*it,abt_collision->getType(), s_dt);
             abt_collision->addAABB(ab);
        }
        abt_collision->updatePointMap(hseList);
        collision_vol = abt_collision->getVolume();
    }
    else
    {
        abt_collision->setTimeStep(s_dt);
        abt_collision->updateAABBTree(hseList);
        abt_collision->updateTreeStructure();
        collision_vol = abt_collision->getVolume();
    }
}

//For sorting Collisions by time of occcurence.
//Still unclear if we want to do this.
static bool CollisionCompare(
        std::unique_ptr<Collision>& A,
        std::unique_ptr<Collision>& B)
{
    return A->dt < B->dt;
}

//TODO: iteratively process collision candidates in
//      Gauess-Seidel fashion
void CollisionSolver3d::processCollisionCandidates()
{
    Collisions.clear();

    std::vector<NodePair>::iterator it;
    for (it = collisionCandidates.begin(); it < collisionCandidates.end(); ++it)
    {
        Node* A = it->first;
        Node* B = it->second;

        CD_HSE* a = A->data->hse;
        CD_HSE* b = B->data->hse;

        std::unique_ptr<Collision> collsn = checkCollision(a,b,s_eps);
        if (collsn)
        {
            collsn->computeImpulse();
            collsn->updateState();
            
            //TODO: check "final position" for proximity
            
            //save so average velocity can be reset if neccessary
            Collisions.push_back(std::move(collsn));
        }
    }

    /*
    //Sort the Collisions vector by time of collision
    std::sort(Collisions.begin(),Collisions.end(),CollisionCompare);

    std::vector<std::unique_ptr<Collision>>::iterator it;
    for (it = Collisions.begin(); it < Collisions.end(); ++it)
    {
        (*it)->computeImpulse();
        (*it)->updateAverageVelocity();
        //TODO: check "final position" for proximity
    }
    */
}

//TODO: rewrite this with new data structures
void CollisionSolver3d::computeImpactZone()
{
    std::cout<<"Starting compute Impact Zone: "<<std::endl;

    const double h = CollisionSolver3d::getRoundingTolerance();

	int niter = 0;
    int numZones = 0;
	bool is_collision = true;

	turnOnImpZone();
	//makeSet(hseList); //this is done in AABBTree::updatePointMap()
    
    //int impzone_counter = 0
    //TODO: This can enter infinite loop
    while(is_collision)
    {
        is_collision = false;

        //start UF alogrithm
        //merge four pts if collision happens

        start_clock("dynamic_AABB_collision");
        aabbCollision();
        abt_collision->query(h);
        stop_clock("dynamic_AABB_collision");

        is_collision = abt_collision->getCollsnState();

        //TODO: Verify Jarret's claim that this should be removed.
        //updateAverageVelocity();

        updateImpactZoneVelocity(numZones);

        std::cout <<"    #"<<niter++ << ": " << abt_collision->getCount() 
                  << " pair of collision tris" << std::endl;
        std::cout <<"     "<< numZones
                  <<" zones of impact" << std::endl;
    }
	
    turnOffImpZone();
}

//TODO: See aftest.cpp airfoil_stat functions which record
//      the strain in a spring mass surface (potential energy)
//      Then look at how to apply strain limiting impulses
void CollisionSolver3d::reduceSuperelast()
{
	bool has_superelas = true;
	int niter = 0;
    int num_edges;
	const int max_iter = 3;
	while(has_superelas && niter++ < max_iter)
    {
	    has_superelas = reduceSuperelastOnce(num_edges);
	}

	if (debugging("strain_limit"))
        printf("    %d edges are over strain limit after %d iterations\n",num_edges,niter);
}

//TODO: Implement this correctly.
//      jacobi iteration style for strain and
//      gauss-seidel iteration style for strain rate.
//      Should be called after collisions have been handled.
bool CollisionSolver3d::reduceSuperelastOnce(int& num_edges)
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

//TODO: This does not work. Fix it.
void CollisionSolver3d::detectDomainBoundaryCollision() {
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
        //TODO: Don't think this is actually in the tangential direction
        //
		//reduce tangential velocity with friction
		double preVt = Mag3d(sl->avgVel);
		if (preVt > MACH_EPS)
		for (int j = 0; j < m_dim; ++j) 
		    sl->avgVel[j] *= std::max(1.0-mu*dv/preVt,0.0); 
	    }
	}
}

void CollisionSolver3d::updateFinalPosition()
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

            std::vector<double> x_old(sl->x_old,sl->x_old+3);
            std::vector<double> avgVel(sl->avgVel,sl->avgVel+3);

            for (int j = 0; j < 3; ++j)
            {
                Coords(pt)[j] = x_old[j] + avgVel[j]*dt;
                sl->x_old[j] = x_old[j] + avgVel[j]*dt;
            }

            /*
            if (std::isnan(Mag3d(Coords(pt))))
            {
                for (int i = 0; i < 3; ++i)
                    printf("nan coords, x_old = %f, avgVel = %f\n",
                            x_old[j],avgVel[j]);
                clean_up(ERROR);
            }
            */
        }
	}
}

//TODO: This is not the correct update.
//      Could try using generic_spring_solver()
//      which is 4th order kunga kutta.
void CollisionSolver3d::updateFinalVelocity()
{
    //avgVel is actually the velocity at t(n+1/2)
    //need to call spring solver to get velocity at t(n+1)
    //for simplicity now set v(n+1) = v(n+1/2)

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
            
            /*
            if (!sl->has_collsn) 
                continue;
            */

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


//TODO: Nothing related to the collision impulses
//      should be in this function.
//      Proximity and Friction impulses only.
//
//      Collision impulses may have a similar function,
//      but they will not consider friction, and will
//      update the "final midstep" velocity instead.
//      
//      We may be able to reuse for the Collision impulses
//      but is better to have a conceptual barrier between
//      the two for the moment.

/*
void CollisionSolver3d::updateAverageVelocity()
{
	POINT *p;
	STATE *sl;
	double maxSpeed = 0;
	double* maxVel = nullptr;

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
                sl->avgVel[k] += sl->collsnImpulse[k]/sl->collsn_num;
                //TODO: does friction impulse also need to divide by num collision?
                sl->avgVel[k] += sl->friction[k]/sl->collsn_num;
                //sl->avgVel[k] += sl->friction[k];
                
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

		// test for RG
		if (sl->collsn_num_RG > 0)
		{
		    sl->has_collsn = true;
		    for (int k = 0; k < 3; ++k)
			sl->avgVel[k] += sl->collsnImpulse_RG[k]/sl->collsn_num_RG;
		    sl->collsn_num_RG = 0;
		}

		if (debugging("average_velocity"))
        {
		    double speed = Mag3d(sl->avgVel);
		    if (speed > maxSpeed) {
			maxVel = sl->avgVel;
                        maxSpeed = speed;
                    }
		}

		sorted(p) = YES;
	    }
	}
	
    //TODO: Consider moving this, and all other non-fabric code
    //      into seperate functions/managers.
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects

	if (debugging("average_velocity"))
	if (maxVel != nullptr)
	    printf("\t\tmax velocity = [%f %f %f]\n",maxVel[0],maxVel[1],maxVel[2]);
	if (debugging("printDebugVariable"))
	    printDebugVariable();
}
*/

void CollisionSolver3d::printDebugVariable(){
	std::cout << "Enter EdgeToEdge " << edg_to_edg 
		  << " times"<< std::endl;
	std::cout << "Enter PointToTri " << pt_to_tri 
		  << " times"<< std::endl;
	std::cout << "Enter isCoplanar " << is_coplanar
		  << " times"<< std::endl;
	moving_edg_to_edg = moving_pt_to_tri = is_coplanar = 0;
	edg_to_edg = pt_to_tri = 0;
}

void unsortHseList(std::vector<CD_HSE*>& hseList)
{
	for (unsigned j = 0; j < hseList.size(); ++j)
	{
	    CD_HSE* hse = hseList[j];
	    int np = hse->num_pts();
	    for (int i = 0; i < np; ++i){
		sorted(hse->Point_of_hse(i)) = NO;
	    }
	}
}

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


/*******************************
* utility functions start here *
*******************************/

/* The followings are helper functions for vector operations. */
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

void Pts2Vec(const POINT* p1, const POINT* p2, double* v){
	for (int i = 0; i < 3; ++i)	
	    v[i] = Coords(p2)[i] - Coords(p1)[i];
}

double distBetweenCoords(double* v1, double* v2)
{
	double dist = 0.0;
	for (int i = 0; i < 3; ++i){
		dist += sqr(v1[i]-v2[i]);
	}
	return std::sqrt(dist);
}

double myDet3d(double a[3][3]){
    return  a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) 
	  - a[0][1]*(a[1][0]*a[2][2] - a[2][0]*a[1][2]) 
	  + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]);
}


