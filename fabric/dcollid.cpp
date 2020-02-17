#include <FronTier.h>
#include "collid.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cfenv>

#include <omp.h>

/*****declaration of static functions starts here********/
//static void makeSet(std::vector<CD_HSE*>&);
static POINT* findSet(POINT*);
static void mergePoint(POINT*,POINT*);
inline POINT*& root(POINT*);
inline POINT*& tail(POINT*);
/*******************end of declaration*******************/

//define default parameters for collision detection
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


double CollisionSolver3d::setVolumeDiff(double vd){vol_diff = vd;}

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

//TODO: fix this
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

    //Reset points to original positions.
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

//this function is needed if collisions still
//present after several iterations;
void CollisionSolver3d::computeImpactZone()
{
    std::cout<<"Starting compute Impact Zone: "<<std::endl;

	int niter = 0;
	bool is_collision = true;

	turnOnImpZone();
    
    while(is_collision)
    {
        is_collision = false;

        //start UF alogrithm
        //merge four pts if collision happens

        start_clock("dynamic_AABB_collision");
        aabbCollision();
        abt_collision->query();
        stop_clock("dynamic_AABB_collision");

        is_collision = abt_collision->getCollsnState();

        updateImpactZoneVelocity();

        if (debugging("collision"))
        {
            std::cout << "    #"<<niter++ << ": "
                      << abt_collision->getCount() 
                      << " collision pairs" << std::endl;
            std::cout << "     " << numImpactZones
                      << " impact zones" << std::endl;
            std::cout << "     " << numImpactZonePoints
                      << " total impact zone points" << std::endl;
        }
    }
	
    turnOffImpZone();
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

void CollisionSolver3d::updateImpactZoneVelocity()
{
	numImpactZonePoints = 0;
	numImpactZones = 0;

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

	start_clock("computeAverageVelocity");
	computeAverageVelocity();
	stop_clock("computeAverageVelocity");

    start_clock("detectProximity");
	detectProximity();
	stop_clock("detectProximity");

	//if (debugging("collision"))
	  //  printDebugVariable();

	//check linear trajectories for collisions
	start_clock("detectCollision");
	detectCollision();
	stop_clock("detectCollision");

	//if (debugging("collision"))
	  //  printDebugVariable();
	
	//start_clock("detectDomainBoundaryCollision");
	//detectDomainBoundaryCollision();
	//stop_clock("detectDomainBoundaryCollision");

	//update position using final midstep velocity
	updateFinalPosition();
    //TODO: verify this should not be called a second time
            //detectProximity();
    
    //TODO: implement this function correctly
	//start_clock("reduceSuperelast");
	    //reduceSuperelast();
	//stop_clock("reduceSuperelast");
	
    //TODO: implement this function correctly
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
            if ((*it)->type == CD_HSE_TYPE::FABRIC_BOND)
                tol = CollisionSolver3d::getStringThickness();

            AABB* ab = new AABB(tol,*it);
            abt_proximity->addAABB(ab);
        }
        abt_proximity->updatePointMap(hseList);
        volume = abt_proximity->getVolume();
    }
    else
    {
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

	updateAverageVelocity();
}

// AABB tree for collision detection process
void CollisionSolver3d::aabbCollision()
{
    if (!abt_collision)
    {
        abt_collision =
            std::unique_ptr<AABBTree>(new AABBTree(MotionState::MOVING));

        for (auto it = hseList.begin(); it != hseList.end(); it++)
        {
            double tol = CollisionSolver3d::getFabricRoundingTolerance();
            if ((*it)->type == CD_HSE_TYPE::FABRIC_BOND)
                tol = CollisionSolver3d::getStringRoundingTolerance();

            AABB* ab = new AABB(tol,*it,s_dt);
            abt_collision->addAABB(ab);
        }
        abt_collision->updatePointMap(hseList);
        volume = abt_collision->getVolume();
    }
    else
    {
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
                << " collision pairs" << std::endl;
        }

        //TODO: don't call when using gauss-seidel update
        //updateAverageVelocity();
	    
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

void CollisionSolver3d::reduceSuperelast()
{
	int niter = 0;
    int num_edges;
	const int max_iter = 3;

	bool has_superelas = true;
    while(has_superelas && niter++ < max_iter)
    {
	    has_superelas = reduceSuperelastOnce(num_edges);
	}

	if (debugging("strain_limiting"))
    {
        printf("   %d edges are over strain limit\n",num_edges);
    }
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

            /*
            double x_cand[2][3];    
            for (int k = 0; k < 2; ++k)
            {
                double tmp[3];
                scalarMult(dt,sl[k]->avgVel,tmp);
                addVec(sl[k]->x_old,tmp,x_cand[k]);
            }

            //double len_new = distance_between_positions(x_cand[0],x_cand[1],3);
            */
            
            double edge_vec[3];
            minusVec(sl[1]->x_old,sl[0]->x_old,edge_vec);
            
            //double v_rel[3];
	        //minusVec(sl[0]->avgVel,sl[1]->avgVel,v_rel);
		    
            double len_old = Mag3d(edge_vec);
            //double len_old = distance_between_positions(sl[1]->x_old,sl[0]->x_old,3);
		
            double k;
		    double len0;

            if (CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(hse))
            {
                k = getFabricSpringConstant();
                len0 = cd_tri->m_tri->side_length0[j];
            }
            else if (CD_BOND* cd_bond = dynamic_cast<CD_BOND*>(hse))
            {
                k = getStringSpringConstant();
                len0 = cd_bond->m_bond->length0;
            }
            else
            {
                std::cout << "Unknown type" << std::endl;
                clean_up(ERROR);
            }
            
            //if (len_old > ROUND_EPS && len_new > ROUND_EPS)
            if (len_old > ROUND_EPS)
            {
                scalarMult(1/len_old,edge_vec,edge_vec); //normalize
                double strain = (len_old-len0)/len0;
                //double strain_rate = (len_new-len_old)/len_old;
                
                //if (fabs(strain) > superelasTol || fabs(strain_rate) > superelasTol)
                if (fabs(strain) > superelasTol)
                {
                    double correction = (fabs(strain)-superelasTol)*len0;
                    double impulse = k*correction*dt;

                    double sign = -1.0;
                    if (strain < 0.0)
                        sign = 1.0;

                    //double v_tmp1[3]; double v_tmp2[3];
                    for (int k = 0; k < 3; ++k)
                    {
                        sl[0]->avgVel[k] += sign*impulse*edge_vec[k];
                        sl[1]->avgVel[k] += -sign*impulse*edge_vec[k];
                        //v_tmp1[i] = sign*impulse*edge_vec[i];
                        //v_tmp2[i] = -sign*impulse*edge_vec[i];
                    }
                    //memcpy((void*)sl[0]->avgVel,(void*)v_tmp1,3*sizeof(double)); 
                    //memcpy((void*)sl[1]->avgVel,(void*)v_tmp2,3*sizeof(double));
                    
                    has_superelas = true;
                    num_edges++;
                }
            }
            else
            {
                printf("ERROR: len_old < ROUND_EPS\n");
                clean_up(ERROR);
                /*
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
                */
            }	
	    
        }
	
    }

	return has_superelas;
}

void CollisionSolver3d::updateFinalPosition()
{
	double dt = getTimeStepSize();
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
    //avgVel is actually the velocity at t(n+1/2)
    //need to call spring solver to get velocity at t(n+1)
    //for simplicity now set v(n+1) = v(n+1/2)
    
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
            if (!sl->has_collsn) //TODO: why is this needed?
                continue;

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

    unsortHseList(hseList);
	for (unsigned i = 0; i < hseList.size(); ++i)
	{
	    CD_HSE* hse = hseList[i];
	    int np = hse->num_pts(); 

	    for (int j = 0; j < np; ++j)
	    {
		p = hse->Point_of_hse(j);
		
        if (sorted(p) || isStaticRigidBody(p))
            continue;

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
	
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects

	if (debugging("average_velocity"))
    {
	    if (maxVel != nullptr)
	        printf("    max velocity = [%f %f %f]\n",
                    maxVel[0],maxVel[1],maxVel[2]);
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
