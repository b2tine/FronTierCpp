#include <armadillo>
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
bool CollisionSolver3d::gs_update = false;

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
Front* CollisionSolver3d::ft;

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
double CollisionSolver3d::getStrainLimit() {return strain_limit;}

void CollisionSolver3d::setCompressiveStrainLimit(double cslim) {compressive_strain_limit = cslim;}
double CollisionSolver3d::getCompressiveStrainLimit() {return compressive_strain_limit;}

void CollisionSolver3d::setStrainRateLimit(double srlim) {strainrate_limit = srlim;}
double CollisionSolver3d::getStrainRateLimit() {return strainrate_limit;}

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

void unsortHseList(std::vector<CD_HSE*>& list)
{
    for (auto it = list.begin(); it < list.end(); ++it)
	{
	    int np = (*it)->num_pts();
	    for (int i = 0; i < np; ++i)
            sorted((*it)->Point_of_hse(i)) = NO;
	}
}

void CollisionSolver3d::initializeSystem(Front* front)
{
    ft = front;
    setStep(front->step);
    setTimeStepSize(front->dt);
    setOutputDirectory(OutName(front));
    assembleFromInterface(front->interf);
    recordOriginalPosition();
    setHseTypeLists();
    initializeImpactZones();
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
	int n_fabric_tri = 0;
	int n_static_rigid_tri = 0;
	int n_movable_rigid_tri = 0;
    int n_bond = 0;
	
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
                n_fabric_tri++;
            }
            else if (wave_type(*s) == NEUMANN_BOUNDARY)
            {
                tag = CD_HSE_TYPE::STATIC_RIGID_TRI;
                n_static_rigid_tri++;
            }
            else if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
            {
                tag = CD_HSE_TYPE::MOVABLE_RIGID_TRI;
                n_movable_rigid_tri++;
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
        if (is_bdry(*c)) continue;
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue; 

        unsort_curve_point(*c);//TODO: can remove?

        CD_HSE_TYPE tag = CD_HSE_TYPE::STRING_BOND;
	    curve_bond_loop(*c,b)
	    {
            hseList.push_back(new CD_BOND(b,tag));
		    n_bond++;
	    }
	}

    setSizeCollisionTimes(hseList.size());
    setDomainBoundary(intfc->table->rect_grid.L,intfc->table->rect_grid.U);
	
	if (debugging("intfc_assembly"))
    {
	    printf("%d num of tris, %d num of bonds\n",n_tri,n_bond);
	    printf("%lu number of elements is assembled\n",hseList.size());
	    printf("%d num fabric tris\n",n_fabric_tri);
	    printf("%d num static rigid tris\n",n_static_rigid_tri);
	    printf("%d num movable rigid tris\n",n_movable_rigid_tri);
	}
}

/*
//NOTE: Must be called before calling the spring solver
void CollisionSolver3d::assembleFromInterface(ELASTIC_SET* geom_set)
{
	clearHseList();

	SURFACE** s;
	CURVE** c;
	TRI *tri;
	BOND *b;

	int n_tri = 0;
	int n_fabric_tri = 0;
	int n_static_rigid_tri = 0;
	int n_movable_rigid_tri = 0;
    int n_bond = 0;
	
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
                n_fabric_tri++;
            }
            else if (wave_type(*s) == NEUMANN_BOUNDARY)
            {
                tag = CD_HSE_TYPE::STATIC_RIGID_TRI;
                n_static_rigid_tri++;
            }
            else if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
            {
                tag = CD_HSE_TYPE::MOVABLE_RIGID_TRI;
                n_movable_rigid_tri++;
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
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue; 

        unsort_curve_point(*c);

        CD_HSE_TYPE tag = CD_HSE_TYPE::STRING_BOND;
	    curve_bond_loop(*c,b)
	    {
            hseList.push_back(new CD_BOND(b,tag));
		    n_bond++;
	    }
	}

    setSizeCollisionTimes(hseList.size());
    setDomainBoundary(intfc->table->rect_grid.L,intfc->table->rect_grid.U);
	
	if (debugging("intfc_assembly"))
    {
	    printf("%d num of tris, %d num of bonds\n",n_tri,n_bond);
	    printf("%lu number of elements is assembled\n",hseList.size());
	    printf("%d num fabric tris\n",n_fabric_tri);
	    printf("%d num static rigid tris\n",n_static_rigid_tri);
	    printf("%d num movable rigid tris\n",n_movable_rigid_tri);
	}
}
*/

void CollisionSolver3d::setDomainBoundary(double* L, double* U)
{
	for (int i = 0; i < m_dim; ++i)
    {
	    Boundary[i][0] = L[i];
	    Boundary[i][1] = U[i];
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
            
            sl->has_proximity = false;
            sl->has_strainlim_prox = false;
            sl->has_collsn = false;
            sl->has_strainlim_collsn = false;
            sl->collsn_dt = -1.0;

            //NOTE: sl->x_old for movable rigid body points
            //      is recorded in rgbody_point_propagate()
            if (isMovableRigidBody(pt)) continue;

            for (int j = 0; j < 3; ++j)
            {
                sl->x_old[j] = Coords(pt)[j];
            
                if (std::isnan(sl->x_old[j]))
                {
                    std::cout << "nan_x_old" << std::endl;
                    LOC(); clean_up(ERROR);
                }
            }
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

void CollisionSolver3d::setHseTypeLists()
{
    fabricTriList = getHseTypeList(CD_HSE_TYPE::FABRIC_TRI);
    staticRigidTriList = getHseTypeList(CD_HSE_TYPE::STATIC_RIGID_TRI);
    movableRigidTriList = getHseTypeList(CD_HSE_TYPE::MOVABLE_RIGID_TRI);
    stringBondList = getHseTypeList(CD_HSE_TYPE::STRING_BOND);
    elasticHseList.assign(stringBondList.begin(),stringBondList.end());
    elasticHseList.insert(elasticHseList.end(),fabricTriList.begin(),fabricTriList.end());
}

void CollisionSolver3d::initializeImpactZones()
{
    makeSet(hseList);
    initRigidBodyImpactZones();
        //initRigidBodyImpactZones(ft->interf);
}

//TODO: Can we remove?
void CollisionSolver3d::initializeImpactZones(const INTERFACE* intfc)
{
    makeSet(hseList);
    initRigidBodyImpactZones();
        //initRigidBodyImpactZones(intfc);
}

void CollisionSolver3d::resolveCollision()
{
	//catch floating point exception: nan/inf
        //feenableexcept(FE_INVALID | FE_OVERFLOW);

	computeAverageVelocity();
    resetPositionCoordinates();
            
        //computeMaxSpeed(); //debug

    // Apply impulses to enforce strain limiting
    // distance/positional constraints on adjacent mesh vertices. 
    if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
    {
        limitStrainPosnJac();
            //computeMaxSpeed(); //debug
    }

    // Static proximity handling
    start_clock("detectProximity");
	detectProximity();
    stop_clock("detectProximity");
    
    /*
    //TODO: updateImpactListVelocity() needs to be modified in order
    //      to treat movable rigid bodies as impact zones.
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects
    */

    saveAverageVelocity();
	
    // Check linear trajectories for collisions
    start_clock("detectCollision");
	detectCollision();
    stop_clock("detectCollision");

    /*
    //TODO: function needs fixing -- do we even need this?
    detectDomainBoundaryCollision();
    */

	//update position using final midstep velocity
	updateFinalPosition();

    //TODO: Don't think we need this ...
    //check for proximity again at end step positions
        //detectProximity();

    // Zero out the relative velocity between adjacent mesh vertices
    // with excess edge strain directed along their connecting edge.
    if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
    {
        limitStrainVelGS();
        //limitStrainVelJAC();
            //computeMaxSpeed(); //debug
    }

	updateFinalVelocity();

    updateFinalForRG();
}

//for debugggin -- should rename to printMaxSpeed()
void CollisionSolver3d::computeMaxSpeed()
{
    POINT* p;
    STATE* sl;
    
    double max_speed = 0.0;
    double* max_vel = nullptr;
    POINT* max_pt = nullptr;

    unsortHseList(elasticHseList);
    
    std::vector<CD_HSE*>::iterator it;
    for (it = elasticHseList.begin(); it < elasticHseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            p = (*it)->Point_of_hse(i);
            if (sorted(p)) continue;
            
            sl = (STATE*)left_state(p); 
            
            if (Mag3d(sl->avgVel) > max_speed)
            {
                max_speed = Mag3d(sl->avgVel);
                max_vel = sl->avgVel;
                max_pt = p;
            }

            sorted(p);
        }
    }


    prev_max_fabric_speed = max_fabric_speed;
    max_fabric_speed = max_speed;


    if (debugging("collision_max_speed"))
    {
        std::cout << "\n    Max fabric speed = " << max_speed << "\n";
        
        if (max_vel)
        {
            std::cout << "      Maximum average velocity is "
                << max_vel[0] << " "
                << max_vel[1] << " "
                << max_vel[2] << "\n"; 
        }
        
        if (max_pt)
        {
            STATE* sl = (STATE*)left_state(max_pt);
            printf("        x_old = [%f %f %f]\t",
                    sl->x_old[0],sl->x_old[1],sl->x_old[2]);
            printf("Gindex(max_pt) = %d\n\n",Gindex(max_pt));
        }
    }
}

void CollisionSolver3d::computeAverageVelocity()
{
    POINT* p;
    STATE* sl;
    
    double max_speed = 0.0;
    double* max_vel = nullptr;
    POINT* max_pt = nullptr;

    double dt = getTimeStepSize();

    unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
    for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            p = (*it)->Point_of_hse(i);
            if (sorted(p) || isStaticRigidBody(p)) continue;
            
            sl = (STATE*)left_state(p); 
            for (int j = 0; j < 3; ++j)
    	    {
                if (dt > ROUND_EPS)
                {
                    sl->avgVel[j] = (Coords(p)[j] - sl->x_old[j])/dt;
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
                            dt,sl->x_old[j],Coords(p)[j]);
                    LOC(); clean_up(ERROR);
                }
		
            }

            if (debugging("average_velocity"))
            {
                if (Mag3d(sl->avgVel) >= max_speed)
                {
                    max_speed = Mag3d(sl->avgVel);
                    max_vel = sl->avgVel;
                    max_pt = p;
                }
            }

            sorted(p);
        
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
}

void CollisionSolver3d::resetPositionCoordinates()
{
    POINT* p;
    STATE* sl;

    unsortHseList(hseList);

    for (auto it = hseList.begin(); it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            p = (*it)->Point_of_hse(i);
            if (sorted(p) || isStaticRigidBody(p)) continue;
            
            STATE* sl = (STATE*)left_state(p);
            for (int j = 0; j < m_dim; ++j)
                Coords(p)[j] =  sl->x_old[j];

            sorted(p);
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
        updateAverageVelocity(MotionState::STATIC);

        //computeMaxSpeed(); //debug    

    if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
    {
        limitStrainRatePosnGS(MotionState::STATIC);
        //limitStrainRatePosnJac(MotionState::STATIC);
            //computeMaxSpeed(); //debug    
    }
}

// function to perform AABB tree building, updating structure
// and query for proximity detection process
void CollisionSolver3d::aabbProximity()
{
    if (!abt_proximity)
    {
        abt_proximity =
            std::unique_ptr<AABBTree>(new AABBTree(MotionState::STATIC));

        for (auto it = hseList.begin(); it != hseList.end(); ++it)
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
	    LOC(); clean_up(ERROR);
	}
}

//For Jacobi velocity update
void CollisionSolver3d::updateAverageVelocity(MotionState mstate)
{
	POINT *p;
	STATE *sl;
	double maxSpeed = 0;
	double* maxVel = nullptr;

    unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        int np = (*it)->num_pts(); 
	    for (int j = 0; j < np; ++j)
	    {
            p = (*it)->Point_of_hse(j);
            if (sorted(p) || isStaticRigidBody(p)) continue;

            sl = (STATE*)left_state(p);
            if (sl->collsn_num > 0)
            {
                if (mstate == MotionState::STATIC)
                    sl->has_proximity = true;
                else
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

            /*
            //NOTE: This is for rigid-rigid body collision.
            if (sl->collsn_num_RG > 0)
            {
                sl->has_collsn = true;
                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] += sl->collsnImpulse_RG[k]/sl->collsn_num_RG;
                }
                sl->collsn_num_RG = 0;
            }
            */

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

void CollisionSolver3d::saveAverageVelocity()
{
	POINT *p;
	STATE *sl;

    unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        int np = (*it)->num_pts(); 
	    for (int j = 0; j < np; ++j)
	    {
            p = (*it)->Point_of_hse(j);
            if (sorted(p)) continue;

            sl = (STATE*)left_state(p);
            for (int k = 0; k < 3; ++k)
            {
                sl->avgVel_postprox[k] = sl->avgVel[k];
            }

            sorted(p) = YES;
	    }
	}
}

void CollisionSolver3d::revertAverageVelocity()
{
	POINT *p;
	STATE *sl;

    unsortHseList(hseList);
    
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
        int np = (*it)->num_pts(); 
	    for (int j = 0; j < np; ++j)
	    {
            p = (*it)->Point_of_hse(j);
            if (sorted(p)) continue;

            sl = (STATE*)left_state(p);
            for (int k = 0; k < 3; ++k)
            {
                sl->avgVel[k] = sl->avgVel_postprox[k];
                sl->has_collsn = false;
                sl->has_strainlim_collsn = false;
            }

            sorted(p) = YES;
	    }
	}

    clearCollisionTimes();
}

void CollisionSolver3d::detectCollision()
{
	std::cout << "Starting collision handling: " << std::endl;
	
	const int MAX_ITER = 8;
	    //const int MAX_ITER = 12;
	
    int niter = 0;
	int cd_count = 0;
    bool is_collision = true; 

    while(is_collision)
    {
        niter++;
	    is_collision = false;
	    
        start_clock("dynamic_AABB_collision");
        aabbCollision();
        abt_collision->turn_on_GS_update();
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
        
                //computeMaxSpeed(); //debug
        }
        
        if (is_collision)
        {
            if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
            {
                limitStrainRatePosnGS(MotionState::MOVING);
                //limitStrainRatePosnJac(MotionState::MOVING);
                    //computeMaxSpeed(); //debug
            }
        }

        if (niter >= MAX_ITER) break;
	}

    start_clock("computeImpactZone");
	if (is_collision) 
    {
        //TODO: Return avg_vel to value before point to point collisions???
        //      See todo in computeImpactZoneJac() regarding a startup step
        //
        revertAverageVelocity();
        computeImpactZoneJac();
            //computeImpactZoneGS();
    }
    stop_clock("computeImpactZone");

	std::cout << "End collision handling. " << std::endl;
}

// AABB tree for kinetic collision detection process
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

//TODO: If impact zone handling is enabled, should all
//      points of the hypersurface elements be added to
//      an impact zone all at once, rather than only adding
//      4 interfering points at a time in the point to tri
//      and edge to edge tests?
bool getCollisionGS(const CD_HSE* a, const CD_HSE* b)
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
	    return MovingTriToTriGS(t1,t2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return MovingBondToBondGS(b1,b2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		     (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBondGS(t1,b1);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
             (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBondGS(t1,b1);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    LOC(); clean_up(ERROR);
	}
}

bool getCollisionJac(const CD_HSE* a, const CD_HSE* b)
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
	    return MovingTriToTriJac(t1,t2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return MovingBondToBondJac(b1,b2);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		     (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBondJac(t1,b1);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
             (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBondJac(t1,b1);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    LOC(); clean_up(ERROR);
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
                    printf("updateFinalPosition() ERROR: nan final position\n");
                    fprint_general_vector(stdout,"x_old",sl->x_old,MAXD,"\n");
                    fprint_general_vector(stdout,"avgVel",sl->avgVel,MAXD,"\n");
                    LOC(); clean_up(EXIT_FAILURE);
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
            if (!sl->has_proximity && !sl->has_strainlim_prox &&
                !sl->has_collsn && !sl->has_strainlim_collsn) continue;

            for (int j = 0; j < 3; ++j)
            {
                pt->vel[j] = sl->avgVel[j];
                sl->vel[j] = sl->avgVel[j];
                
                if (std::isnan(pt->vel[j]))
                {
                    printf("updateFinalVelocity() ERROR: nan final velocity\n");
                    fprint_general_vector(stdout,"coords",Coords(pt),MAXD,"\n");
                    fprint_general_vector(stdout,"vel",pt->vel,MAXD,"\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }
            }
        }
        
    }
    
}

//Factored into separate position and velocity updates
//in order to apply strain limiting impulses to excessively
//strained edges in the final configuration. 
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
            if (!sl->has_proximity && !sl->has_strainlim_prox &&
                !sl->has_collsn && !sl->has_strainlim_collsn) continue;
            
            for (int j = 0; j < 3; ++j)
            {
                pt->vel[j] = sl->avgVel[j];
                sl->vel[j] = sl->avgVel[j];
            }
        }
	}
    
}

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
                    center_of_mass(pt->hs)[j] = (mrg_com[rg_index])[j] + sl->avgVel[j]*dt;
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

//Note: num has default value of 4,
//and first has default value of false
extern void createImpZone(POINT* pts[], int num, bool first)
{
	for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
            //TODO: Should it check isRigidBody() instead?
            //      Otherwise, static rigid bodies can become
            //      part of the impact zone.
            //
            //      However, this appears to be useful for resting
            //      fabric-solid contact -- should find out if
            //      anyone has done this before. Almost certain
            //      this was done unintentionally ... 
            
            if (!first)
            {
                if (isRigidBody(pts[i]) ||
                    isRigidBody(pts[j])) continue;
                /*
                if (isMovableRigidBody(pts[i]) ||
                    isMovableRigidBody(pts[j])) continue;
                */
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

//What is the minimum number of mergePoint() calls?
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
//          mergePoint(pts[0],pts[2]);
//          mergePoint(pts[0],pts[3]);
//
//      The current looping structure makes 6 calls to mergePoint().
//      Each findSet() call within mergePoint() is O(n) in the worst
//      case, i.e. when impact zones have grown large.
	
void createImpactZone(POINT* pts[], int num)
{
    for (int i = 1; i < num; ++i)
	{
        mergePoint(pts[0],pts[i]);
	}
}

//void CollisionSolver3d::initRigidBodyImpactZones(const INTERFACE* intfc)
void CollisionSolver3d::initRigidBodyImpactZones()
{
    POINT* pts[3];
	
    unsortHseList(movableRigidTriList);
	for (auto it = movableRigidTriList.begin(); it != movableRigidTriList.end(); ++it)
    {
        for (int i = 0; i < 3; ++i)
            pts[i] = (*it)->Point_of_hse(i);
            
        createImpactZone(pts,3);
    }

    /*
	SURFACE** s;
	TRI* tri;

	intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    if (!isMovableRigidBody(Point_of_tri(first_tri(*s))[0])) continue;
	        //if (!isRigidBody(Point_of_tri(first_tri(*s))[0])) continue;

        surf_tri_loop(*s, tri)
	    {
            //TODO: can make more efficient by marking points sorted
            //      and skipping if (sorted(pt)). Hold a pointer to
            //      the head that each point is always merged to
            //      mergePoint(head,p);
    		
            createImpactZone(Point_of_tri(tri),3);
    		//createImpactZoneRigidBody(Point_of_tri(tri),3);
	    }
	}
    */
}

//TODO: may not need this
void createImpactZoneRigidBody(POINT* pts[], int num)
{
	for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
            mergePoint(pts[i],pts[j]); 
	    }
	}
}

void CollisionSolver3d::turnOnImpZone(){s_detImpZone = true;}
void CollisionSolver3d::turnOffImpZone(){s_detImpZone = false;}
bool CollisionSolver3d::getImpZoneStatus(){return s_detImpZone;}

void CollisionSolver3d::computeImpactZoneGS()
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
        abt_collision->turn_on_GS_update();
        abt_collision->query();
        stop_clock("dynamic_AABB_collision");

        is_collision = abt_collision->getCollsnState();
        if (is_collision)
        {
            //Update due to connecting nearby impact zones only.
            //Gauss-seidel updates have already been applied.
            if (niter > 5)
            {
                connectNearbyImpactZones();
                updateImpactZoneVelocity();
            }
        }


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

                //computeMaxSpeed(); //debug
        }
        
        //TODO: Appropriate to limit strain rate during impact zone handling?
        if (is_collision)
        {
            if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
            {
                limitStrainRatePosnGS(MotionState::MOVING);
                //limitStrainRatePosnJac(MotionState::MOVING);
                    //computeMaxSpeed(); //debug
            }
        }

        if (niter >= MAXITER)
        {
            printf("computeImpactZoneGS(): ERROR\n\t\
                    maxiters %d without convergence!\n",
                    MAXITER);
            
            debugImpactZones();
            clean_up(EXIT_FAILURE);
        }
    }
	
    turnOffImpZone();
}

//TODO: Try reverting to post proximity, and then perform one round of ordinary
//      point to point collision (i.e. repeat first iteration from detectCollision())
//      as a startup step to impact zone handling. Then start impact zone handling.
//
//      The idea is...
void CollisionSolver3d::computeImpactZoneJac()
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
        abt_collision->turn_off_GS_update();//turn on Jacobi style update
        abt_collision->query();
        stop_clock("dynamic_AABB_collision");

        is_collision = abt_collision->getCollsnState();
        if (is_collision)
        {
            if (niter > 5)
            {
                connectNearbyImpactZones();
            }
            updateImpactZoneVelocity();
        }

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

                //computeMaxSpeed(); //debug
        }
        
        //TODO: Appropriate to limit strain rate during impact zone handling?
        if (is_collision)
        {
            if (debugging("strain_limiting")) //if (!debugging("strainlim_off"))
            {
                limitStrainRatePosnGS(MotionState::MOVING);
                //limitStrainRatePosnJac(MotionState::MOVING);
                    //computeMaxSpeed(); //debug
            }
        }

        if (niter >= MAXITER)
        {
            printf("computeImpactZoneJac(): ERROR\n\t\
                    maxiters %d without convergence!\n",
                    MAXITER);
            
            debugImpactZones();
            clean_up(EXIT_FAILURE);
        }
    }
	
    turnOffImpZone();
}

void CollisionSolver3d::connectNearbyImpactZones()
{
	unsortHseList(elasticHseList);
	for (auto it = elasticHseList.begin(); it != elasticHseList.end(); ++it)
    {
        if ((*it)->type == CD_HSE_TYPE::STRING_BOND)
        {
            POINT* p1 = (*it)->Point_of_hse(0);
            POINT* head1 = findSet(p1);
            if (weight(head1) != 1) continue;

            CD_BOND* cd_bond = dynamic_cast<CD_BOND*>(*it);
            BOND* bond = cd_bond->m_bond;
            BOND* nb_bond = bond->prev;
            if (!nb_bond) continue;
                
            POINT* p0 = nb_bond->start;
            POINT* head0 = findSet(p0);
                
            POINT* p2 = (*it)->Point_of_hse(1);
            POINT* head2 = findSet(p2);

            if (weight(head0) > 1 && weight(head2) > 1)
            {
                POINT* pts[3] = {p0, p1, p2};
                createImpactZone(pts,3);
            }
        }
        else if ((*it)->type == CD_HSE_TYPE::FABRIC_TRI)
        {
            CD_TRI* cd_tri = dynamic_cast<CD_TRI*>(*it);
            TRI* tri = cd_tri->m_tri;
            
            int num_pts = (*it)->num_pts();
            for (int i = 0; i < num_pts; ++i)
            {
                POINT* p = (*it)->Point_of_hse(i);
                POINT* head = findSet(p);
                if (weight(head) != 1) continue;
                
                int npts = 0;
                POINT* ringpts[25];
                PointArrayRing1(p,Hyper_surf_element(tri),Hyper_surf(tri->surf),&npts,ringpts);
                
                std::vector<POINT*> pvec;
                for (int j = 0; j < npts; ++j)
                {
                    if (weight(findSet(ringpts[j])) > 1)
                    {
                        pvec.push_back(ringpts[j]);
                    }
                }

                if (pvec.size() >= 2)
                {
                    pvec.push_back(p);
                    POINT* pts[pvec.size()];
                    std::copy(pvec.begin(),pvec.end(),pts);
                    createImpactZone(pts,pvec.size());
                }
            }
        }
        else
        {
            printf("connectNearbyImpactZones() ERROR: "
                    "unknown CD_HSE_TYPE\n");
            LOC(); clean_up(EXIT_FAILURE);
        }
    }
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
            //skip traversed or isolated pts
		    POINT* pt = (*it)->Point_of_hse(i);
            if (sorted(pt)) continue;

		    POINT* head = findSet(pt);
            if (weight(head) == 1) continue;

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
    
    FT_Save(ft);
    FT_Draw(ft);
}

void CollisionSolver3d::infoImpactZones()
{
	numImpactZones = 0;
	numImpactZonePoints = 0;

	unsortHseList(hseList);
	//unsortHseList(elasticHseList);
    
	//for (auto it = elasticHseList.begin(); it < elasticHseList.end(); ++it)
	for (auto it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            //skip traversed or isolated pts
		    POINT* pt = (*it)->Point_of_hse(i);
            if (sorted(pt)) continue;
		    
            POINT* head = findSet(pt);
            if (weight(head) == 1) continue;
		
            markImpactZonePoints(head);
            numImpactZonePoints += weight(head);
            numImpactZones++;
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
	unsortHseList(movableRigidTriList);

	for (auto it = movableRigidTriList.begin(); it != movableRigidTriList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            //skip traversed or isolated pts
            POINT* pt = (*it)->Point_of_hse(i);
            if (sorted(pt)) continue;
            
		    POINT* head = findSet(pt);
            if (!isMovableRigidBody(pt) || weight(head) == 1)
            {
                sorted(pt) = YES;
                continue;
            }
            
            updateImpactListVelocity(head);
	    }
	}
}

//For jacobi style update of impact zones
void CollisionSolver3d::updateImpactZoneVelocity()
{
	numImpactZones = 0;
	numImpactZonePoints = 0;

	unsortHseList(hseList);
    
	for (auto it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            //skip traversed or isolated pts
		    POINT* pt = (*it)->Point_of_hse(i);
            if (sorted(pt)) continue;

		    POINT* head = findSet(pt);
            if (weight(head) == 1) continue;
            
            updateImpactListVelocity(head);
            numImpactZonePoints += weight(head);
            numImpactZones++;
	    }
	}
}

void updateImpactListVelocity(POINT* head)
{
    //compute impact zone's center of mass position and velocity
    double totalmass = 0.0;
    double avg_dt = 0.0;

    double x_cm[3] = {0.0};
    double v_cm[3] = {0.0};
	
    int num_pts = weight(head);

    POINT* p = head;
	while(p)
    {
		STATE* sl = (STATE*)left_state(p);
        avg_dt += sl->collsn_dt;
        sl->has_collsn = true;

        double m = CollisionSolver3d::getFabricPointMass();
        if (sl->is_stringpt)
            m = CollisionSolver3d::getStringPointMass();
        else if (isStaticRigidBody(p))
        {
            //TODO: Document "impact zone anchoring"
            SURFACE* s = (SURFACE*)Surface_of_hs(p->hs);
            int num_surfpts = I_NumOfSurfInteriorPoints(s);
            m = HUGE/num_pts;
        }
        else if (isMovableRigidBody(p))
        {
            SURFACE* s = (SURFACE*)Surface_of_hs(p->hs);
            int num_surfpts = I_NumOfSurfInteriorPoints(s);
            m = total_mass(p->hs)/num_surfpts;
        }

        totalmass += m;
        for (int i = 0; i < 3; ++i)
        {
		    x_cm[i] += sl->x_old[i]*m; 
		    v_cm[i] += sl->avgVel[i]*m;
		}

        sorted(p) = YES;
        //Still need to mark sorted() for updateImpactZoneVelocityForRG(),
        //and potentially for jacobi style updating
        
        p = next_pt(p);
    }
	
	for (int i = 0; i < 3; ++i)
    {
	    x_cm[i] /= totalmass;
	    v_cm[i] /= totalmass;
	}

    //TODO: Need to use full time step, dt, for movable rigid body
    //      impact zone update.
    //
    double dt = CollisionSolver3d::getTimeStepSize();
    
    //TODO: Is this justified, or just use the full step dt?
    avg_dt += dt; //maybe don't add this...
    avg_dt /= num_pts;
    
        //printf("avg_dt = %g,  dt = %g\n",avg_dt,dt);

	//compute angular momentum
	double L[3] = {0.0};

    p = head;
	while(p)

    {
	    STATE* sl = (STATE*)left_state(p);
        double m = CollisionSolver3d::getFabricPointMass();
        if (sl->is_stringpt)
            m = CollisionSolver3d::getStringPointMass();
        else if (isStaticRigidBody(p))
        {
            SURFACE* s = (SURFACE*)Surface_of_hs(p->hs);
            int num_surfpts = I_NumOfSurfInteriorPoints(s);
            m = HUGE/num_pts;
        }
        else if (isMovableRigidBody(p))
        {
            SURFACE* s = (SURFACE*)Surface_of_hs(p->hs);
            int num_surfpts = I_NumOfSurfInteriorPoints(s);
            m = total_mass(p->hs)/num_pts;
        }
	    
	    double dx[3], dv[3], Li[3];
        minusVec(sl->x_old,x_cm,dx);
	    minusVec(sl->avgVel,v_cm,dv); 	
	    Cross3d(dx,dv,Li);
	    scalarMult(m,Li,Li);
	    addVec(Li,L,L);    
	    p = next_pt(p);
	}

	//compute Inertia tensor
	double I[3][3] = {0.0};

	p = head;
	while(p)
    {
	    STATE* sl = (STATE*)left_state(p);
        double m = CollisionSolver3d::getFabricPointMass();
        if (sl->is_stringpt)
            m = CollisionSolver3d::getStringPointMass();
        else if (isStaticRigidBody(p))
        {
            SURFACE* s = (SURFACE*)Surface_of_hs(p->hs);
            int num_surfpts = I_NumOfSurfInteriorPoints(s);
            m = HUGE/num_pts;
        }
        else if (isMovableRigidBody(p))
        {
            SURFACE* s = (SURFACE*)Surface_of_hs(p->hs);
            int num_surfpts = I_NumOfSurfInteriorPoints(s);
            m = total_mass(p->hs)/num_pts;
        }

	    double dx[3];
        minusVec(sl->x_old,x_cm,dx);
	    double mag_dx = Mag3d(dx);

        for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < 3; ++j)
        {
		    double val = -dx[i]*dx[j];
            if (i == j)
                val += mag_dx*mag_dx; 
    	 	I[i][j] += val*m;
	    }

        p = next_pt(p);
	}

	//compute angular velocity w: I*w = L;
	double w[3];
    
    if (myDet3d(I) < ROUND_EPS)
    {
        //I is non-invertible, calculate pseudoinverse with SVD
        arma::mat arI(3, 3);
        arma::vec arL(3);

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
             arI(i, j) = I[i][j];
        for (int i = 0; i < 3; i++)
             arL(i) = L[i];

        arma::mat arU;
        arma::mat arV;
        arma::vec ars;

        arma::svd(arU, ars, arV, arI);
        for (int i = 0; i < 3; i++)
        {
            if (ars(i)) ars(i) = 1.0/ars(i);
        }
        arma::mat pinvarI = arV*arma::diagmat(ars)*arU.t();
        arma::vec arw = pinvarI*arL;

        for (int i = 0; i < 3; i++)
            w[i] = arw[i];
    }
    else
    {
        double tmp[3][3];
        for (int i = 0; i < 3; ++i)
        {
            memcpy(tmp,I,9*sizeof(double));
            for (int j = 0; j < 3; j++)
                tmp[j][i] = L[j];
            w[i] = myDet3d(tmp)/myDet3d(I);
        }
    }
    double mag_w = Mag3d(w);
	
	//compute average velocity for each point
	p = head;
    while(p)
    {
        if (isStaticRigidBody(p))
        {
            p = next_pt(p);
            continue;
	    }

	    double x_new[3],dx[3];
	    double xF[3], xR[3];
	    double wxR[3],tmpV[3];
	    
        STATE* sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    
        if (mag_w < ROUND_EPS)
        {
	        for (int i = 0; i < 3; ++i)
            {
                xF[i] = dx[i];
                wxR[i] = 0.0;
            }
		    minusVec(dx,xF,xR);
	    }
	    else
        {
	        scalarMult(Dot3d(dx,w)/Dot3d(w,w),w,xF);
	        minusVec(dx,xF,xR);
            scalarMult(sin(avg_dt*mag_w)/mag_w,w,tmpV);
                //scalarMult(sin(dt*mag_w)/mag_w,w,tmpV);
	        Cross3d(tmpV,xR,wxR);
	    }

	    for (int i = 0; i < 3; ++i)
	    {
            /*
            x_new[i] = x_cm[i] + dt*v_cm[i] + xF[i]
                       + cos(dt*mag_w)*xR[i] + wxR[i];
	
		    sl->avgVel[i] = (x_new[i] - sl->x_old[i])/dt;
            */
	    	
            x_new[i] = x_cm[i] + avg_dt*v_cm[i] + xF[i]
                       + cos(avg_dt*mag_w)*xR[i] + wxR[i];

		    sl->avgVel[i] = (x_new[i] - sl->x_old[i])/dt;
		        //sl->avgVel[i] = (x_new[i] - sl->x_old[i])/avg_dt;

            //TODO: see above todo regarding use of full time step
            //      for movable rigid body impact zones
            
            if (std::isnan(sl->avgVel[i]))
            { 
                printf("x_old[3], avgVel[3]\n");
                
                p = head;
                while(p)
                {
                    STATE* sl = (STATE*)left_state(p);
                    printf("%f %f %f %f %f %f;\n",
                            sl->x_old[0],sl->x_old[1],sl->x_old[2],
                            sl->avgVel[0],sl->avgVel[1],sl->avgVel[2]);
                    p = next_pt(p);
                }

                printf("num_pts = %d, weight = %d\n",num_pts,weight(head));
                printf("nan vel, w = %f, mag_w = %f\n",w[i],mag_w);
                printf("L = [%f %f %f]\n",L[0],L[1],L[2]);
                printf("I = [%f %f %f;  %f %f %f; %f %f %f]\n",
                        I[0][0],I[0][1],I[0][2],
                        I[1][0],I[1][1],I[1][2],
                        I[2][0],I[2][1],I[2][2]);
                printf("xF = %f %f %f, xR = %f %f %f\n",
                        xF[0],xF[1],xF[2],xR[0],xR[1],xR[2]);
            
                clean_up(ERROR);
            }
	    }

        p = next_pt(p);
	}
}

void CollisionSolver3d::turnOnGsUpdate() {gs_update = true;}
void CollisionSolver3d::turnOffGsUpdate() {gs_update = false;}
bool CollisionSolver3d::getGsUpdateStatus() {return gs_update;}

//jacobi iteration
void CollisionSolver3d::limitStrainPosnJac()
{
    double dt = getTimeStepSize();
    if (dt < 1.0e-04) return;

	const int MAX_ITER = 2;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        int numBondStrain = computeStrainImpulsesPosn(stringBondList);
        int numTriStrain = computeStrainImpulsesPosn(fabricTriList);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Edges\n",numBondStrain);
            printf("%d TRI Strain Edges\n",numTriStrain);
        }

        if (numBondStrain == 0 && numTriStrain == 0) break;

        applyStrainImpulses(MotionState::STATIC);
	}
}

void CollisionSolver3d::limitStrainPosnGS()
{
    double dt = getTimeStepSize();
    if (dt < 1.0e-04) return;

    turnOnGsUpdate();

	const int MAX_ITER = 2;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        int numBondStrain = computeStrainImpulsesPosn(stringBondList);
        int numTriStrain = computeStrainImpulsesPosn(fabricTriList);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Edges\n",numBondStrain);
            printf("%d TRI Strain Edges\n",numTriStrain);
        }

        if (numBondStrain == 0 && numTriStrain == 0) break;
	}
    
    turnOffGsUpdate();
}

int CollisionSolver3d::computeStrainImpulsesPosn(std::vector<CD_HSE*>& list)
{
    double dt = getTimeStepSize();
    double TOL = getStrainLimit();
    double CTOL = getCompressiveStrainLimit();
    bool gauss_seidel = getGsUpdateStatus();

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
            
            if (delta_len0 > TOL*len0 || delta_len0 < -1.0*CTOL*len0)
            {
                double I;
                if (delta_len0 > TOL*len0) //Tension
                { 
                    I = 0.5*(delta_len0 - TOL*len0)/dt;
                }
                else                       //Compression
                {
                    I = 0.5*(delta_len0 + CTOL*len0)/dt;
                }

                double vec01[MAXD];
                Pts2Vec(p[0],p[1],vec01);
                scalarMult(1.0/lnew,vec01,vec01);
                
                //Do not apply impulses to fixed nodes
                if (!isRigidBody(sl[0]) && !isRigidBody(sl[1]))
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        sl[0]->strainImpulse[j] += I*vec01[j];
                        sl[1]->strainImpulse[j] -= I*vec01[j];
                    }
                    sl[0]->strain_num++;
                    sl[1]->strain_num++;
                } 
                else if (!isRigidBody(sl[0]) && isRigidBody(sl[1]))
                {
                    for (int j = 0; j < 3; ++j)
                        sl[0]->strainImpulse[j] += 2.0*I*vec01[j];
                    sl[0]->strain_num++;
                }
                else if (isRigidBody(sl[0]) && !isRigidBody(sl[1]))
                {
                    for (int j = 0; j < 3; ++j)
                        sl[1]->strainImpulse[j] -= 2.0*I*vec01[j];
                    sl[1]->strain_num++;
                }
                else
                {
                    //NOTE: 2 "registered points" could show up here
                    continue;
                    /*
                    printf("ERROR: \n");
                    LOC(); clean_up(EXIT_FAILURE);
                    */
                }

                numStrainEdges++;
                
                if (gauss_seidel)
                {
                    for (int k = 0; k < 2; ++k)
                    {
                        if (sl[k]->strain_num > 0)
                        {
                            sl[k]->has_strainlim_prox = true;
                            for (int j = 0; j < 3; ++j)
                            {
                                sl[k]->avgVel[j] += sl[k]->strainImpulse[j];
                                sl[k]->strainImpulse[j] = 0.0;
                            }
                            sl[k]->strain_num = 0;
                        }
                    }
                }

            }
        }
    }
    
    return numStrainEdges;
}

//gauss-seidel iteration
void CollisionSolver3d::limitStrainRatePosnGS(MotionState mstate)
{
    turnOnGsUpdate();
	
    const int MAX_ITER = 2;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        int numBondStrainRate = computeStrainRateImpulsesPosn(stringBondList,mstate);
        int numTriStrainRate = computeStrainRateImpulsesPosn(fabricTriList,mstate);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Rate Edges\n",numBondStrainRate);
            printf("%d TRI Strain Rate Edges\n",numTriStrainRate);
        }

        if (numBondStrainRate == 0 && numTriStrainRate == 0) break;
	}

    turnOffGsUpdate();
}

//jacobi iteration
void CollisionSolver3d::limitStrainRatePosnJac(MotionState mstate)
{
	const int MAX_ITER = 2;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        int numBondStrainRate = computeStrainRateImpulsesPosn(stringBondList,mstate);
        int numTriStrainRate = computeStrainRateImpulsesPosn(fabricTriList,mstate);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Rate Edges\n",numBondStrainRate);
            printf("%d TRI Strain Rate Edges\n",numTriStrainRate);
        }

        if (numBondStrainRate == 0 && numTriStrainRate == 0) break;

        applyStrainImpulses(mstate);
	}
}

int CollisionSolver3d::computeStrainRateImpulsesPosn(
        std::vector<CD_HSE*>& list,
        MotionState mstate)
{
    double dt = getTimeStepSize();
    double TOL = getStrainRateLimit();
    bool gauss_seidel = getGsUpdateStatus();

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
                    if (Gindex(tri_nb) < Gindex(tri)) continue;
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

                double vec01[MAXD];
                Pts2Vec(p[0],p[1],vec01);
                scalarMult(1.0/lnew,vec01,vec01);
                
                //Do not apply impulses to fixed nodes
                if (!isRigidBody(sl[0]) && !isRigidBody(sl[1]))
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        sl[0]->strainImpulse[j] += I*vec01[j];
                        sl[1]->strainImpulse[j] -= I*vec01[j];
                    }
                    sl[0]->strain_num++;
                    sl[1]->strain_num++;
                } 
                else if (!isRigidBody(sl[0]) && isRigidBody(sl[1]))
                {
                    for (int j = 0; j < 3; ++j)
                        sl[0]->strainImpulse[j] += 2.0*I*vec01[j];
                    sl[0]->strain_num++;
                }
                else if (isRigidBody(sl[0]) && !isRigidBody(sl[1]))
                {
                    for (int j = 0; j < 3; ++j)
                        sl[1]->strainImpulse[j] -= 2.0*I*vec01[j];
                    sl[1]->strain_num++;
                }
                else
                {
                    //NOTE: two "registered points" could show up here
                    continue;
                    /*
                    printf("ERROR: \n");
                    LOC(); clean_up(EXIT_FAILURE);
                    */
                }
                
                numStrainRateEdges++;

                if (gauss_seidel)
                {
                    for (int k = 0; k < 2; ++k)
                    {
                        if (sl[k]->strain_num > 0)
                        {
                            if (mstate == MotionState::STATIC)
                                sl[k]->has_strainlim_prox = true;
                            else
                                sl[k]->has_strainlim_collsn = true;

                            for (int j = 0; j < 3; ++j)
                            {
                                sl[k]->avgVel[j] += sl[k]->strainImpulse[j];
                                sl[k]->strainImpulse[j] = 0.0;
                            }
                            sl[k]->strain_num = 0;
                        }
                    }
                }

            }
        }
    }
    
    return numStrainRateEdges;
}

//gauss-seidel iteration
void CollisionSolver3d::limitStrainVelGS()
{
    turnOnGsUpdate();

    const int MAX_ITER = 2;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        int numBondStrainVel = computeStrainImpulsesVel(stringBondList);
        int numTriStrainVel = computeStrainImpulsesVel(fabricTriList);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Velocity Edges\n",numBondStrainVel);
            printf("%d TRI Strain Velocity Edges\n",numTriStrainVel);
        }

        if (numBondStrainVel == 0 && numTriStrainVel == 0) break;
	}
    
    turnOffGsUpdate();
}

//jacobi iteration
void CollisionSolver3d::limitStrainVelJAC()
{
    const int MAX_ITER = 2;
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        int numBondStrainVel = computeStrainImpulsesVel(stringBondList);
        int numTriStrainVel = computeStrainImpulsesVel(fabricTriList);
        
        if (debugging("strain_limiting"))
        {
            printf("%d BOND Strain Velocity Edges\n",numBondStrainVel);
            printf("%d TRI Strain Velocity Edges\n",numTriStrainVel);
        }

        if (numBondStrainVel == 0 && numTriStrainVel == 0) break;

        applyStrainImpulses(MotionState::POSTCOLLISION);
	}
}

int CollisionSolver3d::computeStrainImpulsesVel(std::vector<CD_HSE*>& list)
{
    bool gauss_seidel = getGsUpdateStatus();

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
                    if (Gindex(tri_nb) < Gindex(tri)) continue;
                }
            }

            p[0] = (*it)->Point_of_hse(i%np);
            p[1] = (*it)->Point_of_hse((i+1)%np);

            sl[0] = (STATE*)left_state(p[0]);
            sl[1] = (STATE*)left_state(p[1]);
            
            //skip edges that did not get strain limiting impulses
            //in limitStrainPos() or limitStrainRatePos()
            if (!(sl[0]->has_strainlim_prox || sl[0]->has_strainlim_collsn) &&
                !(sl[1]->has_strainlim_prox || sl[1]->has_strainlim_collsn)) continue;

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
            if (fabs(vcomp01) < MACH_EPS) continue;
                
            double I = 0.5*vcomp01;

            //Do not apply impulses to rg_string_nodes
            if (!isRigidBody(sl[0]) && !isRigidBody(sl[1]))
            {
                for (int j = 0; j < 3; ++j)
                {
                    sl[0]->strainImpulse[j] += I*vec01[j];
                    sl[1]->strainImpulse[j] -= I*vec01[j];
                }
                sl[0]->strain_num++;
                sl[1]->strain_num++;
            } 
            else if (!isRigidBody(sl[0]) && isRigidBody(sl[1]))
            {
                for (int j = 0; j < 3; ++j)
                    sl[0]->strainImpulse[j] += 2.0*I*vec01[j];
                sl[0]->strain_num++;
            }
            else if (isRigidBody(sl[0]) && !isRigidBody(sl[1]))
            {
                for (int j = 0; j < 3; ++j)
                    sl[1]->strainImpulse[j] -= 2.0*I*vec01[j];
                sl[1]->strain_num++;
            }
            else
            {
                //NOTE: two "registered points" could show up here
                continue;
                /*
                printf("ERROR: \n");
                LOC(); clean_up(EXIT_FAILURE);
                */
            }
                
            numRelVelStrainEdges++;

            if (gauss_seidel)
            {
                for (int k = 0; k < 2; ++k)
                {
                    if (sl[k]->strain_num > 0)
                    {
                        sl[k]->has_strainlim_collsn = true;

                        /*
                        if (mstate == MotionState::STATIC)
                            sl[k]->has_strainlim_prox = true;
                        else
                            sl[k]->has_strainlim_collsn = true;
                        */

                        for (int j = 0; j < 3; ++j)
                        {
                            sl[k]->avgVel[j] += sl[k]->strainImpulse[j];
                            sl[k]->strainImpulse[j] = 0.0;
                        }
                        sl[k]->strain_num = 0;
                    }
                }
            }
        
        }
    }
    
    return numRelVelStrainEdges;
}

void CollisionSolver3d::applyStrainImpulses(MotionState mstate)
{
	POINT *p;
	STATE *sl;
	double maxSpeed = 0;
	double* maxVel = nullptr;

	unsortHseList(elasticHseList);
    
	for (auto it = elasticHseList.begin(); it < elasticHseList.end(); ++it)
    {
        int np = (*it)->num_pts(); 
	    for (int j = 0; j < np; ++j)
	    {
            p = (*it)->Point_of_hse(j);
            if (sorted(p)) continue;

            sl = (STATE*)left_state(p);
            if (sl->strain_num > 0)
            {
                if (mstate == MotionState::STATIC)
                    sl->has_strainlim_prox = true;
                else
                    sl->has_strainlim_collsn = true;

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

//functions for UF alogrithm
int& weight(POINT* p)
{
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

    for (auto it = hseList.begin(); it < hseList.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sorted(pt) = NO;
            sl = (STATE*)left_state(pt);
            sl->impZone.next_pt = nullptr;
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

bool isStaticRigidBody(const STATE* sl)
{
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

bool isMovableRigidBody(const STATE* sl)
{
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

bool isRigidBody(const STATE* sl)
{
    return sl->is_fixed || sl->is_movableRG;
}

bool isRigidBody(const CD_HSE* hse)
{
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
