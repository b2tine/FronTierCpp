//#include <armadillo>
#include "collid.h"
#include <iostream>

static bool MovingEdgeToEdgeGS(POINT**);
static bool MovingEdgeToEdgeJac(POINT**);

static bool EdgeToEdge(POINT**, double,
        MotionState mstate = MotionState::STATIC, double root = -1.0);

static void EdgeToEdgeImpulse(POINT**,double*,double,double,double,double,MotionState,double);
static void EdgeToEdgeInelasticImpulse(double,POINT**,double*,double*,double*);
static void EdgeToEdgeElasticImpulse(double,double,double,POINT**,double*,double*,double,double,double);

static bool MovingPointToTriGS(POINT**);
static bool MovingPointToTriJac(POINT**);

static bool PointToTri(POINT**, double,
        MotionState mstate = MotionState::STATIC, double root = -1.0);

static void PointToTriImpulse(POINT**,double*,double*,double,double,MotionState,double);
static void PointToTriInelasticImpulse(double,POINT**,double*,double*,double*,double*);
static void PointToTriElasticImpulse(double,double,double,POINT**,double*,double*,double,double,double);

static bool isCoplanar(POINT**,double,double*);

static double getPointMass(POINT* pt);
static double getPointFrictionConstant(POINT* pt);



bool MovingTriToBondGS(const TRI* tri,const BOND* bd)
{
    bool status = false;

	POINT* pts[4];
	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	/* detect collision of start point of bond w.r.t to tri */
	pts[3] = bd->start;
    if (MovingPointToTriGS(pts))
        status = true;
    
    /* detect collision of end point of bond to w.r.t. tri */
	pts[3] = bd->end;
    if (MovingPointToTriGS(pts))
        status = true;

    /* detect collision of each of tri edge w.r.t to bond */
	pts[2] = bd->start;
	pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
        if (MovingEdgeToEdgeGS(pts))
            status = true;
	}

    return status;
}

bool MovingBondToBondGS(const BOND* b1, const BOND* b2)
{
	POINT* pts[4];

	pts[0] = b1->start;
	pts[1] = b1->end;
	pts[2] = b2->start;
	pts[3] = b2->end;

	bool status = false;
    if(MovingEdgeToEdgeGS(pts))
        status = true;
    
    return status;
}

bool MovingTriToTriGS(const TRI* a, const TRI* b)
{
	POINT* pts[4];
	bool status = false;

	//detect point to tri collision
	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
    {
	    const TRI* tmp_tri1 = (k == 0) ? a : b;
	    const TRI* tmp_tri2 = (k == 0) ? b : a;
	    for (int j = 0; j < 3; ++j)
            pts[j] = Point_of_tri(tmp_tri1)[j];
        pts[3] = Point_of_tri(tmp_tri2)[i];

        if (MovingPointToTriGS(pts))
            status = true;
	}

	//detect edge to edge collision
	for (int i = 0; i < 3; ++i)
    {
        pts[0] = Point_of_tri(a)[i];
        pts[1] = Point_of_tri(a)[(i+1)%3];
        for (int j = 0; j < 3; ++j)
        {
            pts[2] = Point_of_tri(b)[j];
            pts[3] = Point_of_tri(b)[(j+1)%3];
		
            if (MovingEdgeToEdgeGS(pts))
                status = true;
	    }
    }

    /*
    if (status)
    {
        collisionPairsList.push_back();
    }
    */

    return status;
}

bool MovingTriToBondJac(const TRI* tri,const BOND* bd)
{
    bool status = false;

	POINT* pts[4];
	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	/* detect collision of start point of bond w.r.t to tri */
	pts[3] = bd->start;
    if (MovingPointToTriJac(pts))
        status = true;
    
    /* detect collision of end point of bond to w.r.t. tri */
	pts[3] = bd->end;
    if (MovingPointToTriJac(pts))
        status = true;

    /* detect collision of each of tri edge w.r.t to bond */
	pts[2] = bd->start;
	pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
        if (MovingEdgeToEdgeJac(pts))
            status = true;
	}

    return status;
}

bool MovingBondToBondJac(const BOND* b1, const BOND* b2)
{
	POINT* pts[4];

	pts[0] = b1->start;
	pts[1] = b1->end;
	pts[2] = b2->start;
	pts[3] = b2->end;

	bool status = false;
    if(MovingEdgeToEdgeJac(pts))
        status = true;
    
    return status;
}

bool MovingTriToTriJac(const TRI* a, const TRI* b)
{
	POINT* pts[4];
	bool status = false;

	//detect point to tri collision
	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
    {
	    const TRI* tmp_tri1 = (k == 0) ? a : b;
	    const TRI* tmp_tri2 = (k == 0) ? b : a;
	    for (int j = 0; j < 3; ++j)
            pts[j] = Point_of_tri(tmp_tri1)[j];
        pts[3] = Point_of_tri(tmp_tri2)[i];

        if (MovingPointToTriJac(pts))
            status = true;
	}

	//detect edge to edge collision
	for (int i = 0; i < 3; ++i)
    {
        pts[0] = Point_of_tri(a)[i];
        pts[1] = Point_of_tri(a)[(i+1)%3];
        for (int j = 0; j < 3; ++j)
        {
            pts[2] = Point_of_tri(b)[j];
            pts[3] = Point_of_tri(b)[(j+1)%3];
		
            if (MovingEdgeToEdgeJac(pts))
                status = true;
	    }
    }

    /*
    if (status)
    {
        collisionPairsList.push_back();
    }
    */

    return status;
}

//For use with jacobi avgVel update.
//Saves impulses to be applied to avgVel.
static bool MovingPointToTriJac(POINT* pts[])
{
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};

    double tol = CollisionSolver3d::getFabricRoundingTolerance();
    /*
    STATE* s = (STATE*)left_state(pts[3]);
    if (s->is_stringpt)
        tol = CollisionSolver3d::getStringRoundingTolerance();
    */  
    
    bool status = false;
	if (isCoplanar(pts,dt,roots))
    {
        for (int i = 0; i < 4; ++i)
        {
            if (roots[i] < 0)
                continue;
    
            for (int j = 0; j < 4; ++j)
            {
                STATE* sl = (STATE*)left_state(pts[j]);
    		    for (int k = 0; k < 3; ++k)
                    Coords(pts[j])[k] = sl->x_old[k] + roots[i]*sl->avgVel[k];
		    }
    
            MotionState mstate = MotionState::MOVING;
            if (PointToTri(pts,tol,mstate,roots[i]))
            {
                status = true;
                for (int j = 0; j < 4; ++j)
                {
                    STATE* sl = (STATE*)left_state(pts[j]);
                    sl->collsn_dt = roots[i];
                }
                
                CollisionSolver3d::addCollisionTime(roots[i]);
                break;
            }
	    }
	}

    for (int j = 0; j < 4; ++j)
    {
        STATE* sl = (STATE*)left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k];
    }
    
    bool rigid_body_point = false;
    if (isRigidBody(pts[0]) || isRigidBody(pts[3]))
        rigid_body_point = true;
    
    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    //if (status && is_detImpZone && !rigid_body_point)
    if (status && is_detImpZone)
    {
        createImpZone(pts,4);//only collects fabric points
    }
    
    return status;
}

//gauss-seidel update
static bool MovingPointToTriGS(POINT* pts[])
{
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};
    double rtol = CollisionSolver3d::getFabricRoundingTolerance();
    double ftol = CollisionSolver3d::getFabricThickness();
    //
    //STATE* s = (STATE*)left_state(pts[3]);
    //if (s->is_stringpt)
    //    rtol = CollisionSolver3d::getStringRoundingTolerance();
    //
        
    bool status = false;
    double collision_dt;

	if (isCoplanar(pts,dt,roots))
    {
        for (int i = 0; i < 4; ++i)
        {
            if (roots[i] < 0) continue;
    
            for (int j = 0; j < 4; ++j)
            {
                STATE* sl = (STATE*)left_state(pts[j]);
    		    for (int k = 0; k < 3; ++k)
                    Coords(pts[j])[k] = sl->x_old[k] + roots[i]*sl->avgVel[k];
		    }

            MotionState mstate = MotionState::MOVING;
            if (PointToTri(pts,rtol,mstate,roots[i]))
            {
                status = true;
                for (int j = 0; j < 4; ++j)
                {
                    STATE* sl = (STATE*)left_state(pts[j]);
                    sl->collsn_dt = roots[i];
                }
                    
                collision_dt = roots[i];
                CollisionSolver3d::addCollisionTime(roots[i]);
                break;
            }
	    }
	}

    //return to original position
    for (int j = 0; j < 4; ++j)
    {
        STATE* sl = (STATE*)left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k];
    }

    bool rigid_body_point = false;
    if (isRigidBody(pts[0]) || isRigidBody(pts[3]))
        rigid_body_point = true;
    
    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    if (status && (!is_detImpZone || rigid_body_point))
    {
        /*
        if (is_detImpZone && rigid_body_point)
        {
            printf("Point To Tri Orphans\n");
            printPointSetCollisionStats(pts,4);
        }
        */
        
        for (int j = 0; j < 4; ++j)
        {
            STATE* sl = (STATE*)left_state(pts[j]);
            if (is_detImpZone && rigid_body_point)
            {
                sl->impzone_orphan = true;
            }

            if (sl->collsn_num > 0)
            {
                sl->has_collsn = true;
                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] += sl->collsnImpulse[k];
                    sl->collsnImpulse[k] = 0.0;
                }
                sl->collsn_num = 0;
            }

            //move to new candidate position
            for (int k = 0; k < 3; ++k)
                Coords(pts[j])[k] = sl->x_old[k] + sl->avgVel[k]*sl->collsn_dt;
        }
            
        //Apply elastic repulsion impulse if still not well separated
        MotionState pmstate = MotionState::POSTCOLLISION;
        if (PointToTri(pts,ftol,pmstate,collision_dt))
        {
            /*
            if (is_detImpZone && rigid_body_point)
            {
                printf("Point To Tri Orphans POSTCOLLISION\n");
                printPointSetCollisionStats(pts,4);
            }
            */
        
            //apply elastic repulsion impulse if still not well separated
            for (int j = 0; j < 4; ++j)
            {
                STATE* sl = (STATE*)left_state(pts[j]);
                if (sl->collsn_num > 0)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        sl->avgVel[k] += sl->collsnImpulse[k];
                        sl->collsnImpulse[k] = 0.0;
                    }
                    sl->collsn_num = 0;
                }

                //compute new candidate position and effective velocity
                double x_new[3];
                for (int k = 0; k < 3; ++k)
                {
                    x_new[k] = sl->x_old[k] + sl->avgVel[k]*sl->collsn_dt;
                    sl->avgVel[k] = (x_new[k] - sl->x_old[k])/dt;
                    //NOTE: Using collsn_dt in avgVel computation can
                    //      cause penetration and loss of side orientation.
                }
            }
        }

        //return to original position
        for (int j = 0; j < 4; ++j)
        {
            STATE* sl = (STATE*)left_state(pts[j]);
            for (int k = 0; k < 3; ++k)
                Coords(pts[j])[k] = sl->x_old[k];
        }

    }
    else if (status && is_detImpZone)
    {
        //Note: no rigid body points
        createImpZone(pts,4);
        POINT* head = findSet(pts[0]);
        updateImpactListVelocity(head);
    }
    
    return status;
}

//For use with jacobi avgVel update.
//Saves impulses to be applied to avgVel.
static bool MovingEdgeToEdgeJac(POINT* pts[])
{
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};

    STATE* s0 = (STATE*)left_state(pts[0]);
    STATE* s2 = (STATE*)left_state(pts[2]);

    double tol = CollisionSolver3d::getFabricRoundingTolerance();
    if (s0->is_stringpt && s2->is_stringpt)
        tol = CollisionSolver3d::getStringRoundingTolerance();

    bool status = false;
	if (isCoplanar(pts,dt,roots))
    {
        for (int i = 0; i < 4; ++i)
        {
            if (roots[i] < 0) continue;
                
            for (int j = 0; j < 4; ++j)
            {
                STATE* sl = (STATE*)left_state(pts[j]);
                for (int k = 0; k < 3; ++k)
                    Coords(pts[j])[k] = sl->x_old[k] + roots[i]*sl->avgVel[k];
            }

            MotionState mstate = MotionState::MOVING;
            if (EdgeToEdge(pts,tol,mstate,roots[i]))
            {
                status = true;
                for (int j = 0; j < 4; ++j)
                {
                    STATE* sl = (STATE*)left_state(pts[j]);
                    sl->collsn_dt = roots[i];
                }

                CollisionSolver3d::addCollisionTime(roots[i]);
                break;
            }
        }
    }

    //return to original position
    for (int j = 0; j < 4; ++j)
    {
        STATE* sl = (STATE*)left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k];
    }
    
    bool rigid_body_point = false;
    if (isRigidBody(pts[0]) || isRigidBody(pts[3]))
        rigid_body_point = true;
    
    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    //if (status && is_detImpZone && !rigid_body_point)
    if (status && is_detImpZone)
    {
        createImpZone(pts,4);//only collects fabric points
        //createImpactZone(pts,4);
    }
    
    return status;
}

//gauss-seidel update
static bool MovingEdgeToEdgeGS(POINT* pts[])
{
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};

    STATE* s0 = (STATE*)left_state(pts[0]);
    STATE* s2 = (STATE*)left_state(pts[2]);

    double rtol = CollisionSolver3d::getFabricRoundingTolerance();
    double ftol = CollisionSolver3d::getFabricThickness();
    if (s0->is_stringpt && s2->is_stringpt)
    {
        rtol = CollisionSolver3d::getStringRoundingTolerance();
        ftol = CollisionSolver3d::getStringThickness();
    }

    bool status = false;
    double collision_dt;

	if (isCoplanar(pts,dt,roots))
    {
        for (int i = 0; i < 4; ++i)
        {
            if (roots[i] < 0) continue;
                
            for (int j = 0; j < 4; ++j)
            {
                STATE* sl = (STATE*)left_state(pts[j]);
                for (int k = 0; k < 3; ++k)
                    Coords(pts[j])[k] = sl->x_old[k] + roots[i]*sl->avgVel[k];
            }

            MotionState mstate = MotionState::MOVING;
            if (EdgeToEdge(pts,rtol,mstate,roots[i]))
            {
                status = true;
                for (int j = 0; j < 4; ++j)
                {
                    STATE* sl = (STATE*)left_state(pts[j]);
                    sl->collsn_dt = roots[i];
                }
                
                collision_dt = roots[i];
                CollisionSolver3d::addCollisionTime(roots[i]);
                break;
            }
        }
    }

    //return to original position
    for (int j = 0; j < 4; ++j)
    {
        STATE* sl = (STATE*)left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k];
    }

    //TODO: ALLOW IMPACT ZONES FOR STRING-STRING INTERACTIONS FOR NOW.
    bool string_string = false;
    //
    //No Impact Zones for string-string interactions
    //bool string_string = false;
    //if (s0->is_stringpt && s2->is_stringpt)
    //    string_string = true;
    //

    bool rigid_body_point = false;
    if (isRigidBody(pts[0]) || isRigidBody(pts[3]))
        rigid_body_point = true;
    
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    //if (status && (!is_detImpZone || string_string))
    if (status && (!is_detImpZone || rigid_body_point || string_string))
    {
        /*
        if (is_detImpZone && rigid_body_point)
        {
            printf("Edge To Edge Orphans\n");
            printPointSetCollisionStats(pts,4);
        }
        */
        
        for (int j = 0; j < 4; ++j)
        {
            STATE* sl = (STATE*)left_state(pts[j]);
            
            /*
            if (is_detImpZone && rigid_body_point)
            {
                sl->impzone_orphan = true;
            }
            */

            if (sl->collsn_num > 0)
            {
                sl->has_collsn = true;
                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] += sl->collsnImpulse[k];
                    sl->collsnImpulse[k] = 0.0;
                }
                sl->collsn_num = 0;
            }

            //move to new candidate position
            for (int k = 0; k < 3; ++k)
                Coords(pts[j])[k] = sl->x_old[k] + sl->avgVel[k]*sl->collsn_dt;
        }
        
        //Apply elastic repulsion impulse if still not well separated
        MotionState pmstate = MotionState::POSTCOLLISION;
        if (EdgeToEdge(pts,ftol,pmstate,collision_dt))
        {
            /*
            if (is_detImpZone && rigid_body_point)
            {
                printf("Edge To Edge Orphans POSTCOLLISION\n");
                printPointSetCollisionStats(pts,4);
            }
            */
        
            for (int j = 0; j < 4; ++j)
            {
                STATE* sl = (STATE*)left_state(pts[j]);
                if (sl->collsn_num > 0)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        sl->avgVel[k] += sl->collsnImpulse[k];
                        sl->collsnImpulse[k] = 0.0;
                    }
                    sl->collsn_num = 0;
                }
                
                //compute new candidate position and effective velocity
                double x_new[3];
                for (int k = 0; k < 3; ++k)
                {
                    x_new[k] = sl->x_old[k] + sl->avgVel[k]*sl->collsn_dt;
                    sl->avgVel[k] = (x_new[k] - sl->x_old[k])/dt;
                    //NOTE: Using collsn_dt in avgVel computation can
                    //      cause penetration and loss of side orientation.
                }
            }
        }
        
        //return to original position
        for (int j = 0; j < 4; ++j)
        {
            STATE* sl = (STATE*)left_state(pts[j]);
            for (int k = 0; k < 3; ++k)
                Coords(pts[j])[k] = sl->x_old[k];
        }

    }
    else if (status && is_detImpZone && !string_string)
    {
        //Note: no rigid body points
        createImpZone(pts,4);
        //createImpactZone(pts,4);
        POINT* head = findSet(pts[0]);
        updateImpactListVelocity(head);
    }

    return status;
}

static void isCoplanarHelper(double* s[], double v[][3])
{
    v[0][0] = s[0][0];         v[0][1] = s[0][1];         v[0][2] = s[0][2];
	v[1][0] = s[1][0]-s[0][0]; v[1][1] = s[1][1]-s[0][1]; v[1][2] = s[1][2]-s[0][2];
	v[2][0] = s[2][0]-s[0][0]; v[2][1] = s[2][1]-s[0][1]; v[2][2] = s[2][2]-s[0][2];
	v[3][0] = s[3][0]-s[0][0]; v[3][1] = s[3][1]-s[0][1]; v[3][2] = s[3][2]-s[0][2];
}

static bool isCoplanar(POINT* pts[], const double dt, double roots[])
{
    if (debugging("collision"))
        CollisionSolver3d::is_coplanar++;

	double v[4][3] = {0.0};
    double x[4][3] = {0.0};
	
    double* tmp[4] = {nullptr};

	//for performance, unrolling the loop
	tmp[0] = ((STATE*)left_state(pts[0]))->avgVel;
	tmp[1] = ((STATE*)left_state(pts[1]))->avgVel;
	tmp[2] = ((STATE*)left_state(pts[2]))->avgVel;
	tmp[3] = ((STATE*)left_state(pts[3]))->avgVel;
	isCoplanarHelper(tmp, v);

	tmp[0] = ((STATE*)left_state(pts[0]))->x_old;
	tmp[1] = ((STATE*)left_state(pts[1]))->x_old;
	tmp[2] = ((STATE*)left_state(pts[2]))->x_old;
	tmp[3] = ((STATE*)left_state(pts[3]))->x_old;
	isCoplanarHelper(tmp, x);

	//get roots "t" of a cubic equation
	//(x12 + t*v12) x (x13 + t*v13)*(x14 + t*v14) = 0
	//transform to a*t^3 + b*t^2 + c*t + d = 0
	double a, b, c, d;
	double vv[3], vx[3], xx[3];
	vv[0] = v[1][1]*v[2][2]-v[1][2]*v[2][1];
	vv[1] = v[1][0]*v[2][2]-v[1][2]*v[2][0];
	vv[2] = v[1][0]*v[2][1]-v[1][1]*v[2][0];
	
	vx[0] = v[1][1]*x[2][2]-v[1][2]*x[2][1]-v[2][1]*x[1][2]+v[2][2]*x[1][1];
	vx[1] = v[1][0]*x[2][2]-v[1][2]*x[2][0]-v[2][0]*x[1][2]+v[2][2]*x[1][0];
	vx[2] = v[1][0]*x[2][1]-v[1][1]*x[2][0]-v[2][0]*x[1][1]+v[2][1]*x[1][0];

	xx[0] = x[1][1]*x[2][2]-x[1][2]*x[2][1];
	xx[1] = x[1][0]*x[2][2]-x[1][2]*x[2][0];
	xx[2] = x[1][0]*x[2][1]-x[1][1]*x[2][0];

	a = v[3][0]*vv[0] - v[3][1]*vv[1] + v[3][2]*vv[2];

	b = x[3][0]*vv[0] - x[3][1]*vv[1] + x[3][2]*vv[2] + 
	    v[3][0]*vx[0] - v[3][1]*vx[1] + v[3][2]*vx[2];

	c = x[3][0]*vx[0] - x[3][1]*vx[1] + x[3][2]*vx[2] +
            v[3][0]*xx[0] - v[3][1]*xx[1] + v[3][2]*xx[2];

	d = x[3][0]*xx[0] - x[3][1]*xx[1] + x[3][2]*xx[2]; 
	
    //solve equation using method from "Art of Scientific Computing"
	//transform equation to t^3+at^2+bt+c = 0
	if (fabs(a) > MACH_EPS)
    {
	    b /= a; c /= a; d /= a;
	    a = b; b = c; c = d;
	    double Q, R, theta;
	    double Q3, R2;
	    Q = (a*a-3.0*b)/9.0;
	    R = (2.0*a*a*a-9.0*a*b+27.0*c)/54.0;
	    Q3 = Q*Q*Q;
	    R2 = R*R;

	    if (R2 < Q3)
        {
	        double Qsqrt = sqrt(Q);
            theta = acos(R/sqrt(Q3));
            roots[0] = -2.0*Qsqrt*cos(theta/3.0)-a/3.0;
            roots[1] = -2.0*Qsqrt*cos((theta+2.0*M_PI)/3.0)-a/3.0;
            roots[2] = -2.0*Qsqrt*cos((theta-2.0*M_PI)/3.0)-a/3.0;
        }
	    else
        {
            double A, B;
            double sgn = (R > 0) ? 1.0 : -1.0;
            A = -sgn*pow(fabs(R)+sqrt(R2-Q3),1.0/3.0);
            B = (fabs(A) < ROUND_EPS) ? 0.0 : Q/A;
            roots[0] = (A+B)-a/3.0;
            if (fabs(A-B) < ROUND_EPS)
                roots[1] = roots[2] = -0.5*(A+B)-a/3.0; //multiple roots
	    }
	}
	else
    {
		a = b; b = c; c = d;
	   	double delta = b*b-4.0*a*c;
	   	if (fabs(a) > ROUND_EPS && delta > 0)
        {
		    double delta_sqrt = sqrt(delta);
		    roots[0] = (-b+delta_sqrt)/(2.0*a);
            roots[1] = (-b-delta_sqrt)/(2.0*a);
	   	}
		else if (fabs(a) < ROUND_EPS && fabs(b) > ROUND_EPS)
		{
		    roots[0] = -c/b;
        }
	}

	//elimiate invalid roots;
	for (int i = 0; i < 3; ++i)
    {
        roots[i] -= MACH_EPS;
        if (roots[i] < 0.0 || roots[i] > dt)
            roots[i] = -1;
	}

	//sort the roots
	if (roots[0] > roots[1])
	    std::swap(roots[0], roots[1]);
	if (roots[0] > roots[2])
	    std::swap(roots[0], roots[2]);
	if (roots[1] > roots[2])
	    std::swap(roots[1], roots[2]);

	if (roots[0] > MACH_EPS || roots[1] > MACH_EPS || roots[2] > MACH_EPS)
	    return true;
	else
	    return false;
}

bool TriToBond(const TRI* tri,const BOND* bd)
{
	bool status = false;

    double tol = CollisionSolver3d::getFabricThickness();
    STATE* sl = (STATE*)left_state(bd->start);
    if (sl->is_stringpt)
        tol = CollisionSolver3d::getStringThickness();

	POINT* pts[4];
	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	pts[3] = bd->start;
	if (PointToTri(pts,tol))
        status = true;

	pts[3] = bd->end;
    if (PointToTri(pts,tol))
        status = true;
	
	pts[2] = bd->start;
    pts[3] = bd->end;
    for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
	    if (EdgeToEdge(pts,tol))
            status = true;
	}

	return status;
}

bool BondToBond(const BOND* b1, const BOND* b2)
{
	POINT* pts[4];

	pts[0] = b1->start;
	pts[1] = b1->end;
	pts[2] = b2->start;
	pts[3] = b2->end;

    double tol = CollisionSolver3d::getStringThickness();

	bool status = false;
	if (EdgeToEdge(pts,tol))
		status = true;

	return status;
}

bool TriToTri(const TRI* tri1, const TRI* tri2)
{
	bool status = false;
    double tol = CollisionSolver3d::getFabricThickness();

	POINT* pts[4];

	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
	{
	    const TRI* tmp_tri1 = (k == 0) ? tri1 : tri2;
	    const TRI* tmp_tri2 = (k == 0) ? tri2 : tri1;
	    
        for (int j = 0; j < 3; ++j)
            pts[j] = Point_of_tri(tmp_tri2)[j];
	    pts[3] = Point_of_tri(tmp_tri1)[i];
	
	    if (PointToTri(pts,tol))
            status = true;
	}

	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri1)[i];
	    pts[1] = Point_of_tri(tri1)[(i+1)%3];
        for (int j = 0; j < 3; ++j)
        {
            pts[2] = Point_of_tri(tri2)[j];
            pts[3] = Point_of_tri(tri2)[(j+1)%3];
            
            if (EdgeToEdge(pts,tol))
                status = true;
	    }  
	}

	return status;
}

//For details of this implementation see:
//http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
//
//Note that mstate has default value of MotionState::STATIC,
//and root has default value of -1.0
static bool EdgeToEdge(
        POINT** pts,
        double tol,
        MotionState mstate,
        double root)
{
	double x12[3], x34[3], x31[3];
	Pts2Vec(pts[0],pts[1],x12);    
	Pts2Vec(pts[2],pts[3],x34);
	Pts2Vec(pts[2],pts[0],x31);

    //Matrix entries
    double a = Dot3d(x12,x12);
    double b = Dot3d(x12,x34);
    double c = Dot3d(x34,x34);

    //RHS
    double d = Dot3d(x12,x31);
    double e = Dot3d(x34,x31);
	
    //Matrix Determinant
    double D = fabs(a*c - b*b);

    //Solution, and solution numerators and denominators
    double sC = 0;  double sN = 0;  double sD = D;    
    double tC = 0;  double tN = 0;  double tD = D;    
    
    //The solution is: sC = sN/sD and tC = tN/tD (Cramer's Rule).
    //Seperation of the numerator and denominator allows us to
    //efficiently analyze the boundary of the constrained domain,
    //(s,t) in [0,1]x[0,1], when the global minimum does not occur
    //within this region of parameter space.

    double vec[3];
	Cross3d(x12,x34,vec);

    if (D < ROUND_EPS || Mag3d(vec) < ROUND_EPS)
    {
        //Lines containing the edges are nearly parallel.
        //Setting sC = 0, and solving for tC yields tC = e/c.
        double sN = 0.0;
        double sD = 1.0;
        double tN = e;
        double tD = c;
    }
    else
    {
        //Compute the closest pair of points on the infinite lines.
        sN = b*e - c*d;
        tN = a*e - b*d;
        
        if( sN < 0.0 )
        {
            //Implies sC < 0 and the s = 0 edge is visible.
            sN = 0.0;
            tN = e;
            tD = c;

        }
        else if( sN > sD )
        {
            //Implies sC > 1 and the s = 1 edge is visible.
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if( tN < 0.0 )
    {
        //Implies tC < 0 and the t = 0 edge visible.
        tN = 0.0;
        
        //Recompute sC for this edge
        if (-1.0*d < 0.0)
            sN = 0.0;
        else if (-1.0*d > a)
            sN = sD;
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD)
    {
        //Implies tC > 1 and the t = 1 edge visible.
        tN = tD;
        
        //Recompute sC for this edge
        if ((b - d)  < 0.0)
            sN = 0.0;
        else if ((b - d) > a)
            sN = sD;
        else
        {
            sN = b - d;
            sD = a;
        }
    }

    //Compute the closest pair of points
    if (sN == sD)
        sC = 1.0;
    else
        sC = fabs(sN) < ROUND_EPS ? 0.0 : sN/sD;

    if (tN == tD)
        tC = 1.0;
    else
        tC = fabs(tN) < ROUND_EPS ? 0.0 : tN/tD;
    
	double x13[3];
    Pts2Vec(pts[0],pts[2],x13);
    
    scalarMult(tC,x34,x34);
    addVec(x13,x34,vec);

    scalarMult(sC,x12,x12);
    minusVec(vec,x12,vec);

    double dist = Mag3d(vec);
    if (dist > tol) return false;

    //TODO: handle another way -- restart with smaller dt for example
    if (dist > 0)
        scalarMult(1.0/dist,vec,vec);
    else if (dist == 0 && mstate == MotionState::STATIC)
    {
        printf("\n\tEdgeToEdge() ERROR: dist == 0 in proximity detection\n");
        printf("\t vec = %g %g %g",vec[0],vec[1],vec[2]);
        printf(",\t dist = %g\n\n",dist);
        printf("\tPOINTS:\n");
        for (int i = 0; i < 4; ++i)
        {
            double* coords = Coords(pts[i]);
            printf("\t\tpts[%d]: %g %g %g\t Gindex = %ld\n",
                    i,coords[0],coords[1],coords[2],Gindex(pts[i]));
        }

        //For debugging, comment out clean_up() below to print all
        //violating edge points.
        static int ecount = 0;
        std::string fname = CollisionSolver3d::getOutputDirectory();
        fname += "/EdgeToEdge_error-" + std::to_string(ecount);
        ecount++;

        std::vector<POINT*> edge_pts(pts,pts+4);
        vtk_write_pointset(edge_pts,fname,ERROR);

        LOC(); clean_up(ERROR);
    }

    //TODO: ALLOW IMPACT ZONES FOR STRING-STRING POINTS FOR NOW.
    bool string_string = false;
    
    /*
    STATE* s0 = (STATE*)left_state(pts[0]);
    STATE* s2 = (STATE*)left_state(pts[2]);
    if (s0->is_stringpt && s2->is_stringpt)
        string_string = true;
    */

    bool rigid_body_point = false;
    if (isRigidBody(pts[0]) || isRigidBody(pts[3]))
        rigid_body_point = true;

    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    //if (!is_detImpZone || string_string)
    if (!is_detImpZone || rigid_body_point || string_string)
    {
        double dt = root;
        if (mstate == MotionState::STATIC)
            dt = CollisionSolver3d::getTimeStepSize();
        EdgeToEdgeImpulse(pts,vec,sC,tC,dist,tol,mstate,dt);
    }

	return true;
}

//The "normal" vector, nor, points from the
//closest point on edge01 to the closest point on edge23 
static void EdgeToEdgeImpulse(
        POINT** pts,
        double* nor,
        double a,
        double b,
        double dist,
        double tol,
        MotionState mstate,
        double dt)
{
    if (debugging("collision"))
        CollisionSolver3d::edg_to_edg++;

	double v_rel[3] = {0.0, 0.0, 0.0};
    double vn = 0.0;
    double vt = 0.0;

    double rigid_impulse[2] = {0.0};
	double inelastic_impulse[2] = {0.0};
    double elastic_impulse[2] = {0.0};
	
	double wa[2] = {1.0 - a, a};
    double wb[2] = {1.0 - b, b};
    double wab[4] = {wa[0], wa[1], wb[0], wb[1]};

	double h = CollisionSolver3d::getFabricThickness();
	double k = CollisionSolver3d::getFabricSpringConstant();
	double m = CollisionSolver3d::getFabricPointMass();
    double mu = CollisionSolver3d::getFabricFrictionConstant(); 
    double overlap_coef = 0.1;
	
	STATE *sl[4];
	for (int i = 0; i < 4; ++i)
	    sl[i] = (STATE*)left_state(pts[i]);

    if (sl[0]->is_stringpt || sl[2]->is_stringpt)
    {
        h = CollisionSolver3d::getStringThickness();
        k = CollisionSolver3d::getStringSpringConstant();
        m = CollisionSolver3d::getStringPointMass();
        mu = CollisionSolver3d::getStringFrictionConstant();
    }

    double overlap = h - dist;
    //double overlap = tol - dist;

	//apply impulses to the average velocity (linear trajectory)
	for (int j = 0; j < 3; ++j)
	{
	    v_rel[j]  = (1.0-b) * sl[2]->avgVel[j] + b * sl[3]->avgVel[j];
	    v_rel[j] -= (1.0-a) * sl[0]->avgVel[j] + a * sl[1]->avgVel[j];
	}
    double mag_vrel = Mag3d(v_rel);
	
    vn = Dot3d(v_rel, nor);
	if (Dot3d(v_rel, v_rel) > sqr(vn))
	    vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	else
	    vt = 0.0;
    
    if (debugging("CollisionImpulse"))
    {
        printf("vn = %g\n",vn);
        printf("vn * dt = %g  %s  %g * overlap = %g\n",
                vn*dt, (vn*dt < overlap_coef*overlap) ? "<" : ">",
                overlap_coef, overlap_coef*overlap);
    }
	
    //Edges are approaching each other (vn < 0.0):
    //      apply inelastic impulse
    
    //Edges are seperating from each other (vn > 0.0):
    //      apply elastic impulse

    if (mstate == MotionState::STATIC)
    {
        // May apply both for repulsion.
        // Zero the normal component of relative velocity with inelastic impulse.
        if (vn < 0.0)
            EdgeToEdgeInelasticImpulse(vn,pts,inelastic_impulse,rigid_impulse,wab);
        if (vn * dt <  overlap_coef * overlap)
            EdgeToEdgeElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }
    else if (mstate == MotionState::MOVING)
    {
        // Apply one or the other for collision, NOT BOTH.
        // Zero the relative velocity with inelastic impulse.
        if (vn < 0.0)
            EdgeToEdgeInelasticImpulse(vn,pts,inelastic_impulse,rigid_impulse,wab);
        else if (vn * dt <  overlap_coef * overlap)
            EdgeToEdgeElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }
    else //mstate == MOTION_TYPE::POSTCOLLISION
    {
        if (vn * dt <  overlap_coef * overlap)
            EdgeToEdgeElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }

    
    double impulse[2];
    double m_impulse[2];
    double f_impulse[2];

    for (int i = 0; i < 2; ++i)
    {
        impulse[i] = inelastic_impulse[i] + elastic_impulse[i];

        /*
        impulse[i] = inelastic_impulse[i];
        if (mstate == MotionState::MOVING)
            impulse[i] += elastic_impulse[i];
        */

        if (wab[0] + wab[1] < MACH_EPS || wab[2] + wab[3] < MACH_EPS)
        {
            m_impulse[i] = impulse[i];
            f_impulse[i] = impulse[i];
        }
        else
        {
            double wabs_sqr = sqr(wab[0]) + sqr(wab[1])
                            + sqr(wab[2]) + sqr(wab[3]);
        
            m_impulse[i] = 2.0*impulse[i]/wabs_sqr;
            f_impulse[i] = 2.0*impulse[i]/wabs_sqr;
        }
    }

    ////////////////////////////////////////////////////////////////////
    if (debugging("CollisionImpulse"))
    {
        if (fabs(m_impulse[0] + m_impulse[1]) > 0.0)
        {
            printf("\tEdgeToEdgeImpulse():\n");
            printf("dt = %e, step_dt = %e\n",dt,CollisionSolver3d::getTimeStepSize());
            printf("h = %e, dist = %e, overlap = %e\n",h,dist,overlap);
            printf("inelastic_impulse[0] = %g, inelastic_impulse[1] = %g\n",
                    inelastic_impulse[0],inelastic_impulse[1]);
            printf("elastic_impulse[0] = %g, elastic_impulse[1] = %g\n",
                    elastic_impulse[0],elastic_impulse[1]);
            printf("impulse[0] = %g, impulse[1] = %g\n",impulse[0],impulse[1]);
            printf("m_impulse[0] = %g, m_impulse[1] = %g\n",m_impulse[0],m_impulse[1]);
            printf("k = %g, m = %g, mu = %g\n",k,m,mu);
            printf("vn = %g, vt = %g\n",vn,vt);
            printf("v_rel = %g %g %g\n",v_rel[0],v_rel[1],v_rel[2]);
            printf("nor = %g %g %g\n\n",nor[0],nor[1],nor[2]);
            
            printf("a = %g, b = %g\n",a,b);
            printf("wa[0] = %g, wa[1] =  %g, wb[0] = %g, wb[1] = %g\n\n",
                    wa[0],wa[1],wb[0],wb[1]);
            
            printf("x_old:\n");
            for (int i = 0; i < 4; ++i)
            {
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%g %g %g   Gindex(pts[%d]) = %lu\n",
                        sl1->x_old[0],sl1->x_old[1],sl1->x_old[2],
                        i,Gindex(pts[i]));
            }
            printf("x_new:\n");
            for (int i = 0; i < 4; ++i)
            {
                printf("%g %g %g\n",
                Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
            }
            printf("avgVel:\n");
            for (int i = 0; i < 4; ++i)
            {
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%g %g %g\n",sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
            }
            printf("\n");
        }
    }
    ////////////////////////////////////////////////////////////////////

	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
        //TODO: current nor vector not correct for rigid-rigid collision
	    if (isMovableRigidBody(pts[0]))
            SpreadImpactZoneImpulse(pts[0], -1.0*rigid_impulse[0], nor);
        if (isMovableRigidBody(pts[2]))
            SpreadImpactZoneImpulse(pts[2], rigid_impulse[1], nor);
	    return;
	}
	
    double max_friction = 0.5*vt;
    if ((isRigidBody(pts[0]) && isRigidBody(pts[1])) ||
        (isRigidBody(pts[2]) && isRigidBody(pts[3])))
    {
        max_friction = vt;
    }

    std::vector<double> W = {-wab[0],-wab[1],wab[2],wab[3]};
    std::vector<double> M = {m_impulse[0],m_impulse[0],
                             m_impulse[1],m_impulse[1]};
    std::vector<double> F = {f_impulse[0],f_impulse[0],
                             f_impulse[1],f_impulse[1]};
    std::vector<double> R = {rigid_impulse[0],rigid_impulse[0],
                             rigid_impulse[1],rigid_impulse[1]};

    for (int i = 0; i < 4; ++i)
    {
        if (!isStaticRigidBody(pts[i]))
        {
            double t_impulse = M[i];
            if (isMovableRigidBody(pts[i]))
                t_impulse = R[i];
            
            for (int j = 0; j < 3; ++j)
                sl[i]->collsnImpulse[j] += W[i]*t_impulse*nor[j];
       
            if (mstate == MotionState::STATIC)
            {
                double friction_impulse = F[i];
                if (fabs(vt) > ROUND_EPS)
                {
                    double delta_vt = max_friction;
                    if (fabs(mu*friction_impulse) < max_friction)
                        delta_vt = fabs(mu*friction_impulse);

                    for (int j = 0; j < 3; ++j)
                        sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
                }
            }
                
            sl[i]->collsn_num++;
        }
        else
        {
            for (int j = 0; j < 3; ++j)
            {
                sl[i]->friction[j] = 0.0;
                sl[i]->collsnImpulse[j] = 0.0;
            }
        }
    }

	if (debugging("CollisionImpulse"))
    {
	    for (int i = 0; i < 4; ++i)
        {
            printf("pt[%d]:",i);
            printf(" collsnImpulse = [%g %g %g]", sl[i]->collsnImpulse[0],
                    sl[i]->collsnImpulse[1],sl[i]->collsnImpulse[2]);
            printf(" friction = [%g %g %g]\n",sl[i]->friction[0],
                    sl[i]->friction[1],sl[i]->friction[2]);
        }
        printf("\n");
        fflush(stdout);
    }

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 3; ++j)
    {
	    if (std::isnan(sl[i]->collsnImpulse[j]) ||
            std::isinf(sl[i]->collsnImpulse[j]))
        {
		    printf("EdgeToEdge: sl[%d]->collsnImpulse[%d] = nan\n",i,j);
		    printf("a b = %g %g, nor = [%g %g %g], dist = %g\n",
                    a,b,nor[0],nor[1],nor[2],dist);
	        LOC(); clean_up(ERROR);
	    }
	}
}

static void EdgeToEdgeInelasticImpulse(
        double vn,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double* W)
{
    /*
    double Ip0 = isStaticRigidBody(pts[0]) ? 0.0 : W[0]*W[0]/getPointMass(pts[0]);
    double Ip1 = isStaticRigidBody(pts[1]) ? 0.0 : W[1]*W[1]/getPointMass(pts[1]);
    double Ip2 = isStaticRigidBody(pts[2]) ? 0.0 : W[2]*W[2]/getPointMass(pts[2]);
    double Ip3 = isStaticRigidBody(pts[3]) ? 0.0 : W[3]*W[3]/getPointMass(pts[3]);
    double I = fabs(vn)/(Ip0 + Ip1 +Ip2 +Ip3);
    */

    vn = fabs(vn);

    if ((isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])) ||
        (isStaticRigidBody(pts[2]) && isStaticRigidBody(pts[3])))
    {
        impulse[0] = vn;
        impulse[1] = vn;
        rigid_impulse[0] = vn;
        rigid_impulse[1] = vn;
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1])
            && isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        double m1 = total_mass(pts[0]->hs);
        double m2 = total_mass(pts[2]->hs);
        rigid_impulse[0] = vn * m2 / (m1 + m2);
        rigid_impulse[1] = vn * m1 / (m1 + m2);
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]))
    {
        impulse[1] = 0.5 * vn; 
        rigid_impulse[0] = 0.5 * vn;
    }
    else if (isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        impulse[0] = 0.5 * vn;
        rigid_impulse[1] = 0.5 * vn;
    }
    else
    {
        impulse[0] = 0.5 * vn;
        impulse[1] = 0.5 * vn;
    }

    //multiply inelastic impulses by relaxation parameter in [0,1]
        //impulse[0] *= 0.25;
        //impulse[1] *= 0.25;

    if (isStaticRigidBody(pts[0])) W[0] = 0.0;
    if (isStaticRigidBody(pts[1])) W[1] = 0.0;
    if (isStaticRigidBody(pts[2])) W[2] = 0.0;
    if (isStaticRigidBody(pts[3])) W[3] = 0.0;
}

static void EdgeToEdgeElasticImpulse(
        double vn,
        double overlap_coef,
        double overlap,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double dt,
        double m,
        double k)
{
    if (isRigidBody(pts[0]) && isRigidBody(pts[1]) &&
        isRigidBody(pts[2]) && isRigidBody(pts[3]))
    {
        double cor = CollisionSolver3d::getRestitutionCoef();
        rigid_impulse[0] *= 1.0 + cor;
        rigid_impulse[1] *= 1.0 + cor;
    }
    else
    {
        if (debugging("CollisionImpulse"))
        {
            printf("dt*k*overlap/m = %g,  %g*overlap/dt - vn = %g\n",
                    dt*k*overlap/m, overlap_coef, overlap_coef*overlap/dt - vn);
        }

        double I = std::min(dt*k*overlap/m, (overlap_coef*overlap/dt - vn));

        impulse[0] += 0.5*I;
        impulse[1] += 0.5*I;
        rigid_impulse[0] += 0.5*I;
        rigid_impulse[1] += 0.5*I;
    }
}

//Note that mstate has default value of MotionState::STATIC,
//and root has default value of -1.0
static bool PointToTri(
        POINT** pts,
        double tol,
        MotionState mstate,
        double root)
{
/*	x1
 *  	/\     x4 *
 *     /  \
 * x2 /____\ x3
 *
 * solve equation
 * x13*x13*w1 + x13*x23*w2 = x13*x43
 * x13*x23*w1 + x23*x23*w2 = x23*x43
 */
    double nor[3];
	double w[3] = {0.0};
	double x13[3], x23[3], x43[3], x34[3];
    double tri_nor[3] = {0.0};

	Pts2Vec(pts[0],pts[2],x13);
	Pts2Vec(pts[1],pts[2],x23);
	Pts2Vec(pts[3],pts[2],x43);
    scalarMult(-1.0,x43,x34);
	
    //unit normal vector of the plane of the triangle
    Cross3d(x13,x23,tri_nor);
    double mag_tnor = Mag3d(tri_nor);
    scalarMult(1.0/mag_tnor,tri_nor,tri_nor);

    //correct the triangle's normal direction to point to same
    //side as the point (not used right now, but may need at some
    //for detecting/correcting interpenetration etc.)
    double side = Dot3d(x34,tri_nor);
    if (side < 0.0)
    {
        scalarMult(-1.0,tri_nor,tri_nor);
    }
	
    double dist = fabs(side);
    if (dist > tol) return false;
	
	double det = Dot3d(x13,x13)*Dot3d(x23,x23) - Dot3d(x13,x23)*Dot3d(x13,x23);
	if (fabs(det) < MACH_EPS)
    {   
        printf("\n\tPointToTri() WARNING: degenerate TRI detected,\n \
                \t\t\t (fabs(det) < MACH_EPS)\n\n");
        
        printf("\tPOINTS:\n");
        for (int i = 0; i < 4; ++i)
        {
            double* coords = Coords(pts[i]);
            printf("\t\tpts[%d]: %g %g %g\t Gindex = %ld\n",
                    i,coords[0],coords[1],coords[2],Gindex(pts[i]));
        }

        //For debugging, comment out clean_up() below to print all instances.
        static int ecount = 0;
        std::string dname = CollisionSolver3d::getOutputDirectory();
        std::string fname = dname + "/PointToTri_error-" + std::to_string(ecount);
        ecount++;

        std::vector<POINT*> pt2tri_pts(pts,pts+4);
        vtk_write_pointset(pt2tri_pts,fname,ERROR);

        double BBL[3], BBU[3];
        set_point_list_bounding_box(pts,4,BBL,BBU,NO,YES);
        gview_plot_vertices(dname.c_str(),"PointToTri_error",pts,4,BBL,BBU);

        CollisionSolver3d::saveFront();
        CollisionSolver3d::drawFront();
        LOC(); clean_up(ERROR);
        //return false;
	}
	else
    {
	    w[0] = (Dot3d(x23,x23)*Dot3d(x13,x43)-Dot3d(x13,x23)*Dot3d(x23,x43))/det;
	    w[1] = (Dot3d(x13,x13)*Dot3d(x23,x43)-Dot3d(x13,x23)*Dot3d(x13,x43))/det;
	    w[2] = 1.0 - w[0] - w[1];
	    
        //TODO: Experiment with this alternate definition for the
        //      characteristic length of a triangle.
        /*
        double c_len = 0.0;	
        for (int i = 0; i < 3; ++i)
        {
            double tmp_dist = distance_between_positions(Coords(pts[i]),
                                    Coords(pts[(i+1)%3]),3);
            if (tmp_dist > c_len)
                c_len = tmp_dist;
        }
        double eps = tol/c_len;
        */

        double tri_area = 0.5*mag_tnor;
        double eps = tol/sqrt(tri_area);
        for (int i = 0; i < 3; ++i)
        {
            if (w[i] < -1.0*eps || w[i] > 1.0 + eps)
                return false;
        }
    }

    bool rigid_body_point = false;
    if (isRigidBody(pts[0]) || isRigidBody(pts[3]))
        rigid_body_point = true;

    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    //if (!is_detImpZone)
    if (!is_detImpZone || rigid_body_point)
    {
        double dt = root;
        if (mstate == MotionState::STATIC)
            dt = CollisionSolver3d::getTimeStepSize();
        PointToTriImpulse(pts,tri_nor,w,dist,tol,mstate,dt);
    }

	return true;
}

//The normal vector, nor, is the normal vector of the
//triangle tri012 pointing to the side of p3
static void PointToTriImpulse(
        POINT** pts,
        double* nor,
        double* w,
        double dist,
        double tol,
        MotionState mstate,
        double dt)
{
    if (debugging("collision"))
        CollisionSolver3d::pt_to_tri++;
    
    double vn = 0.0;
    double vt = 0.0;
    double v_rel[3] = {0.0};

	double rigid_impulse[2] = {0.0};
	double inelastic_impulse[2] = {0.0};
    double elastic_impulse[2] = {0.0};
    
    double sum_w = 0.0;

	double h = CollisionSolver3d::getFabricThickness();
	double k = CollisionSolver3d::getFabricSpringConstant();
	double m = CollisionSolver3d::getFabricPointMass();
	double mu = CollisionSolver3d::getFabricFrictionConstant();
   
	STATE *sl[4];
	for (int i = 0; i < 4; ++i)
	    sl[i] = (STATE*)left_state(pts[i]);

    double overlap_coef = 0.1;
    /*
    if (sl[3]->is_stringpt)
    {
        h = CollisionSolver3d::getStringThickness();
        k = CollisionSolver3d::getStringSpringConstant();
        m = CollisionSolver3d::getStringPointMass();
        mu = CollisionSolver3d::getStringFrictionConstant();
        overlap_coef = 0.1;
    }
    */
	
    double overlap = h - dist;
    //double overlap = tol - dist;

	//apply impulses to the average (linear trajectory) velocity
	for (int i = 0; i < 3; ++i)
	{
	    v_rel[i] = sl[3]->avgVel[i];
	    for (int j = 0; j < 3; ++j)
            v_rel[i] -= w[j] * sl[j]->avgVel[i];
	}
    double mag_vrel = Mag3d(v_rel);

	vn = Dot3d(v_rel, nor);
	if (Dot3d(v_rel, v_rel) > sqr(vn))
	    vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	else
	    vt = 0.0;

    if (debugging("CollisionImpulse"))
    {
        printf("vn = %g\n",vn);
        printf("vn * dt = %g  %s  %g * overlap = %g\n",
                vn*dt, (vn*dt < overlap_coef*overlap) ? "<" : ">",
                overlap_coef, overlap_coef*overlap);
    }
	
    if (mstate == MotionState::STATIC)
    {
        // May apply both for repulsion.
        // Zero the normal component of relative velocity with inelastic impulse.
        if (vn < 0.0)
            PointToTriInelasticImpulse(vn,pts,inelastic_impulse,rigid_impulse,w,&sum_w);
        if (vn * dt < overlap_coef*overlap)
            PointToTriElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }
    else if (mstate == MotionState::MOVING)
    {
        // Apply one or the other for collision, NOT BOTH.
        // Zero the relative velocity with inelastic impulse.
        if (vn < 0.0)
            PointToTriInelasticImpulse(vn,pts,inelastic_impulse,rigid_impulse,w,&sum_w);
        else if (vn*dt < overlap_coef*overlap)
            PointToTriElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }
    else //mstate == MOTION_TYPE::POSTCOLLISION
    {
        if (vn * dt < overlap_coef*overlap)
            PointToTriElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }

    double impulse[2];
    double m_impulse[2];
    double f_impulse[2];

    for (int i = 0; i < 2; ++i)
    {
        impulse[i] = inelastic_impulse[i] + elastic_impulse[i];
        
        /*
        impulse[i] = inelastic_impulse[i];
        if (mstate == MotionState::MOVING)
            impulse[i] += elastic_impulse[i];
        */

        if (fabs(sum_w) < MACH_EPS)
        {
            m_impulse[i] = impulse[i];
            f_impulse[i] = impulse[i];
        }
        else
        {
            m_impulse[i] = 2.0*impulse[i]/(1.0 + Dot3d(w, w));
            f_impulse[i] = 2.0*impulse[i]/(1.0 + Dot3d(w, w));
        }
    }

    ////////////////////////////////////////////////////////////////////
    if (debugging("CollisionImpulse"))
    {
        if (fabs(m_impulse[0] + m_impulse[1]) > 0.0)
        {
            printf("\tPointToTriImpulse():\n");
            printf("dt = %e, step_dt = %e\n",dt,CollisionSolver3d::getTimeStepSize());
            printf("h = %e, dist = %e, overlap = %e\n",h,dist,overlap);
            printf("inelastic_impulse[0] = %g, inelastic_impulse[1] = %g\n",
                    inelastic_impulse[0],inelastic_impulse[1]);
            printf("elastic_impulse[0] = %g, elastic_impulse[1] = %g\n",
                    elastic_impulse[0],elastic_impulse[1]);
            printf("impulse[0] = %g, impulse[1] = %g\n",impulse[0],impulse[1]);
            printf("m_impulse[0] = %g, m_impulse[1] = %g\n",m_impulse[0],m_impulse[1]);
            printf("k = %g, m = %g, mu = %g\n",k,m,mu);
            printf("vn = %g, vt = %g\n",vn,vt);
            printf("v_rel = %g %g %g\n",v_rel[0],v_rel[1],v_rel[2]);
            printf("nor = %g %g %g\n\n",nor[0],nor[1],nor[2]);

            printf("w[0] = %g, w[1] = %g, w[2] = %g\n\n",w[0],w[1],w[2]);
            
            printf("x_old:\n");
            for (int i = 0; i < 4; ++i)
            {
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%g %g %g   Gindex(pts[%d]) = %lu\n",
                        sl1->x_old[0],sl1->x_old[1],
                        sl1->x_old[2],i,Gindex(pts[i]));
            }
            printf("x_new:\n");
            for (int i = 0; i < 4; ++i)
            {
                printf("%g %g %g\n",
                        Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
            }
            printf("avgVel:\n");
            for (int i = 0; i < 4; ++i)
            {
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%g %g %g\n",
                        sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
            }
            printf("\n");
            fflush(stdout);
        }
    }
    ////////////////////////////////////////////////////////////////////

	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
        //TODO: current nor vector not correct for rigid-rigid collision
	    if (isMovableRigidBody(pts[0]))
            SpreadImpactZoneImpulse(pts[0], -1.0*rigid_impulse[0], nor);
	    if (isMovableRigidBody(pts[3]))
            SpreadImpactZoneImpulse(pts[3], rigid_impulse[1], nor);
	    return;
	}

    double max_friction = 0.5*vt;
    if (isRigidBody(pts[3]) ||
       (isRigidBody(pts[0]) && isRigidBody(pts[1]) && isRigidBody(pts[2])))
    {
        max_friction = vt;
    }

    std::vector<double> W = {-w[0],-w[1],-w[2],1.0};
    std::vector<double> M = {m_impulse[0],m_impulse[0],
                             m_impulse[0],m_impulse[1]};
    std::vector<double> F = {f_impulse[0],f_impulse[0],
                             f_impulse[0],f_impulse[1]};
    std::vector<double> R = {rigid_impulse[0],rigid_impulse[0],
                             rigid_impulse[0],rigid_impulse[1]};

	for (int i = 0; i < 4; ++i)
	{
        if (!isStaticRigidBody(pts[i]))
        {
            double t_impulse = M[i];
            if (isMovableRigidBody(pts[i]))
                t_impulse = R[i];
            
            for (int j = 0; j < 3; ++j)
                sl[i]->collsnImpulse[j] += W[i]*t_impulse*nor[j];
            
            if (mstate == MotionState::STATIC)
            {
                double friction_impulse = F[i];
                if (fabs(vt) > ROUND_EPS)
                {
                    double delta_vt = max_friction;
                    if (fabs(mu*friction_impulse) < max_friction)
                        delta_vt = fabs(mu*friction_impulse);
                    
                    for (int j = 0; j < 3; ++j)
                        sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
                }
            }
            
            sl[i]->collsn_num++;
        }
        else
        {
            for (int j = 0; j < 3; ++j)
            {
                sl[i]->friction[j] = 0.0;
                sl[i]->collsnImpulse[j] = 0.0;
            }
        }
	}

	if (debugging("CollisionImpulse"))
    {
	    for (int i = 0; i < 4; ++i)
        {
            printf("pt[%d]:",i);
            printf(" collsnImpulse = [%g %g %g]", sl[i]->collsnImpulse[0],
                    sl[i]->collsnImpulse[1],sl[i]->collsnImpulse[2]);
            printf(" friction = [%g %g %g]\n",sl[i]->friction[0],
                    sl[i]->friction[1],sl[i]->friction[2]);
        }
        printf("\n");
        fflush(stdout);
    }

	for (int kk = 0; kk < 4; kk++)
	for (int j = 0; j < 3; ++j)
    {
        if (std::isnan(sl[kk]->collsnImpulse[j]) ||
		    std::isinf(sl[kk]->collsnImpulse[j]))
        {
            printf("PointToTri: sl[%d]->collsnImpulse[%d] = nan\n",kk,j);
            for (int i = 0; i < 4; ++i)
            {
                printf("points[%d] = %p\n",i,(void*)pts[i]);
                printf("coords = [%g %g %g]\n",Coords(pts[i])[0],
                    Coords(pts[i])[1],Coords(pts[i])[2]);
            }
            printf("w = [%g %g %g]\nnor = [%g %g %g]\ndist = %g\n",
            w[0],w[1],w[2],nor[0],nor[1],nor[2],dist);
            printf("v_rel = [%g %g %g]\n",v_rel[0],v_rel[1],v_rel[2]);
            LOC(); clean_up(ERROR);
	    }
	}
}

static void PointToTriInelasticImpulse(
        double vn,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double* w,
        double* sum_w)
{
    vn = fabs(vn);

    if (isStaticRigidBody(pts[3]) ||
       (isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])
        && isStaticRigidBody(pts[2])))
    {
        impulse[0] = vn;
        impulse[1] = vn;
        rigid_impulse[0] = vn;
        rigid_impulse[1] = vn;
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]) 
            && isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        double m1 = total_mass(pts[0]->hs);
        double m2 = total_mass(pts[3]->hs);
        rigid_impulse[0] = vn * m2 / (m1 + m2);
        rigid_impulse[1] = vn * m1 / (m1 + m2);
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1])
            && isMovableRigidBody(pts[2]))
    {
        impulse[1] = 0.5 * vn;
        rigid_impulse[0] = 0.5 * vn;
    }
    else if (isMovableRigidBody(pts[3]))
    {
        impulse[0] = 0.5 * vn;
        rigid_impulse[1] = 0.5 * vn;
    }
    else
    {
        impulse[0] = 0.5 * vn;
        impulse[1] = 0.5 * vn;
    }

    //multiply inelastic impulses by relaxation parameter in [0,1]
        //impulse[0] *= 0.25;
        //impulse[1] *= 0.25;
    
    for (int i = 0; i < 3; ++i)
    {
        if (isStaticRigidBody(pts[i])) w[i] = 0.0;
        *sum_w += w[i];
    }

    if (fabs(*sum_w) > MACH_EPS)
        scalarMult(1.0/(*sum_w),w,w);
}

static void PointToTriElasticImpulse(
        double vn,
        double overlap_coef,
        double overlap,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double dt,
        double m,
        double k)
{
    if (isRigidBody(pts[0]) && isRigidBody(pts[1]) &&
        isRigidBody(pts[2]) && isRigidBody(pts[3]))
    {
        double cor = CollisionSolver3d::getRestitutionCoef();
        rigid_impulse[0] *= 1.0 + cor;
        rigid_impulse[1] *= 1.0 + cor;
    }
    else
    {
        if (debugging("CollisionImpulse"))
        {
            printf("dt*k*overlap/m = %g,  %g*overlap/dt - vn = %g\n",
                    dt*k*overlap/m, overlap_coef, overlap_coef*overlap/dt - vn);
        }

        double I = std::min(dt*k*overlap/m, (overlap_coef*overlap/dt - vn));
        
        impulse[0] += 0.5*I;
        impulse[1] += 0.5*I;
        rigid_impulse[0] += 0.5*I;
        rigid_impulse[1] += 0.5*I;
    }
}

static double getPointMass(POINT* pt)
{
	double m = CollisionSolver3d::getFabricPointMass();
    STATE* sl = (STATE*)left_state(pt);
    if (sl->is_stringpt)
        m = CollisionSolver3d::getStringPointMass();
    return m;
}

static double getPointFrictionConstant(POINT* pt)
{
    double mu = CollisionSolver3d::getFabricFrictionConstant();
    STATE* sl = (STATE*)left_state(pt);
    if (sl->is_stringpt)
        mu = CollisionSolver3d::getStringFrictionConstant();
    return mu;
}

void CollisionSolver3d::printDebugVariable()
{
	std::cout << "Enter EdgeToEdge " << edg_to_edg << " times\n";
	std::cout << "Enter PointToTri " << pt_to_tri << " times\n";
	std::cout << "Enter isCoplanar " << is_coplanar << " times\n";
	
    moving_edg_to_edg = 0;
    moving_pt_to_tri = 0;
    is_coplanar = 0;
	edg_to_edg = 0;
    pt_to_tri = 0;
}

void printPointSetCollisionStats(POINT* pts[], int npts)
{
    printf("\n");
    for (int i = 0; i < npts; ++i)
    {
        printf("#%d ",i);
        printPointCollisionStats(pts[i]);
    }
    printf("\n");
}

void printPointCollisionStats(POINT* pt)
{
    STATE* sl = (STATE*)left_state(pt);
    printf("    Coords = [%f %f %f]\t",
        Coords(pt)[0],Coords(pt)[1],Coords(pt)[2]);
    printf(" Gindex = %d\n",Gindex(pt));
    printf("    avgVel = [%g %g %g]\t",
            sl->avgVel[0],sl->avgVel[1],sl->avgVel[2]);
    printf(" collsnImpulse = [%g %g %g]\n",
            sl->collsnImpulse[0],sl->collsnImpulse[1],sl->collsnImpulse[2]);
}

