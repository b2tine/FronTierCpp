#include <armadillo>
#include "collid.h"

static bool MovingEdgeToEdgeJac(POINT**);
static bool MovingEdgeToEdgeGS(POINT**);

static bool EdgeToEdge(POINT**, double,
        MotionState mstate = MotionState::STATIC, double root = -1.0);

static void EdgeToEdgeImpulse(POINT**,double*,double,double,double,MotionState,double);
static void EdgeToEdgeInelasticImpulse(double,POINT**,double*,double*,double*);
static void EdgeToEdgeElasticImpulse(double,double,double,POINT**,double*,double*,double,double,double);

static bool MovingPointToTriJac(POINT**);
static bool MovingPointToTriGS(POINT**);

static bool PointToTri(POINT**, double,
        MotionState mstate = MotionState::STATIC, double root = -1.0);

static void PointToTriImpulse(POINT**,double*,double*,double,MotionState,double);
static void PointToTriInelasticImpulse(double,POINT**,double*,double*,double*,double*);
static void PointToTriElasticImpulse(double,double,double,POINT**,double*,double*,double,double,double);

static bool isCoplanar(POINT**,double,double*);

// test function for creating impact zone for each movable RG
void CollisionSolver3d::createImpZoneForRG(const INTERFACE* intfc)
{
	SURFACE** s;
	TRI* tri;

	intfc_surface_loop(intfc, s)
	{
	    if (is_bdry(*s)) continue;
	    if (!isMovableRigidBody(Point_of_tri(first_tri(*s))[0])) continue;

        surf_tri_loop(*s, tri)
	    {
    		createImpZone(Point_of_tri(tri), 3, YES);
	    }
	}
}

//TODO: Optimize for when impact zone exists and we are
//      adding points to it. Should be able to  update
//      x_cm, v_cm, and avg_dt without looping through
//      every point of the zone.
void updateImpactListVelocity(POINT* head)
{
    POINT* p = nullptr;
	
    //compute impact zone's center of mass position and velocity
    int num_pts = 0;
    double totalmass = 0.0;
    double avg_dt = 0.0;

    double x_cm[3] = {0.0};
    double v_cm[3] = {0.0};
	
    p = head;
	while(p)
    {
		STATE* sl = (STATE*)left_state(p);
        avg_dt += sl->collsn_dt;
        sl->has_collsn = true;

        double m = CollisionSolver3d::getFabricPointMass();
        if (sl->is_stringpt)
            m = CollisionSolver3d::getStringPointMass();

        totalmass += m;
        for (int i = 0; i < 3; ++i)
        {
		    x_cm[i] += sl->x_old[i]*m; 
		    v_cm[i] += sl->avgVel[i]*m;
		}

        //Still need to mark sorted() for updateImpactZoneVelocityForRG()
        sorted(p) = YES;
        
        p = next_pt(p);
		num_pts++;
    }
	
    //TODO: Is this justified, or just use the full step dt?
    avg_dt += CollisionSolver3d::getTimeStepSize();
    avg_dt /= (double)num_pts;
    
    //temp debug
        //double dt = getTimeStepSize();
        //printf("avg_dt = %g,  dt = %g\n",avg_dt,dt);


	for (int i = 0; i < 3; ++i)
    {
	    x_cm[i] /= totalmass;
	    v_cm[i] /= totalmass;
	}

	//compute angular momentum
	double L[3] = {0.0};

    p = head;
	while(p)
    {
	    STATE* sl = (STATE*)left_state(p);
        double m = CollisionSolver3d::getFabricPointMass();
        if (sl->is_stringpt)
            m = CollisionSolver3d::getStringPointMass();
	    
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
            if (ars(i))
                ars(i) = 1.0/ars(i);
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
        
        //double dt = CollisionSolver3d::getTimeStepSize();
    
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
            x_new[i] = x_cm[i] + avg_dt*v_cm[i] + xF[i]
                       + cos(avg_dt*mag_w)*xR[i] + wxR[i];

		    sl->avgVel[i] = (x_new[i] - sl->x_old[i])/avg_dt;

            /*
            x_new[i] = x_cm[i] + dt*v_cm[i] + xF[i]
                       + cos(dt*mag_w)*xR[i] + wxR[i];
	
		    sl->avgVel[i] = (x_new[i] - sl->x_old[i])/dt;
            */
	    	
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

bool MovingTriToBond(const TRI* tri,const BOND* bd)
{
    bool status = false;

	POINT* pts[4];
	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	/* detect collision of start point of bond w.r.t to tri */
	pts[3] = bd->start;
    if (MovingPointToTriGS(pts))
        status = true;
    
    //if (status && is_detImpZone)
      //  createImpZone(pts,4);
	
    /* detect collision of end point of bond to w.r.t. tri */
	pts[3] = bd->end;
    if (MovingPointToTriGS(pts))
        status = true;

    //if (status && is_detImpZone)
      //  createImpZone(pts,4);
	
    /* detect collision of each of tri edge w.r.t to bond */
	pts[2] = bd->start;
	pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
        if (MovingEdgeToEdgeGS(pts))
            status = true;

        //if (status && is_detImpZone)
          //  createImpZone(pts,4);
	}

    return status;
}

bool MovingBondToBond(const BOND* b1, const BOND* b2)
{
	POINT* pts[4];

	pts[0] = b1->start;
	pts[1] = b1->end;
	pts[2] = b2->start;
	pts[3] = b2->end;

	bool status = false;
    if(MovingEdgeToEdgeGS(pts))
        status = true;

    //if (status && is_detImpZone)
      //  createImpZone(pts,4);
    
    return status;
}

bool MovingTriToTri(const TRI* a,const TRI* b)
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

        if(MovingPointToTriGS(pts))
            status = true;

        //if (status && is_detImpZone)
          //  createImpZone(pts,4);
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
		
            if(MovingEdgeToEdgeGS(pts))
                status = true;
                
            //if (status && is_detImpZone)
              //  createImpZone(pts,4);
	    }
    }

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
                    sl->has_collsn = true;
                }
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
    
    return status;
}

//gauss-seidel update
static bool MovingPointToTriGS(POINT* pts[])
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

    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    if (status && !is_detImpZone)
    {
        for (int j = 0; j < 4; ++j)
        {
            STATE* sl = (STATE*)left_state(pts[j]);
            if (sl->collsn_num > 0)
            {
                sl->has_collsn = true;
                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] += sl->collsnImpulse[k]/((double)sl->collsn_num);
                    sl->collsnImpulse[k] = 0.0;
                }
                sl->collsn_num = 0;
            }
        }
    }
    else if (status && is_detImpZone)
    {
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
            if (roots[i] < 0)
                continue;
                
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
                    sl->has_collsn = true;
                }
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
    
    return status;
}

//gauss-seidel update
static bool MovingEdgeToEdgeGS(POINT* pts[])
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
            if (roots[i] < 0)
                continue;
                
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

    for (int j = 0; j < 4; ++j)
    {
        STATE* sl = (STATE*)left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k];
    }

    //TODO: ALLOW IMPACT ZONES FOR STRING-STRING INTERACTIONS FOR NOW.
    bool string_string = false;

    /*
    //No Impact Zones for string-string interactions
    bool string_string = false;
    if (s0->is_stringpt && s2->is_stringpt)
        string_string = true;
    */

	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    if (status && (!is_detImpZone || string_string))
    {
        for (int j = 0; j < 4; ++j)
        {
            STATE* sl = (STATE*)left_state(pts[j]);
            if (sl->collsn_num > 0)
            {
                sl->has_collsn = true;
                for (int k = 0; k < 3; ++k)
                {
                    sl->avgVel[k] += sl->collsnImpulse[k]/((double)sl->collsn_num);
                    sl->collsnImpulse[k] = 0.0;
                }
                sl->collsn_num = 0;
            }
        }
    }
    else if (status && is_detImpZone && !string_string)
    {
        createImpZone(pts,4);
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
            B = (fabs(A) < MACH_EPS) ? 0.0 : Q/A;
            roots[0] = (A+B)-a/3.0;
            if (fabs(A-B) < MACH_EPS)
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
        //TODO: necessary to subtract off MACH_EPS valid here?
        roots[i] -= MACH_EPS;
        if (roots[i] < 0 || roots[i] > dt)
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
    if (dist > tol)
        return false;

    //TODO: handle another way -- restart with smaller dt for example
    if (dist > 0)
        scalarMult(1.0/dist,vec,vec);
    else if (dist == 0 && mstate == MotionState::STATIC)
    {
        //printf("\n\tEdgeToEdge() WARNING: dist < 0\n\n");
        printf("\n\tEdgeToEdge() ERROR: dist == 0 in proximity detection\n");
        printf("\t vec = %g %g %g",vec[0],vec[1],vec[2]);
        printf(",\t dist = %g\n\n",dist);
        clean_up(ERROR);
    }

    //TODO: ALLOW IMPACT ZONES FOR STRING-STRING POINTS FOR NOW.
    bool string_string = false;
    
    /*
    STATE* s0 = (STATE*)left_state(pts[0]);
    STATE* s2 = (STATE*)left_state(pts[2]);
    if (s0->is_stringpt && s2->is_stringpt)
        string_string = true;
    */

    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    if (!is_detImpZone || string_string)
    {
        double dt = root;
        if (mstate == MotionState::STATIC)
            dt = CollisionSolver3d::getTimeStepSize();
        EdgeToEdgeImpulse(pts,vec,sC,tC,dist,mstate,dt);
    }

	return true;
}

static void EdgeToEdgeImpulse(
        POINT** pts,
        double* nor,
        double a,
        double b,
        double dist,
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

    if (mstate == MotionState::MOVING)
    {
        // Apply one or the other for collision, NOT BOTH.
        // Zero the relative velocity with inelastic impulse.
        if (vn < 0.0)
            EdgeToEdgeInelasticImpulse(mag_vrel,pts,inelastic_impulse,rigid_impulse,wab);
        else if (vn * dt <  overlap_coef * overlap)
            EdgeToEdgeElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }
    else
    {
        // May apply both for repulsion.
        // Zero the normal component of relative velocity with inelastic impulse.
        if (vn < 0.0)
            EdgeToEdgeInelasticImpulse(vn,pts,inelastic_impulse,rigid_impulse,wab);
        if (fabs(vn)*dt < overlap_coef * overlap)
            EdgeToEdgeElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }
    
    double impulse[2];
    double m_impulse[2];
    double f_impulse[2];

    for (int i = 0; i < 2; ++i)
    {
        impulse[i] = inelastic_impulse[i] + elastic_impulse[i];

        if (wab[0] + wab[1] < MACH_EPS || wab[2] + wab[3] < MACH_EPS)
        {
            m_impulse[i] = impulse[i];
            f_impulse[i] = inelastic_impulse[i];
        }
        else
        {
            double wabs_sqr = sqr(wab[0]) + sqr(wab[1])
                              + sqr(wab[2]) + sqr(wab[3]);
        
            m_impulse[i] = 2.0*impulse[i]/wabs_sqr;
            f_impulse[i] = 2.0*inelastic_impulse[i]/wabs_sqr;
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

    //TODO: Do this correctly using impulse method
	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
	    if (isMovableRigidBody(pts[0]))
            SpreadImpactZoneImpulse(pts[0], -1.0*rigid_impulse[0], nor);
        if (isMovableRigidBody(pts[2]))
            SpreadImpactZoneImpulse(pts[2], rigid_impulse[1], nor);
	    return;
	}
	
    double max_friction = 0.5*vt;
    if ((isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])) ||
        (isStaticRigidBody(pts[2]) && isStaticRigidBody(pts[3])))
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
            
            if (t_impulse < 0) continue;
    
            if (mstate == MotionState::STATIC)
            {
                for (int j = 0; j < 3; ++j)
                    sl[i]->collsnImpulse[j] += W[i]*t_impulse*nor[j];
       
                //TODO: only use inelastic component correct?
                //
                //double friction_impulse = M[i];
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
            else
            {
                for (int j = 0; j < 3; ++j)
                    sl[i]->collsnImpulse[j] -= W[i]*t_impulse*v_rel[j]/mag_vrel;
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
	        clean_up(ERROR);
	    }
	}
}

static void EdgeToEdgeInelasticImpulse(
        double v,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double* W)
{
    v = fabs(v);

    if ((isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])) ||
        (isStaticRigidBody(pts[2]) && isStaticRigidBody(pts[3])))
    {
        impulse[0] = v;
        impulse[1] = v;
        rigid_impulse[0] = v;
        rigid_impulse[1] = v;
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1])
            && isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        double m1 = total_mass(pts[0]->hs);
        double m2 = total_mass(pts[2]->hs);
        rigid_impulse[0] = v * m2 / (m1 + m2);
        rigid_impulse[1] = v * m1 / (m1 + m2);
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]))
    {
        impulse[1] = 0.5 * v; 
        rigid_impulse[0] = 0.5 * v;
    }
    else if (isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        impulse[0] = 0.5 * v;
        rigid_impulse[1] = 0.5 * v;
    }
    else
    {
        impulse[0] = 0.5 * v;
        impulse[1] = 0.5 * v;
    }

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
        double dt, double m,
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
            printf("dt*k*overlap/m = %g,  %g*overlap/dt - fabs(vn) = %g\n",
                    dt*k*overlap/m, overlap_coef, overlap_coef*overlap/dt - fabs(vn));
        }

        double I = std::min(dt*k*overlap/m, (overlap_coef*overlap/dt - fabs(vn)));

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
    double dist;
    double nor[3];
	double w[3] = {0.0};
	double x13[3], x23[3], x43[3], x34[3];
    double tri_nor[3] = {0.0};

	Pts2Vec(pts[0],pts[2],x13);
	Pts2Vec(pts[1],pts[2],x23);
	Pts2Vec(pts[3],pts[2],x43);
    scalarMult(-1.0,x43,x34);
	
	double det = Dot3d(x13,x13)*Dot3d(x23,x23)-Dot3d(x13,x23)*Dot3d(x13,x23);
	if (fabs(det) < MACH_EPS)
    {   
        //TODO: Create cgal interface construction library
        //      and remove when this case can be safely omitted.
        printf("\n\tPointToTri() ERROR: fabs(det) < MACH_EPS\n\n");
        clean_up(ERROR);
	}
	else
    {
        //unit normal vector of the plane of the triangle
	    Cross3d(x13,x23,tri_nor);
	    double tri_nor_mag = Mag3d(tri_nor);

        scalarMult(1.0/tri_nor_mag,tri_nor,tri_nor);
        double tri_area = 0.5*tri_nor_mag;

	    //correct the triangle's normal direction to point to same
        //side as the point (not used right now, but may need at some
        //for detecting/correcting interpenetration etc.)
        dist = Dot3d(x34,tri_nor);
        if (dist < 0.0)
        {
            scalarMult(-1.0,tri_nor,tri_nor);
        }
	
        dist = fabs(dist);
        if (dist > tol)
	        return false;
	
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

        double eps = tol/sqrt(tri_area);
        for (int i = 0; i < 3; ++i)
        {
            if (w[i] < -1.0*eps || w[i] > 1.0 + eps)
                return false;
        }
    }

    bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
    if (!is_detImpZone)
    {
        double dt = root;
        if (mstate == MotionState::STATIC)
            dt = CollisionSolver3d::getTimeStepSize();
        PointToTriImpulse(pts,tri_nor,w,dist,mstate,dt);
    }

	return true;
}

static void PointToTriImpulse(
        POINT** pts,
        double* nor,
        double* w,
        double dist,
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

	//apply impulses to the average (linear trajectory) velocity
	for (int i = 0; i < 3; ++i)
	{
	    v_rel[i] += sl[3]->avgVel[i];
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
	
    //Point and Triangle are approaching each other (vn < 0.0):
    //      apply inelastic impulse
    
    //Point and Triangle are seperating from each other (vn > 0.0):
    //      apply elastic impulse
    
    if (mstate == MotionState::MOVING)
    {
        // Apply one or the other for collision, NOT BOTH.
        // Zero the relative velocity with inelastic impulse.
        if (vn < 0.0)
            PointToTriInelasticImpulse(mag_vrel,pts,inelastic_impulse,rigid_impulse,w,&sum_w);
        else if (vn*dt < overlap_coef*overlap)
            PointToTriElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }
    else
    {
        // May apply both for repulsion.
        // Zero the normal component of relative velocity with inelastic impulse.
        if (vn < 0.0)
            PointToTriInelasticImpulse(vn,pts,inelastic_impulse,rigid_impulse,w,&sum_w);
        if (fabs(vn)*dt < overlap_coef*overlap)
            PointToTriElasticImpulse(vn,overlap_coef,overlap,pts,
                    elastic_impulse,rigid_impulse,dt,m,k);
    }

    double impulse[2];
    double m_impulse[2];
    double f_impulse[2];

    for (int i = 0; i < 2; ++i)
    {
        impulse[i] = inelastic_impulse[i] + elastic_impulse[i];

        if (fabs(sum_w) < MACH_EPS)
        {
            m_impulse[i] = impulse[i];
            f_impulse[i] = inelastic_impulse[i];
        }
        else
        {
            m_impulse[i] = 2.0*impulse[i]/(1.0 + Dot3d(w, w));
            f_impulse[i] = 2.0*inelastic_impulse[i]/(1.0 + Dot3d(w, w));
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

    //TODO: Do this correctly using impulse method
	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
	    if (isMovableRigidBody(pts[0]))
            SpreadImpactZoneImpulse(pts[0], -1.0*rigid_impulse[0], nor);
	    if (isMovableRigidBody(pts[3]))
            SpreadImpactZoneImpulse(pts[3], rigid_impulse[1], nor);
	    return;
	}

    double max_friction = 0.5*vt;
    if (isStaticRigidBody(pts[3]) ||
       (isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])
        && isStaticRigidBody(pts[2])))
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
            
            if (t_impulse < 0) continue;
                
            if (mstate == MotionState::STATIC)
            {
                for (int j = 0; j < 3; ++j)
                    sl[i]->collsnImpulse[j] += W[i]*t_impulse*nor[j];

                // Apply friction for static proximity repulsions
                
                //TODO: only use inelastic component correct?
                //
                //double friction_impulse = M[i];
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
            else
            {
                for (int j = 0; j < 3; ++j)
                    sl[i]->collsnImpulse[j] -= W[i]*t_impulse*v_rel[j]/mag_vrel;
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
            clean_up(ERROR);
	    }
	}
}

static void PointToTriInelasticImpulse(
        double v,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double* w,
        double* sum_w)
{
    v = fabs(v);

    if (isStaticRigidBody(pts[3]) ||
       (isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])
        && isStaticRigidBody(pts[2])))
    {
        impulse[0] = v;
        impulse[1] = v;
        rigid_impulse[0] = v;
        rigid_impulse[1] = v;
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]) 
            && isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        double m1 = total_mass(pts[0]->hs);
        double m2 = total_mass(pts[3]->hs);
        rigid_impulse[0] = v * m2 / (m1 + m2);
        rigid_impulse[1] = v * m1 / (m1 + m2);
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1])
            && isMovableRigidBody(pts[2]))
    {
        impulse[1] = 0.5 * v;
        rigid_impulse[0] = 0.5 * v;
    }
    else if (isMovableRigidBody(pts[3]))
    {
        impulse[0] = 0.5 * v;
        rigid_impulse[1] = 0.5 * v;
    }
    else
    {
        impulse[0] = 0.5 * v;
        impulse[1] = 0.5 * v;
    }

    for (int i = 0; i < 3; ++i)
    {
        if (isStaticRigidBody(pts[i]))
            w[i] = 0.0;
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
        double dt, double m,
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
            printf("dt*k*overlap/m = %g,  %g*overlap/dt - fabs(vn) = %g\n",
                    dt*k*overlap/m, overlap_coef, overlap_coef*overlap/dt - fabs(vn));
        }

        double I = std::min(dt*k*overlap/m, (overlap_coef*overlap/dt - fabs(vn)));
        
        impulse[0] += 0.5*I;
        impulse[1] += 0.5*I;
        rigid_impulse[0] += 0.5*I;
        rigid_impulse[1] += 0.5*I;
    }
}

