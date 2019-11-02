#include <armadillo>
#include "collid.h"

//static bool EdgeToEdge(POINT**, double,
  //      MotionState mstate = MotionState::STATIC, double root = 0.0);
static Proximity* EdgeToEdge(POINT**, double);
static bool MovingEdgeToEdge(POINT**,double);

static void EdgeToEdgeImpulse(POINT**,double*,double,double,double,MotionState,double);
static void EdgeToEdgeInelasticImpulse(double,POINT**,double*,double*,double*);
static void EdgeToEdgeElasticImpulse(double,double,POINT**,double*,double*,
                                        double,double,double,double);

//static bool PointToTri(POINT**, double,
  //      MotionState mstate = MotionState::STATIC, double root = 0.0);
static Proximity* PointToTri(POINT**, double);
static bool MovingPointToTri(POINT**,double);

//void PointToTriImpulse(POINT**,double*,double*,double,MotionState,double);
static void PointToTriImpulse(POINT**,double*,double*,double,MotionState,double);
static void PointToTriInelasticImpulse(double,POINT**,double*,double*,double*,double*);
static void PointToTriElasticImpulse(double,double,POINT**,double*,double*,
                                        double,double,double,double);

static bool isCoplanar(POINT**,double,double*);
static void unsort_surface_point(SURFACE *surf);

//functions in CollisionSolver3d
void CollisionSolver3d::assembleFromInterface(
	const INTERFACE* intfc,const double dt)
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

void CollisionSolver3d::updateImpactListVelocity(POINT* head)
{
	STATE* sl = NULL;
	POINT* p = head;

	double m = getPointMass();
    double x_cm[3] = {0.0};
    double v_cm[3] = {0.0};
	double L[3] = {0.0}; //angular momentum
	double I[3][3] = {0.0}; //inertia tensor
	double tmp[3][3];
	int num_pts = 0;

    //compute center of mass position and velocity
	while(p)
    {
		num_pts++;
		sorted(p) = YES;
		sl = (STATE*)left_state(p);
	
        for (int i = 0; i < m_dim; ++i)
        {
		    x_cm[i] += sl->x_old[i]; 
		    v_cm[i] += sl->avgVel[i];
		}

        p = next_pt(p);
    }

    if (debugging("collision"))
	    printf("%d number of points in this zone\n",num_pts);

	for (int i = 0; i < m_dim; ++i)
    {
	    x_cm[i] /= num_pts;
	    v_cm[i] /= num_pts;
	}

	//compute angular momentum
	p = head;
	while(p)
    {
	    double dx[3], dv[3], Li[3];
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    minusVec(sl->avgVel,v_cm,dv); 	
	    Cross3d(dx,dv,Li);
	    scalarMult(m,Li,Li);
	    addVec(Li,L,L);    
	    p = next_pt(p);
	}

	//compute Inertia tensor
	p = head;
	while(p)
    {
	    double dx[3], mag_dx = 0.0;
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    mag_dx = Mag3d(dx);
	   
        for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < 3; ++j)
        {
		    tmp[i][j] = -dx[i]*dx[j];
		    if (i == j)
                tmp[i][j] += mag_dx*mag_dx; 
	
            I[i][j] += tmp[i][j]*m;
	    }

	    p = next_pt(p);
	}

	//compute angular velocity w: I*w = L;
    double mag_w = 0;
	double w[3] = {0.0};

    if (myDet3d(I) > ROUND_EPS)
    {
        for (int i = 0; i < 3; ++i)
        {
            memcpy(tmp,I,9*sizeof(double));
            for (int j = 0; j < 3; j++)
                tmp[j][i] = L[j];
            w[i] = myDet3d(tmp)/myDet3d(I);
        }
    }
    else
    {
        //I is non-invertible, calculate pseudoinverse with SVD
        arma::vec arL(3);
        arma::mat arI(3,3);

        for (int i = 0; i < 3; i++)
        {
            arL(i) = L[i];
            for (int j = 0; j < 3; j++)
                 arI(i,j) = I[i][j];
        }

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

    mag_w = Mag3d(w);
	
	//compute average velocity for each point
	double dt = getTimeStepSize();
	
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
        sl = (STATE*)left_state(p);
        
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
            scalarMult(sin(dt*mag_w)/mag_w,w,tmpV);
            Cross3d(tmpV,xR,wxR);
        }
    
        //TODO: move x_new calculation to above if else block,
        for (int i = 0; i < 3; ++i)
        {
            x_new[i] = x_cm[i] + dt*v_cm[i]
                + xF[i] + cos(dt*mag_w)*xR[i] + wxR[i];
    
            sl = (STATE*)left_state(p);
            sl->avgVel[i] = (x_new[i] - sl->x_old[i])/dt;
   
            if (std::isnan(sl->avgVel[i]))
            { 
                printf("coords[3], vel[3]\n");
                p = head;
                while(p)
                {
                    sl = (STATE*)left_state(p);
                    printf("%f %f %f %f %f %f;\n",
                    sl->x_old[0],sl->x_old[1],sl->x_old[2],
                    sl->avgVel[0],sl->avgVel[1],sl->avgVel[2]);
                    p = next_pt(p);
                }

                printf("num_pts = %d, weight = %d\n",
                num_pts,weight(head));
                printf("nan vel, w = %f, mag_w = %f\n",
                w[i],mag_w);
                printf("L = [%f %f %f]\n",L[0],L[1],L[2]);
                printf("I = [%f %f %f;  %f %f %f; %f %f %f]\n",
                I[0][0],I[0][1],I[0][2],I[1][0],I[1][1],I[1][2],
                I[2][0],I[2][1],I[2][2]);
                printf("xF = %f %f %f, xR = %f %f %f\n",
                xF[0],xF[1],xF[2],xR[0],xR[1],xR[2]);
                clean_up(ERROR);
            }
        }
    
        p = next_pt(p);
    }
}

bool MovingTriToBond(const TRI* tri,const BOND* bd, double h)
{
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
	POINT* pts[4];
	bool status = false;

	/* do not consider bond point that is a tri vertex */
	for (int i = 0; i < 3; ++i)
	{
	    if (Point_of_tri(tri)[i] == bd->start
            || Point_of_tri(tri)[i] == bd->end)
            return false;
	}

	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	/* detect collision of start point of bond w.r.t to tri */
	pts[3] = bd->start;
    if (MovingPointToTri(pts,h))
        status = true;
    
    if (status && is_detImpZone)
        createImpZone(pts,4);
	
    /* detect collision of end point of bond to w.r.t. tri */
	pts[3] = bd->end;
    if (MovingPointToTri(pts,h))
        status = true;

    if (status && is_detImpZone)
        createImpZone(pts,4);
	
    /* detect collision of each of tri edge w.r.t to bond */
	pts[2] = bd->start;
	pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
        if (MovingEdgeToEdge(pts,h))
            status = true;

        if (status && is_detImpZone)
            createImpZone(pts,4);
	}

    return status;
}

bool MovingBondToBond(const BOND* b1, const BOND* b2, double h)
{
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
	POINT* pts[4];
	bool status = false;

	pts[0] = b1->start;
	pts[1] = b1->end;
	pts[2] = b2->start;
	pts[3] = b2->end;

	/* do not consider two bonds that share a common point */
	if (pts[0] == pts[2] || pts[0] == pts[3]
        || pts[1] == pts[2] || pts[1] == pts[3])
    {
        return false;
    }

    /*
    for (int i = 0; i < 4; ++i)
	{
	    for (int j = i + 1; j < 4; ++j)
	    {
            if (pts[i] == pts[j])
                return false;
	    }
	}
    */

    /* detect collision between two bonds */	
    if(MovingEdgeToEdge(pts,h))
        status = true;

    if (status && is_detImpZone)
        createImpZone(pts,4);

    return status;
}

bool MovingTriToTri(const TRI* a,const TRI* b, double h)
{
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
	POINT* pts[4];
	bool status = false;

	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    if (Point_of_tri(a)[i] == Point_of_tri(b)[j])
            return false;
	}

	//detect point to tri collision
	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
    {
	    const TRI* tmp_tri1 = (k == 0) ? a : b;
	    const TRI* tmp_tri2 = (k == 0) ? b : a;

	    for (int j = 0; j < 3; ++j)
            pts[j] = Point_of_tri(tmp_tri1)[j];
        pts[3] = Point_of_tri(tmp_tri2)[i];

	    //Don't consider point against the triangle it belongs to
	    if (pts[3] == pts[0] || pts[3] == pts[1]
                || pts[3] == pts[2])
            continue; 
        
        if(MovingPointToTri(pts,h))
            status = true;

        if (status && is_detImpZone)
            createImpZone(pts,4);
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
		
		    //Don't consider edges with a shared enpoint
            if (pts[0] == pts[2] || pts[0] == pts[3]
            || pts[1] == pts[2] || pts[1] == pts[3])
                continue;
    
            if(MovingEdgeToEdge(pts,h))
                status = true;
                
            if (status && is_detImpZone)
                createImpZone(pts,4);
	    }
    }

	return status;
}

//NOTE: Changes coords of pts involved
static bool MovingPointToTri(POINT* pts[],const double h)
{
	STATE* sl;
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};
    MotionState mstate = MotionState::MOVING;
    bool status = false;

	if (isCoplanar(pts,dt,roots))
    {
	    for (int i = 0; i < 4; ++i)
        {
            if (roots[i] < 0)
                continue;
    
            for (int j = 0; j < 4; ++j)
            {
                sl = (STATE*)left_state(pts[j]);
    		    for (int k = 0; k < 3; ++k)
                    Coords(pts[j])[k] = sl->x_old[k] + roots[i]*sl->avgVel[k];
		    }
    
            if (PointToTri(pts,h,mstate,roots[i]))
            {
                status = true;
                break;
                //return true;
            }
	    }

        for (int j = 0; j < 4; ++j)
        {

            sl = (STATE*)left_state(pts[j]);
            for (int k = 0; k < 3; ++k)
                Coords(pts[j])[k] = sl->x_old[k];
        }
	}

    return status;
    //return false;
}

static bool MovingEdgeToEdge(POINT* pts[],const double h)
{
	STATE* sl;
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};
    MotionState mstate = MotionState::MOVING;
    bool status = false;

	if (isCoplanar(pts,dt,roots))
    {
        for (int i = 0; i < 4; ++i)
        {
            if (roots[i] < 0)
                continue;
                
            for (int j = 0; j < 4; ++j)
            {
                sl = (STATE*)left_state(pts[j]);
                for (int k = 0; k < 3; ++k)
                    Coords(pts[j])[k] = sl->x_old[k] + roots[i]*sl->avgVel[k];
            }

            if (EdgeToEdge(pts,h,mstate,roots[i]))
            {
                status = true;
                break;
            }
        }

        for (int j = 0; j < 4; ++j)
        {
            sl = (STATE*)left_state(pts[j]);
            for (int k = 0; k < 3; ++k)
                Coords(pts[j])[k] = sl->x_old[k];
        }
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
	//(x1+tv1)x(x2+tv2)*(x3+tv3) = 0
	//transform to at^3+bt^2+ct+d = 0
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
	if (fabs(a) > MACH_EPS){
	    b /= a; c /= a; d /= a;
	    a = b; b = c; c = d;
	    double Q, R, theta;
	    double Q3, R2;
	    Q = (a*a-3*b)/9;
	    R = (2*a*a*a-9*a*b+27*c)/54;
	    Q3 = Q*Q*Q;
	    R2 = R*R;
	    if (R2 < Q3){
	        double Qsqrt = sqrt(Q);
		theta = acos(R/sqrt(Q3));
		roots[0] = -2*Qsqrt*cos(theta/3)-a/3;
		roots[1] = -2*Qsqrt*cos((theta+2*M_PI)/3)-a/3;
		roots[2] = -2*Qsqrt*cos((theta-2*M_PI)/3)-a/3;	
	    }
	    else{
		double A, B;
		double sgn = (R > 0) ? 1.0 : -1.0;
		A = -sgn*pow(fabs(R)+sqrt(R2-Q3),1.0/3.0);
		B = (fabs(A) < ROUND_EPS) ? 0.0 : Q/A;
		roots[0] = (A+B)-a/3.0;
		if (fabs(A-B) < ROUND_EPS)
		    roots[1] = roots[2] = -0.5*(A+B)-a/3.0; //multiple roots
	    }
	}
	else{
		a = b; b = c; c = d;
	   	double delta = b*b-4.0*a*c;
	   	if (fabs(a) > ROUND_EPS && delta > 0){
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
	for (int i = 0; i < 3; ++i){
	        roots[i] = roots[i]-MACH_EPS;
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

bool TriToBond(const TRI* tri,const BOND* bd, double h)
{
	/* do not consider bond point that is a tri vertex */
	for (int i = 0; i < 3; ++i)
	{
	    if (Point_of_tri(tri)[i] == bd->start ||
            Point_of_tri(tri)[i] == bd->end)
            return false;
	}

	POINT* pts[4];
    std::vector<Proximity*> proximities;
	
    for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	//detect proximity of triangle to bond end points.
    for (int i = 0; i < 2; ++i)
    {
        pts[3] = (i == 0) ? bd->start : bd->end;
        
        Proximity* proximity = PointToTri(pts,h);
        if (proximity)
            proximities.push_back(proximity);
    }

	//detect proximity of each triangle edge with bond edge
	pts[2] = bd->start;
    pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];

        Proximity* proximity = EdgeToEdge(pts,h);
        if (proximity)
            proximities.push_back(proximity);
	}

    Proximity* closest;
    double min_dist = HUGE;

    std::vector<Proximity*>::iterator it;
    for (it = proximities.begin(); it < proximities.end(); ++it)
    {
        if ((*it)->dist < min_dist)
        {
            min_dist = (*it)->dist;
            closest = *it;
        }
    }

    if (closest)
    {
        closest->applyImpulse();
        for (it = proximities.begin(); it < proximities.end(); ++it)
        {
            delete *it;
        }

        return true;
    }

    return false;
}

//NOTE: this function restores coords of points to x_old
bool BondToBond(const BOND* b1, const BOND* b2, double h)
{
	POINT* pts[4];

	pts[0] = b1->start; pts[1] = b1->end;
	pts[2] = b2->start; pts[3] = b2->end;

	//Don't consider two bonds that share a common point
	if (pts[0] == pts[2] || pts[0] == pts[3]
        || pts[1] == pts[2] || pts[1] == pts[3])
    {
        return false;
    }
    
    Proximity* proximity = EdgeToEdge(pts,h);
    if (proximity)
    {
        proximity->applyImpulse();
        delete proximity;
        return true;
    }

    return false;
}

bool TriToTri(const TRI* tri1, const TRI* tri2, double h)
{
	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    if (Point_of_tri(tri1)[i] == Point_of_tri(tri2)[j])
            return false;
	}

	POINT* pts[4];
    std::vector<Proximity*> proximities;

	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
	{
	    const TRI* tmp_tri1 = (k == 0) ? tri1 : tri2;
	    const TRI* tmp_tri2 = (k == 0) ? tri2 : tri1;

	    pts[3] = Point_of_tri(tmp_tri1)[i];
        for (int j = 0; j < 3; ++j)
            pts[j] = Point_of_tri(tmp_tri2)[j];
	
	    //Don't check a point against the triangle that contains it
        if (pts[0] == pts[3] ||
            pts[1] == pts[3] || pts[2] == pts[3])
            continue;

        Proximity* proximity = PointToTri(pts,h);
        if (proximity)
            proximities.push_back(proximity);
	}

	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri1)[i];
	    pts[1] = Point_of_tri(tri1)[(i+1)%3];
	
        for (int j = 0; j < 3; ++j)
        {
		    pts[2] = Point_of_tri(tri2)[j];
		    pts[3] = Point_of_tri(tri2)[(j+1)%3];
		
            //Don't check edges that share an endpoint
            if (pts[0] == pts[2] || pts[0] == pts[3] ||
                pts[1] == pts[2] || pts[1] == pts[3])
                continue;

            Proximity* proximity = EdgeToEdge(pts,h);
            if (proximity)
                proximities.push_back(proximity);
	    }
    }

    Proximity* closest;
    double min_dist = HUGE;

    std::vector<Proximity*>::iterator it;
    for (it = proximities.begin(); it < proximities.end(); ++it)
    {
        if ((*it)->dist < min_dist)
        {
            min_dist = (*it)->dist;
            closest = *it;
        }
    }

    if (closest)
    {
        closest->applyImpulse();
        for (it = proximities.begin(); it < proximities.end(); ++it)
        {
            delete *it;
        }

        return true;
    }

    return false;
}

static void PointToLine(POINT* pts[],double &a)
{
/*
*	x1 -----projP---- x2
*		  |
*		  |  dist
*		  *x3
*/
	double x12[3], x13[3];
	Pts2Vec(pts[0],pts[1],x12);
    Pts2Vec(pts[0],pts[2],x13);
    a = Dot3d(x13,x12)/Dot3d(x12,x12);
}

/*
static bool EdgeToEdge(POINT** pts, double h, double root)
{
//     x1	  x3
//	    /	   \
//     /	    \
// x2 /		     \ x4
// solve equation
// x21*x21*a - x21*x43*b = x21*x31
// -x21*x43*a + x43*x43*b = -x43*x31
//
	double x21[3], x43[3], x31[3];
	double a, b;
	double tmp[3];
	double v1[3],v2[3];
	double nor[3], nor_mag, dist;

	Pts2Vec(pts[1],pts[0],x21);    
	Pts2Vec(pts[3],pts[2],x43);
	Pts2Vec(pts[2],pts[0],x31);
	Cross3d(x21,x43,tmp);
	if (Mag3d(tmp) < ROUND_EPS)
	{
        //TODO: is this working?
	     return false; //ignore the case where two edges are parallel??
	    
        //degenerate cases to parallel line segments
        if (Mag3d(x21) > ROUND_EPS || Mag3d(x43) > ROUND_EPS){
	    	POINT* plist[3];
		double tmp_min_dist = HUGE;
		for (int i = 0; i < 4; ++i){
		    // p0-> p2--p3
		    // p1-> p2--p3
		    // p2-> p0--p1
		    // p3-> p0--p1
		    plist[0] = (i%2 == 0) ? pts[(i+2)%4] : pts[(i+1)%4];
		    plist[1] = (i%2 == 0) ? pts[(i+3)%4] : pts[(i+2)%4];
		    plist[2] = pts[i];
		    double tmp_vec[3], tmp_nor[3], tmp_dist, tmp_a;
		    minusVec(Coords(plist[1]),Coords(plist[0]),tmp_vec);
		    if (Mag3d(tmp_vec) < ROUND_EPS) continue;

		    PointToLine(plist,tmp_a);
		    tmp_a = std::max(std::min(tmp_a,1.0),0.0);
		    scalarMult(tmp_a,tmp_vec,tmp_vec);
		    addVec(Coords(plist[0]),tmp_vec,tmp_vec);
		    if (i/2 == 0)
		        minusVec(Coords(plist[2]),tmp_vec,tmp_nor);
		    else 
		        minusVec(tmp_vec,Coords(plist[2]),tmp_nor);
		    tmp_dist = distance_between_positions(
				tmp_vec,Coords(plist[2]),3);
		    if (tmp_dist < tmp_min_dist){
			memcpy((void*)nor,(void*)tmp_nor,3*sizeof(double));
			dist = tmp_min_dist = tmp_dist;
			if	(i == 0){a = 0.0; b = tmp_a;}
			else if (i == 1){a = 1.0; b = tmp_a;}
			else if (i == 2){a = tmp_a; b = 0.0;}
			else if (i == 3){a = tmp_a; b = 1.0;}
		    }
		}		

		if (dist < ROUND_EPS){
		    memcpy((void*)nor,
			   (Mag3d(x21) < ROUND_EPS) ? (void*)x43 : (void*)x21,
			    3*sizeof(double));
		}

	    }
	    else{
	   	//both x21 and x43 degenerate to points
		minusVec(Coords(pts[0]),Coords(pts[2]),nor);
		dist = Mag3d(nor);
		if (dist < ROUND_EPS){
		    nor[0] = nor[1] = nor[2] = 1.0;
		}
		a = 0.0; b = 0.0;
	    }
	}
	else
	{
	    a = (Dot3d(x43,x43)*Dot3d(x21,x31)-Dot3d(x21,x43)*Dot3d(x43,x31))/
	        (Dot3d(x21,x21)*Dot3d(x43,x43)-Dot3d(x21,x43)*Dot3d(x21,x43)); 
	    b = (Dot3d(x21,x43)*Dot3d(x21,x31)-Dot3d(x21,x21)*Dot3d(x43,x31))/
                (Dot3d(x21,x21)*Dot3d(x43,x43)-Dot3d(x21,x43)*Dot3d(x21,x43));
	    a = std::max(std::min(a,1.0),0.0);	
	    b = std::max(std::min(b,1.0),0.0);	
	    scalarMult(a,x21,v1);
	    scalarMult(b,x43,v2);
	    addVec(Coords(pts[0]),v1,v1);
	    addVec(Coords(pts[2]),v2,v2);
	    minusVec(v2,v1,nor);
	    nor_mag = Mag3d(nor);
	    //if (nor_mag < ROUND_EPS)
	    if (nor_mag < 1000 * MACH_EPS)
	    {
		//v1 == v2;
                //two edges intersect with each other
                //normal direction is calculated with old position
		STATE* sl[4];
		for (int i = 0; i < 4; ++i)
		    sl[i] = (STATE*)left_state(pts[i]);
		for (int j = 0; j < 3; ++j)
		{
		    nor[j]  = (1.0-b) * sl[2]->x_old[j] + b * sl[3]->x_old[j];
		    nor[j] -= (1.0-a) * sl[0]->x_old[j] + a * sl[1]->x_old[j];
		}
	    }
	    dist = distBetweenCoords(v1,v2);
	}
	if (dist > h) 
	    return false;
	nor_mag = Mag3d(nor);
	if (nor_mag < MACH_EPS)
	{
            printf("Normal vector is degenerate:\t");
            printf("a = %f, b = %f\n",a,b);
            printf("x_old:\n");
            for (int i = 0; i < 4; ++i){
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%f %f %f\n",sl1->x_old[0],sl1->x_old[1],sl1->x_old[2]);
            }
            printf("x_new:\n");
            for (int i = 0; i < 4; ++i){
                printf("%f %f %f\n",Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
            }
            printf("avgVel:\n");
            for (int i = 0; i < 4; ++i){
                STATE* sl1 = (STATE*)left_state(pts[i]);
                printf("%f %f %f\n",sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
            }
            clean_up(ERROR);
	}
	else
        for (int i = 0; i < 3; ++i)
            nor[i] /= nor_mag;
    
    EdgeToEdgeImpulse(pts, nor, a, b, dist, root);
	return true;
}
*/

//For details of this implementation see:
//http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
//
//Note that mstate has default value of MotionState::STATIC,
//and root has default value of 0.0

/*
static bool EdgeToEdge(
        POINT** pts,
        double h,
        MotionState mstate,
        double root)
{
//	   x1	x3
//	    /	 \
//     /	  \
// x2 /		   \ x4
//
	double x12[3], x34[3], x13[3];
	Pts2Vec(pts[0],pts[1],x12);    
	Pts2Vec(pts[2],pts[3],x34);
	Pts2Vec(pts[0],pts[2],x13);

    //Matrix entries
    double a = Dot3d(x12,x12);
    double b = -Dot3d(x12,x34);
    double c = Dot3d(x34,x34);

    //RHS
    double d = Dot3d(x12,x13);
    double e = -Dot3d(x34,x13);
	
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


    double dist, vec[3];
	Cross3d(x12,x34,vec);
	//if (Mag3d(vec) < ROUND_EPS)
	//{
        //return false;
    //}


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

    //Compute the solution and obtain the closest pair of points.
    if (sN == sD)
        sC = 1.0;
    else
        sC = fabs(sN) < ROUND_EPS ? 0.0 : sN/sD;

    if (tN == tD)
        tC = 1.0;
    else
        tC = fabs(tN) < ROUND_EPS ? 0.0 : tN/tD;
    
    
    //
    //if (std::isinf(sC) || std::isinf(tC))
    //{
    //    LOC();
    //    printf("sC or tC is inf\n");
    //    clean_up(ERROR);
    //}
    //
    
    //"normal vector" always points from
    //the x12 edge to the x34 edge
    double nor[3];
	for (int j = 0; j < 3; ++j)
	{
	    nor[j]  = (1.0 - tC)*Coords(pts[2])[j] + tC*Coords(pts[3])[j];
	    nor[j] -= (1.0 - sC)*Coords(pts[0])[j] + sC*Coords(pts[1])[j];
    }
    
    //TODO: Make sure nor is pointing in the correct direction.

    //
    //scalarMult(tC,x34,x34);
    //scalarMult(sC,x12,x12);
    //minusVec(x34,x12,nor);
    //minusVec(nor,x31,nor);
    //

    dist = Mag3d(nor);
    if (dist > h)
        return false;

    if (dist > ROUND_EPS)
        scalarMult(1.0/dist,nor,nor);
    else
    {
        std::cout << "dist (nor_mag) < ROUND_EPS" << std::endl;
        std::cout << "sC = " << sC << " tC = " << tC << std::endl;
        printPointList(pts,4);
	    double p12[3]; double p34[3];
        for (int i = 0; i < 3; ++i)
        {
            p34[i] = (1.0 - tC)*Coords(pts[2])[i] + tC*Coords(pts[3])[i];
	        p12[i] = (1.0 - sC)*Coords(pts[0])[i] + sC*Coords(pts[1])[i];
        }
    
        printf("nor \t p34 \t p12\n");
        for (int i = 0; i < 3; ++i)
        {
            printf("%g \t %g \t %g\n",nor[i],p34[i],p12[i]);
        }
        clean_up(ERROR);
    }
    
    EdgeToEdgeImpulse(pts,nor,sC,tC,dist,mstate,root);
	return true;
}
*/

static Proximity* EdgeToEdge(POINT** pts, double h)
{
//	   x1	x3
//	    /	 \
//     /	  \
// x2 /		   \ x4
//
	double x12[3], x34[3], x13[3];
	Pts2Vec(pts[0],pts[1],x12);    
	Pts2Vec(pts[2],pts[3],x34);
	Pts2Vec(pts[0],pts[2],x13);

    //Matrix entries
    double a = Dot3d(x12,x12);
    double b = -Dot3d(x12,x34);
    double c = Dot3d(x34,x34);

    //RHS
    double d = Dot3d(x12,x13);
    double e = -Dot3d(x34,x13);
	
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


    double dist, vec[3];
	Cross3d(x12,x34,vec);
	//if (Mag3d(vec) < ROUND_EPS)
	//{
        //return false;
    //}


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

    //Compute the solution and obtain the closest pair of points.
    if (sN == sD)
        sC = 1.0;
    else
        sC = fabs(sN) < ROUND_EPS ? 0.0 : sN/sD;

    if (tN == tD)
        tC = 1.0;
    else
        tC = fabs(tN) < ROUND_EPS ? 0.0 : tN/tD;
    
    
    /*
    if (std::isinf(sC) || std::isinf(tC))
    {
        LOC();
        printf("sC or tC is inf\n");
        clean_up(ERROR);
    }
    */
    
    
    //TODO: Make sure nor is pointing in the correct direction.
    
    //"normal vector" always points from
    //the x12 edge to the x34 edge
    
    double nor[3];
	for (int j = 0; j < 3; ++j)
	{
	    nor[j]  = (1.0 - tC)*Coords(pts[2])[j] + tC*Coords(pts[3])[j];
	    nor[j] -= (1.0 - sC)*Coords(pts[0])[j] + sC*Coords(pts[1])[j];
    }

    dist = Mag3d(nor);
    if (dist > h)
        return {};

    if (dist > ROUND_EPS)
        scalarMult(1.0/dist,nor,nor);
    else
    {
        std::cout << "dist (nor_mag) < ROUND_EPS" << std::endl;
        std::cout << "sC = " << sC << " tC = " << tC << std::endl;
        printPointList(pts,4);
	    double p12[3]; double p34[3];
        for (int i = 0; i < 3; ++i)
        {
            p34[i] = (1.0 - tC)*Coords(pts[2])[i] + tC*Coords(pts[3])[i];
	        p12[i] = (1.0 - sC)*Coords(pts[0])[i] + sC*Coords(pts[1])[i];
        }
    
        printf("nor \t p34 \t p12\n");
        for (int i = 0; i < 3; ++i)
        {
            printf("%g \t %g \t %g\n",nor[i],p34[i],p12[i]);
        }
        clean_up(ERROR);
    }

    Proximity* proximity = new EdgeEdgeProximity(pts,sC,tC,w,dist);
    return proximity;
}

//Note that mstate has default value of MotionState::STATIC,
//and root has default value 0.0

void EdgeToEdgeProximityImpulse(
        POINT** pts,
        double* nor,
        double a,
        double b,
        double dist)
{
    MotionState mstate = MotionState::STATIC;
    EdgeToEdgeImpulse(pts,nor,a,b,dist,mstate,-1.0);
}

static void EdgeToEdgeImpulse(
        POINT** pts,
        double* nor,
        double a,
        double b,
        double dist,
        MotionState mstate,
        double root)
{
	if (debugging("collision"))
	    CollisionSolver3d::edg_to_edg++;

	double k      = CollisionSolver3d::getSpringConstant();
	double m      = CollisionSolver3d::getPointMass();
	double dt     = CollisionSolver3d::getTimeStepSize();
	double mu     = CollisionSolver3d::getFrictionConstant(); 
	double h      = CollisionSolver3d::getFabricThickness();
	double cr     = CollisionSolver3d::getRestitutionCoef();
	
    dist   = h - dist;

	STATE *sl[4];
	for (int i = 0; i < 4; ++i)
	    sl[i] = (STATE*)left_state(pts[i]);

	double v_rel[3];
	for (int i = 0; i < 3; ++i)
	{
	    v_rel[i]  = (1.0 - b)*sl[2]->avgVel[i] + b*sl[3]->avgVel[i];
	    v_rel[i] -= (1.0 - a)*sl[0]->avgVel[i] + a*sl[1]->avgVel[i];
	}
	
    double vn = Dot3d(v_rel,nor);
    double vt = sqrt(Dot3d(v_rel,v_rel) - sqr(vn));
	
    /*
    if (Dot3d(v_rel, v_rel) > sqr(vn))
	    vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	else
	    vt = 0.0;
    */

    //Edges are approaching each other (vn < 0.0):
    //      apply inelastic impulse
    //Edges are seperating from each other (vn > 0.0):
    //      apply elastic impulse
    
	double impulse = 0.0;
	double friction_impulse = 0.0;
    double rigid_impulse[2] = {0.0};
    
    double wab[4] = {1.0 - a, a, 1.0 - b, b};

    if (mstate == MotionState::MOVING)
    {
        //Apply one or the other for collision, NOT BOTH
        dt = root;
        if (vn < 0.0)
            EdgeToEdgeInelasticImpulse(vn,pts,&impulse,rigid_impulse,wab);
        else if (vn * dt < 0.1 * dist)
            EdgeToEdgeElasticImpulse(vn,dist,pts,&impulse,rigid_impulse,dt,m,k,cr);
    }
    else
    {
        //Apply both if needed for proximity
        if (vn < 0.0)
            EdgeToEdgeInelasticImpulse(vn,pts,&impulse,rigid_impulse,wab);
        if (vn * dt < 0.1 * dist)
            EdgeToEdgeElasticImpulse(vn,dist,pts,&impulse,rigid_impulse,dt,m,k,cr);
        if (vt > ROUND_EPS)
        {
            double delta_vt = 0.5*vt;
            if (fabs(mu*impulse) < vt)
                delta_vt = fabs(mu*impulse);

            //TODO: check sign of update
            friction_impulse = delta_vt;
            //sl[i]->friction[j] += W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
            //sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
        }
    }

    double m_impulse; double f_impulse;
	if (wab[0] + wab[1] < MACH_EPS || wab[2] + wab[3] < MACH_EPS)
    {
	    m_impulse = impulse;
	    f_impulse = friction_impulse;
    }
    else
    {
        double wabs_sqr = sqr(wab[0]) + sqr(wab[1])
                            + sqr(wab[2]) + sqr(wab[3]);
        
        m_impulse = 2.0*impulse/wabs_sqr;
        f_impulse = 2.0*friction_impulse/wabs_sqr;
    }

    //uncomment the following for debugging
    if (debugging("CollisionImpulse"))
    if (fabs(m_impulse) > 0.0){
        printf("real EdgeToEdge collision\n");
        printf("vt = %f, vn = %f, dist = %f\n",vt,vn,dist);
        printf("v_rel = %f %f %f\n",v_rel[0],v_rel[1],v_rel[2]);
        printf("nor = %f %f %f\n",nor[0],nor[1],nor[2]);
        printf("m_impulse = %f, impulse = %f, a = %f, b = %f\n",
            m_impulse,impulse,a,b);
        printf("root = %e,h = %e, dt = %e\n",root,h,dt);
        printf("x_old:\n");
        for (int i = 0; i < 4; ++i){
            STATE* sl1 = (STATE*)left_state(pts[i]);
            printf("%f %f %f\n",sl1->x_old[0],sl1->x_old[1],sl1->x_old[2]);
        }
        printf("x_new:\n");
        for (int i = 0; i < 4; ++i){
            printf("%f %f %f\n",Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
        }
        printf("avgVel:\n");
        for (int i = 0; i < 4; ++i){
            STATE* sl1 = (STATE*)left_state(pts[i]);
            printf("%f %f %f\n",sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
        }
        printf("\n");
    }
    ////////////////////////////////////////////////

    //TODO: This doesn't look right.
    //      Also should be its own function.
	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
	    if (isMovableRigidBody(pts[0]))
            SpreadImpactZoneImpulse(pts[0], rigid_impulse[0], nor);
        if (isMovableRigidBody(pts[2]))
            SpreadImpactZoneImpulse(pts[2], -1.0 * rigid_impulse[1], nor);
	    return;
	}
	

    std::vector<double> W = {wab[0],wab[1],-wab[2],-wab[3]};
    std::vector<double> R = {rigid_impulse[0],rigid_impulse[0],
                             rigid_impulse[1],rigid_impulse[1]};

    //TODO: should probably be in own function together with
    //      the elastic and inelastic impulse functions above
    for (int i = 0;  i < 4; ++i)
    {
        if (!isStaticRigidBody(pts[i]))
        {
            //TODO: distinguish between collsn_impulse and proximity_impulse
            sl[i]->collsn_num++;

            double t_impulse = m_impulse;
            double t_friction_impulse = f_impulse;
            if (isMovableRigidBody(pts[i]))
            {
                t_impulse = R[i];
            }

            for (int j = 0; j < 3; ++j)
            {
                sl[i]->collsnImpulse[j] += W[i]*t_impulse*nor[j];
        
                if (vt > ROUND_EPS)
                    sl[i]->friction[j] += W[i]*t_friction_impulse*(v_rel[j] - vn*nor[j])/vt;
       
                /*
                //Friction only applied for proximities, not collisions.
                if (mstate == MotionState::MOVING)
                    continue;

                if (fabs(vt) > ROUND_EPS)
                {
                    double delta_vt = vt;
                    if (fabs(lambda*t_impulse) < vt)
                        delta_vt = fabs(lambda*t_impulse);
   
                    //TODO: check sign of update
                    //sl[i]->friction[j] += W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
                    sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
                   
                    //double delta_vt = vt;
                    //if (fabs(lambda*t_impulse) < vt)
                      //  delta_vt = lambda*t_impulse;

                    //sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;

                    //double delta_vt = std::max(-fabs(lambda*W[i]*t_impulse/vt), -1.0);
                    //sl[i]->friction[j] += delta_vt*(v_rel[j] - vn*nor[j]);
                    //double delta_vt = std::min(1.0,fabs(lambda*delta_vn/vt));
                }
                */
            }
        }
        else
        {
            for (int j = 0; j < 3; ++j)
            {
                sl[i]->collsnImpulse[j] = 0.0;
                sl[i]->friction[j] = 0.0;
            }
        }
    }

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 3; ++j)
    {
	    if (std::isnan(sl[i]->collsnImpulse[j]) ||
            std::isinf(sl[i]->collsnImpulse[j]))
        {
		    printf("EdgeToEdge: sl[%d]->impl[%d] = nan\n",i,j);
		    printf("a b = %f %f, nor = [%f %f %f], dist = %f\n",
                    a,b,nor[0],nor[1],nor[2],dist);
	        clean_up(ERROR);
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
    if ((isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])) ||
        (isStaticRigidBody(pts[2]) && isStaticRigidBody(pts[3])))
    {
        *impulse += vn;
        rigid_impulse[0] += vn;
        rigid_impulse[1] += vn;
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1])
            && isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        double m1 = total_mass(pts[0]->hs);
        double m2 = total_mass(pts[2]->hs);
        rigid_impulse[0] += vn * m2 / (m1 + m2);
        rigid_impulse[1] += vn * m1 / (m1 + m2);
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]))
    {
        *impulse += 0.5 * vn; 
        rigid_impulse[0] += 0.5 * vn;
    }
    else if (isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        *impulse += 0.5 * vn;
        rigid_impulse[1] += 0.5 * vn;
    }
    else
    {
        *impulse += 0.5 * vn;
    }

    if (isStaticRigidBody(pts[0])) W[0] = 0.0;
    if (isStaticRigidBody(pts[1])) W[1] = 0.0;
    if (isStaticRigidBody(pts[2])) W[2] = 0.0;
    if (isStaticRigidBody(pts[3])) W[3] = 0.0;
}

static void EdgeToEdgeElasticImpulse(
        double vn,
        double dist,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double dt, double m,
        double k, double cor)
{
    if (isRigidBody(pts[0]) && isRigidBody(pts[1])
        && isRigidBody(pts[2]) && isRigidBody(pts[3]))
    {
        rigid_impulse[0] *= 1.0 + cor;
        rigid_impulse[1] *= 1.0 + cor;
    }
    else
    {
        double tmp = -1.0*std::min(dt*k*dist/m, (0.1*dist/dt - vn));
        *impulse += tmp;
        rigid_impulse[0] += tmp;
        rigid_impulse[1] += tmp;
    }
}

//Note that mstate has default value of MotionState::STATIC,
//and root has default value of 0.0
/*
static bool PointToTri(
        POINT** pts,
        double h,
        MotionState mstate,
        double root)
{
//	x1
 *  	/\     x4 *
 *     /  \
 * x2 /____\ x3
 *
 * solve equation
 * x13*x13*w1 + x13*x23*w2 = x13*x43
 * x13*x23*w1 + x23*x23*w2 = x23*x43
 //
    double dist;
    //double nor[3];
	double w[3] = {0.0};
	double x13[3], x23[3], x43[3], x34[3];
    double tri_nor[3] = {0.0};

	Pts2Vec(pts[0],pts[2],x13);
	Pts2Vec(pts[1],pts[2],x23);
	Pts2Vec(pts[3],pts[2],x43);
    scalarMult(-1.0,x43,x34);
	
	double det = Dot3d(x13,x13)*Dot3d(x23,x23)-Dot3d(x13,x23)*Dot3d(x13,x23);
	if (fabs(det) < 1000 * MACH_EPS)
    {   
	    // ignore cases where tri reduces to a line or point
         
        //TODO: let's not ignore it
        return false;

        //
	    //consider the case when det = 0
	    //x13 and x23 are collinear

	    POINT* tmp_pts[3]; 
	    double max_len = 0;
	    //find longest side
	    for (int i = 0; i < 3; ++i){
		double tmp_dist = distance_between_positions(Coords(pts[i]),
				Coords(pts[(i+1)%3]),3);
		if (tmp_dist > max_len){
		    tmp_pts[0] = pts[i];
		    tmp_pts[1] = pts[(i+1)%3];
		    max_len = tmp_dist;
		}
	    }
	    if (max_len < ROUND_EPS){
		//triangle degenerate to a point
		minusVec(Coords(pts[0]),Coords(pts[3]),nor);
		dist = Mag3d(nor);
		if (dist < ROUND_EPS)
		for (int i = 0; i < 3; ++i){
		    w[i] = 1.0/3.0;
		    nor[i] = 1.0; 
		}
	    }
	    else{
		//triangle degnerate to a line between tmp_pts
		tmp_pts[2] = pts[3];
		double a;
		PointToLine(tmp_pts,a);
		for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 2; ++j){
		    if (tmp_pts[j] == pts[i]){
			w[i] = (j == 0) ? a : 1.0-a;
		    }
		}
		double v[3];
		minusVec(Coords(tmp_pts[1]),Coords(tmp_pts[0]),v);
		scalarMult(a,v,v);
		addVec(Coords(tmp_pts[0]),v,v);
		minusVec(Coords(pts[3]),v,nor);
		dist = Mag3d(nor);
		//define nor vec if dist == 0
		if (Mag3d(nor) < ROUND_EPS)
		    nor[0] = nor[1] = nor[2] = 1.0;
	    }
        //
	}
	else
    {
	    //det != 0
	    //x13 and x23 are non-collinear

        //unit normal vector of the plane of the triangle
	    Cross3d(x13,x23,tri_nor);
	    double tri_nor_mag = Mag3d(tri_nor);

        scalarMult(1.0/tri_nor_mag,tri_nor,tri_nor);
        double tri_area = 0.5*tri_nor_mag;

        //
	    //get the old direction
	    //double x34_old[3];
	    //STATE* tmp_sl[2];
	    //tmp_sl[0] = (STATE*)left_state(pts[2]);
	    //tmp_sl[1] = (STATE*)left_state(pts[3]);
	    //minusVec(tmp_sl[1]->x_old,tmp_sl[0]->x_old,x34_old);
        //

	    //correct the triangle's normal direction to point to same
        //side as the point (not used right now, but may need at some
        //for detecting/correcting interpenetration etc.)

        //dist = Dot3d(x34_old,tri_nor);
        dist = Dot3d(x34,tri_nor);
        if (dist < 0.0)
        {
            scalarMult(-1.0,tri_nor,tri_nor);
        }
	
        dist = fabs(dist);
        if (dist > h)
	        return false;
	
	    w[0] = (Dot3d(x23,x23)*Dot3d(x13,x43)-Dot3d(x13,x23)*Dot3d(x23,x43))/det;
	    w[1] = (Dot3d(x13,x13)*Dot3d(x23,x43)-Dot3d(x13,x23)*Dot3d(x13,x43))/det;
	    w[2] = 1.0 - w[0] - w[1];
	    
        //
        //double c_len = 0.0;	
        //for (int i = 0; i < 3; ++i)
        //{
        //    double tmp_dist = distance_between_positions(Coords(pts[i]),
        //                            Coords(pts[(i+1)%3]),3);
        //    if (tmp_dist > c_len)
        //        c_len = tmp_dist;
        //}
        //

        //double eps = h/c_len;
        double eps = h/sqrt(tri_area);
        for (int i = 0; i < 3; ++i)
        {
            if (w[i] < -1.0*eps || w[i] > 1.0 + eps)
                return false;
        }

        //
        //compute the "normal vector" pointing from the
        //projected point of the triangle's plane to the point
        //for (int i = 0; i < 3; ++i)
        //{
          //  nor[i] = x34[i] + w[0]*x13[i] + w[1]*x23[i];
        //}
        //
        
        STATE* tmp_sl = (STATE*)left_state(pts[3]);
        if (fabs(w[0]) < ROUND_EPS || fabs(w[1]) < ROUND_EPS
                || fabs(w[2]) < ROUND_EPS)
        {
        
        for (int j = 0; j < 3; ++j)
            nor[j] = tmp_sl->x_old[j];

        for (int i = 0; i < 3; ++i)
        {
            tmp_sl = (STATE*)left_state(pts[i]);
            for (int j = 0; j < 3; ++j)
                nor[j] -= w[i] * tmp_sl->x_old[j];
        }
            
        double nor_mag = Mag3d(nor);
        if (nor_mag > ROUND_EPS)
        {
            scalarMult(1.0/nor_mag,nor,nor);
        }
        else
        {
            std::cout << "nor_mag < ROUND_EPS" << std::endl;
            printPointList(pts,4);
            clean_up(ERROR);
        }
    }

	PointToTriImpulse(pts,tri_nor,w,dist,mstate,root);
	return true;
}
*/

static Proximity* PointToTri(POINT** pts, double h)
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
    
	double x13[3], x23[3], x43[3], x34[3];
	
    Pts2Vec(pts[0],pts[2],x13);
	Pts2Vec(pts[1],pts[2],x23);
	Pts2Vec(pts[3],pts[2],x43);
    scalarMult(-1.0,x43,x34);
	
	//double det = Dot3d(x13,x13)*Dot3d(x23,x23)-Dot3d(x13,x23)*Dot3d(x13,x23);

    double tri_nor[3] = {0.0};
    Cross3d(x13,x23,tri_nor);
    double tri_nor_mag = Mag3d(tri_nor);

    scalarMult(1.0/tri_nor_mag,tri_nor,tri_nor);
    double tri_area = 0.5*tri_nor_mag;

    double dist = Dot3d(x34,tri_nor);
    if (dist < 0.0)
    {
        dist = fabs(dist);
        scalarMult(-1.0,tri_nor,tri_nor);
    }

    if (dist > h)
        return {};

	double w[3];
    w[0] = (Dot3d(x23,x23)*Dot3d(x13,x43)-Dot3d(x13,x23)*Dot3d(x23,x43))/det;
    w[1] = (Dot3d(x13,x13)*Dot3d(x23,x43)-Dot3d(x13,x23)*Dot3d(x13,x43))/det;
    w[2] = 1.0 - w[0] - w[1];
    
    /*
    double c_len = 0.0;	
    for (int i = 0; i < 3; ++i)
    {
        double tmp_dist = distance_between_positions(Coords(pts[i]),
                                Coords(pts[(i+1)%3]),3);
        if (tmp_dist > c_len)
            c_len = tmp_dist;
    }
    */

    //double eps = h/c_len;
    double eps = h/sqrt(tri_area);
    for (int i = 0; i < 3; ++i)
    {
        if (w[i] < -1.0*eps || w[i] > 1.0 + eps)
            return {};
            //return false;
    }

    /*
    double nor_mag = Mag3d(nor);
    if (nor_mag > ROUND_EPS)
    {
        scalarMult(1.0/nor_mag,nor,nor);
    }
    else
    {
        std::cout << "nor_mag < ROUND_EPS" << std::endl;
        printPointList(pts,4);
        clean_up(ERROR);
    }
    */

    Proximity* proximity = new PointToTriProximity(pts,tri_nor,w,dist);
    return proximity;
}

void PointToTriProximityImpulse(
        POINT** pts,
        double* nor,
        double* w,
        double dist)
{
    MotionState mstate = MotionState::STATIC;
    PointToTriImpulse(pts,nor,w,dist,mstate,-1.0);
}

static void PointToTriImpulse(
        POINT** pts,
        double* nor,
        double* w,
        double dist,
        MotionState mstate,
        double root)
{
	if (debugging("collision"))
	    CollisionSolver3d::pt_to_tri++;
    
    double k      = CollisionSolver3d::getSpringConstant();
	double m      = CollisionSolver3d::getPointMass();
	double dt     = CollisionSolver3d::getTimeStepSize();
	double mu     = CollisionSolver3d::getFrictionConstant(); 
	double h      = CollisionSolver3d::getFabricThickness();
	double cr     = CollisionSolver3d::getRestitutionCoef();

	dist   = h - dist; //overlap with fabric thickness

	STATE *sl[4];
	for (int i = 0; i < 4; ++i)
	    sl[i] = (STATE*)left_state(pts[i]);

    double v_rel[3];
	for (int i = 0; i < 3; ++i)
	{
	    v_rel[i] = sl[3]->avgVel[i];
        v_rel[i] -= w[0]*sl[0]->avgVel[i] + w[1]*sl[1]->avgVel[i] + w[2]*sl[2]->avgVel[i]; 
	}

	double vn = Dot3d(v_rel,nor);
	double vt = sqrt(Dot3d(v_rel,v_rel) - sqr(vn));
    
    /*
	if (Dot3d(v_rel, v_rel) > sqr(vn))
	    vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	else
	    vt = 0.0;
    */

    //Point and Triangle are approaching each other (vn < 0.0):
    //      apply inelastic impulse
    //Point and Triangle are seperating from each other (vn > 0.0):
    //      apply elastic impulse
    
    double sum_w = 0.0;
	double impulse = 0.0;
    double friction_impulse = 0.0;
	double rigid_impulse[2] = {0.0};

    if (mstate == MotionState::MOVING)
    {
        //Apply one or the other for collisions, NOT BOTH
        dt = root;
        if (vn < 0.0)
            PointToTriInelasticImpulse(vn,pts,&impulse,rigid_impulse,w,&sum_w);
        else if (vn * dt < 0.1 * dist)
            PointToTriElasticImpulse(vn,dist,pts,&impulse,rigid_impulse,dt,m,k,cr);
    }
    else
    {
        //Apply both if needed for proximity
        if (vn < 0.0)
            PointToTriInelasticImpulse(vn,pts,&impulse,rigid_impulse,w,&sum_w);
        if (vn * dt < 0.1 * dist)
            PointToTriElasticImpulse(vn,dist,pts,&impulse,rigid_impulse,dt,m,k,cr);
        if (vt > ROUND_EPS)
        {
            double delta_vt = 0.5*vt;
            if (fabs(mu*impulse) < vt)
                delta_vt = fabs(mu*impulse);

            //TODO: check sign of update
            friction_impulse = delta_vt;
            //sl[i]->friction[j] += W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
            //sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
        }
    }


    double m_impulse;
    double f_impulse;

    if (fabs(sum_w) < MACH_EPS)
    {
	    m_impulse = impulse;
        f_impulse = friction_impulse;
    }
    else
    {
	    m_impulse = 2.0 * impulse / (1.0 + Dot3d(w,w));
	    f_impulse = 2.0 * friction_impulse / (1.0 + Dot3d(w,w));
    }

    //uncomment the following the debugging purpose
    if (debugging("CollisionImpulse"))
    if (fabs(m_impulse) > 0.0){
        printf("real PointToTri collision, dist = %e\n",dist);
        printf("vt = %f, vn = %f, dist = %f\n",vt,vn,dist);
        printf("v_rel = %f %f %f\n",v_rel[0],v_rel[1],v_rel[2]);
        printf("nor = %f %f %f\n",nor[0],nor[1],nor[2]);
        printf("m_impulse = %f, impulse = %f, w = [%f %f %f]\n",
            m_impulse,impulse,w[0],w[1],w[2]);
        printf("dt = %f, root = %f\n",dt,root);
        printf("k = %f, m = %f\n",k,m);
        printf("x_old:\n");
        for (int i = 0; i < 4; ++i){
            STATE* sl1 = (STATE*)left_state(pts[i]);
            printf("%f %f %f\n",sl1->x_old[0],sl1->x_old[1],sl1->x_old[2]);
        }
        printf("x_new:\n");
        for (int i = 0; i < 4; ++i){
            printf("%f %f %f\n",Coords(pts[i])[0],Coords(pts[i])[1],Coords(pts[i])[2]);
        }
        printf("avgVel:\n");
        for (int i = 0; i < 4; ++i){
            STATE* sl1 = (STATE*)left_state(pts[i]);
            printf("%f %f %f\n",sl1->avgVel[0],sl1->avgVel[1],sl1->avgVel[2]);
        }
    }
    ////////////////////////////////////////////////////

    //TODO: This doesn't look right.
    //      Also should be its own function.
	if (isRigidBody(pts[0]) && isRigidBody(pts[1]) && 
	    isRigidBody(pts[2]) && isRigidBody(pts[3]))
	{
	    if (isMovableRigidBody(pts[0]))
            SpreadImpactZoneImpulse(pts[0], rigid_impulse[0], nor);
	    if (isMovableRigidBody(pts[3]))
            SpreadImpactZoneImpulse(pts[3], -1.0 * rigid_impulse[1], nor);
	    return;
	}

    std::vector<double> W = {w[0],w[1],w[2],-1.0};
    std::vector<double> R = {rigid_impulse[0],rigid_impulse[0],
                             rigid_impulse[1],rigid_impulse[1]};

    //TODO: should probably be in own function together with
    //      the elastic and inelastic impulse functions above
	for (int i = 0; i < 4; ++i)
	{
        if (!isStaticRigidBody(pts[i]))
        {
            sl[i]->collsn_num++;

            double t_impulse = m_impulse;
            double t_friction_impulse = f_impulse;
            if (isMovableRigidBody(pts[i]))
            {
                t_impulse = R[i];
            }

            for(int j = 0; j < 3; ++j)
            {
                sl[i]->collsnImpulse[j] += W[i]*t_impulse*nor[j];

                if (vt > ROUND_EPS)
                    sl[i]->friction[j] += W[i]*t_friction_impulse*(v_rel[j] - vn*nor[j])/vt;

                /*
                //Friction only applied for proximity, not collisions.
                if (mstate == MotionState::MOVING)
                    continue;

                if (fabs(vt) > ROUND_EPS)
                {
                    double delta_vt = vt;
                    if (fabs(lambda*t_impulse) < vt)
                        delta_vt = fabs(lambda*t_impulse);
   
                    //TODO: check sign of update
                    //sl[i]->friction[j] += W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;
                    sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;

                    //double delta_vt = vt;
                    //if (fabs(lambda*t_impulse) < vt)
                      //  delta_vt = lambda*t_impulse;

                    //sl[i]->friction[j] -= W[i]*delta_vt*(v_rel[j] - vn*nor[j])/vt;

                    //double frcoef = std::max(-fabs(lambda*W[i]*t_impulse/vt), -1.0);
                    //sl[i]->friction[j] += frcoef*(v_rel[j] - vn*nor[j]);
                    //double delta_vt = std::min(1.0,fabs(lambda*delta_vn/vt));
                    //sl[i]->friction[j] -= delta_vt*(v_rel[j] - vn*nor[j]);
                }
                */
            }
        }
        else
        {
            for (int j = 0; j < 3; ++j)
            {
                sl[i]->collsnImpulse[j] = 0.0;
                sl[i]->friction[j] = 0.0;
            }
        }
	}

	if (debugging("CollisionImpulse"))
	for (int i = 0; i < 4; ++i)
    {
	    printf("pt[%d], collsnImp = [%f %f %f], friction = [%f %f %f]\n", i,
        sl[i]->collsnImpulse[0],sl[i]->collsnImpulse[1],sl[i]->collsnImpulse[2],
        sl[i]->friction[0],sl[i]->friction[1],sl[i]->friction[2]);
	}

	for (int kk = 0; kk < 4; kk++)
	for (int j = 0; j < 3; ++j)
    {
        if (std::isnan(sl[kk]->collsnImpulse[j]) ||
		std::isinf(sl[kk]->collsnImpulse[j]))
        {
            printf("PointToTri: sl[%d]->impl[%d] = nan\n",kk,j);
            for (int i = 0; i < 4; ++i)
            {
                printf("points[%d] = %p\n",i,(void*)pts[i]);
                printf("coords = [%f %f %f]\n",Coords(pts[i])[0],
                    Coords(pts[i])[1],Coords(pts[i])[2]);
            }
            printf("w = [%f %f %f]\nnor = [%f %f %f]\ndist = %f\n",
            w[0],w[1],w[2],nor[0],nor[1],nor[2],dist);
            printf("v_rel = [%f %f %f]\n",v_rel[0],v_rel[1],v_rel[2]);
            clean_up(ERROR);
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
    if (isStaticRigidBody(pts[3]) ||
       (isStaticRigidBody(pts[0]) && isStaticRigidBody(pts[1])
        && isStaticRigidBody(pts[2])))
    {
        *impulse += vn;
        rigid_impulse[0] += vn;
        rigid_impulse[1] += vn;
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1]) 
            && isMovableRigidBody(pts[2]) && isMovableRigidBody(pts[3]))
    {
        double m1 = total_mass(pts[0]->hs);
        double m2 = total_mass(pts[3]->hs);
        rigid_impulse[0] += vn * m2 / (m1 + m2);
        rigid_impulse[1] += vn * m1 / (m1 + m2);
    }
    else if (isMovableRigidBody(pts[0]) && isMovableRigidBody(pts[1])
            && isMovableRigidBody(pts[2]))
    {
        *impulse += 0.5 * vn;
        rigid_impulse[0] += 0.5 * vn;
    }
    else if (isMovableRigidBody(pts[3]))
    {
        *impulse += 0.5 * vn;
        rigid_impulse[1] += 0.5 * vn;
    }
    else
    {
        *impulse += 0.5 * vn;
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
        double dist,
        POINT** pts,
        double* impulse,
        double* rigid_impulse,
        double dt, double m,
        double k, double cor)
{
    if (isRigidBody(pts[0]) && isRigidBody(pts[1]) &&
    isRigidBody(pts[2]) && isRigidBody(pts[3]))
    {
        rigid_impulse[0] *= 1.0 + cor;
        rigid_impulse[1] *= 1.0 + cor;
    }
    else
    {
        double tmp = -1.0*std::min(dt*k*dist/m, (0.1*dist/dt - vn));
        *impulse += tmp;
        rigid_impulse[0] += tmp;
        rigid_impulse[1] += tmp;
    }
}

static void unsort_surface_point(SURFACE *surf)
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

