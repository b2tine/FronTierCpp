#include "Proximity.h"
#include "collid.h"

static std::unique_ptr<Proximity> TriToTri(const TRI*,const TRI*,double);
static std::unique_ptr<Proximity> TriToBond(const TRI*,const BOND*,double);
static std::unique_ptr<Proximity> BondToBond(const BOND*,const BOND*,double);

static std::unique_ptr<Proximity> StaticPointToTri(POINT**,double);
static std::unique_ptr<Proximity> StaticEdgeToEdge(POINT**,double);

static std::unique_ptr<Collision> MovingTriToTri(const TRI*,const TRI*,double);
static std::unique_ptr<Collision> MovingTriToBond(const TRI*,const BOND*,double);
static std::unique_ptr<Collision> MovingBondToBond(const BOND*,const BOND*,double);

static std::unique_ptr<Collision> MovingPointToTri(POINT**,double);
static std::unique_ptr<Collision> MovingEdgeToEdge(POINT**,double);

static bool isCoplanar(POINT**,double,double*);

static std::unique_ptr<Collision> KineticPointToTri(POINT**,double,double,double);
static std::unique_ptr<Collision> KineticEdgeToEdge(POINT**,double,double,double);

static bool PointToTri(POINT**,double,double*,double*,double*);
static bool EdgeToEdge(POINT**,double,double*,double*,double*,double*);


static void UpdateAverageVelocity(POINT** pts);
static void UpdatePostCollisionState(POINT** pts, double dt);

//static void SaveState(POINT** pts);
static void UpdateState(POINT** pts, double dt);
static void RestorePrevState(POINT** pts);



std::unique_ptr<Proximity> checkProximity(const CD_HSE* a, const CD_HSE* b, double tol)
{
	const CD_BOND *cd_b1, *cd_b2;
	const CD_TRI  *cd_t1, *cd_t2;

	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
            return {};
	    return TriToTri(t1,t2,tol);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return BondToBond(b1,b2,tol);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1,tol);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return TriToBond(t1,b1,tol);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
}

std::unique_ptr<Proximity> TriToTri(const TRI* tri1, const TRI* tri2, double tol)
{
    //TODO: ensure AABB::getCandidates() eliminates
    //      this possibility, and remove.
	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
        //tris share a point
	    if (Point_of_tri(tri1)[i] == Point_of_tri(tri2)[j])
            return {};
	}

	POINT* pts[4];
    std::vector<std::unique_ptr<Proximity>> proximities;

	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
	{
	    const TRI* tmp_tri1 = (k == 0) ? tri1 : tri2;
	    const TRI* tmp_tri2 = (k == 0) ? tri2 : tri1;

	    pts[3] = Point_of_tri(tmp_tri1)[i];
        for (int j = 0; j < 3; ++j)
            pts[j] = Point_of_tri(tmp_tri2)[j];
	
        //TODO: ensure AABB::getCandidates() eliminates
        //      this possibility, and remove.
        //
	    //Don't check a point against the triangle that contains it
        if (pts[0] == pts[3] ||
            pts[1] == pts[3] || pts[2] == pts[3])
            continue;

        std::unique_ptr<Proximity> proximity = StaticPointToTri(pts,tol);
        if (proximity)
            proximities.push_back(std::move(proximity));
	}

	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri1)[i];
	    pts[1] = Point_of_tri(tri1)[(i+1)%3];
	
        for (int j = 0; j < 3; ++j)
        {
		    pts[2] = Point_of_tri(tri2)[j];
		    pts[3] = Point_of_tri(tri2)[(j+1)%3];
		
            //TODO: ensure AABB::getCandidates() eliminates
            //      this possibility, and remove.
            //
            //Don't check edges that share an endpoint
            if (pts[0] == pts[2] || pts[0] == pts[3] ||
                pts[1] == pts[2] || pts[1] == pts[3])
                continue;

            std::unique_ptr<Proximity> proximity = StaticEdgeToEdge(pts,tol);
            if (proximity)
                proximities.push_back(std::move(proximity));
	    }
    }

    double min_dist = HUGE;
    std::unique_ptr<Proximity> closest;

    std::vector<unique_ptr<Proximity>>::iterator it;
    for (it = proximities.begin(); it < proximities.end(); ++it)
    {
        if ((*it)->dist < min_dist)
        {
            min_dist = (*it)->dist;
            closest = std::move(*it);
        }
    }

    return closest;
}

std::unique_ptr<Proximity> TriToBond(const TRI* tri,const BOND* bd, double h)
{
    //TODO: ensure AABB::getCandidates() eliminates
    //      this possibility, and remove.
    //
	/* do not consider bond point that is a tri vertex */
	for (int i = 0; i < 3; ++i)
	{
	    if (Point_of_tri(tri)[i] == bd->start ||
            Point_of_tri(tri)[i] == bd->end)
            return {};
	}

	POINT* pts[4];
    std::vector<std::unique_ptr<Proximity>> proximities;
	
    for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

	//detect proximity of triangle to bond end points.
    for (int i = 0; i < 2; ++i)
    {
        pts[3] = (i == 0) ? bd->start : bd->end;
        
        st::unique_ptr<Proximity> proximity = StaticPointToTri(pts,h);
        if (proximity)
            proximities.push_back(std::move(proximity));
    }

	//detect proximity of each triangle edge with the bond edge
	pts[2] = bd->start;
    pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];

        std::unique_ptr<Proximity> proximity = StaticEdgeToEdge(pts,h);
        if (proximity)
            proximities.push_back(std::move(proximity));
	}

    double min_dist = HUGE;
    std::unique_ptr<Proximity> closest;

    std::vector<unique_ptr<Proximity>>::iterator it;
    for (it = proximities.begin(); it < proximities.end(); ++it)
    {
        if ((*it)->dist < min_dist)
        {
            min_dist = (*it)->dist;
            closest = std::move(*it);
        }
    }

    return closest;
}

std::unique_ptr<Proximity> BondToBond(const BOND* b1, const BOND* b2, double tol)
{
	POINT* pts[4];

	pts[0] = b1->start; pts[1] = b1->end;
	pts[2] = b2->start; pts[3] = b2->end;

    //TODO: ensure AABB::getCandidates() eliminates
    //      this possibility, and remove.
	if (pts[0] == pts[2] || pts[0] == pts[3]
        || pts[1] == pts[2] || pts[1] == pts[3])
    {
        return {};
    }
    
    std::uniqe_ptr<Proximity> proximity = StaticEdgeToEdge(pts,tol);
    return proximity;
}

static std::unique_ptr<Proximity> StaticPointToTri(POINT** pts, double h)
{
    double nor[3];
    double w[3];
    double dist;

    std::unique_ptr<Proximity> proximity;

    if (PointToTri(pts,h,nor,w,&dist))
    {
        proximity =
            std::unique_ptr<Proximity>(new PointTriProximity(pts,nor,a,b,dist));
    }

    return proximity;
}

static std::unique_ptr<Proximity> StaticEdgeToEdge(POINT** pts, double h)
{
    double nor[3];
    double a, b, dist;
    
    std::unique_ptr<Proximity> proximity;

    if (EdgeToEdge(pts,h,nor,&a,&b,&dist))
    {
        proximity =
            std::unique_ptr<Proximity>(new EdgeEdgeProximity(pts,nor,a,b,dist));
    }

    return proximity;
}

std::unique_ptr<Collision> checkCollision(const CD_HSE* a, const CD_HSE* b, double tol)
{
	const CD_TRI  *cd_t1, *cd_t2;
	const CD_BOND *cd_b1, *cd_b2;

	if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) && 
	    (cd_t2 = dynamic_cast<const CD_TRI*>(b)))
	{
	    TRI* t1 = cd_t1->m_tri;
	    TRI* t2 = cd_t2->m_tri;
	    if ((t1->surf == t2->surf) && isRigidBody(a))
            return {};
	    return MovingTriToTri(t1,t2,tol);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) && 
	         (cd_b2 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    BOND* b2 = cd_b2->m_bond;
	    return MovingBondToBond(b1,b2,tol);
	}
	else if ((cd_b1 = dynamic_cast<const CD_BOND*>(a)) &&
		 (cd_t1 = dynamic_cast<const CD_TRI*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1,tol);
	}
	else if ((cd_t1 = dynamic_cast<const CD_TRI*>(a)) &&
                 (cd_b1 = dynamic_cast<const CD_BOND*>(b)))
	{
	    BOND* b1 = cd_b1->m_bond;
	    TRI* t1  = cd_t1->m_tri;
	    return MovingTriToBond(t1,b1,tol);
	}
	else
	{
	    std::cout<<"This case has not been implemented"<<std::endl;
	    clean_up(ERROR);
	}
}

std::unique_ptr<Collision> MovingTriToTri(const TRI* a,const TRI* b, double h)
{
    //TODO: ensure AABB::getCandidates() eliminates
    //      this possibility, and remove.
	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    if (Point_of_tri(a)[i] == Point_of_tri(b)[j])
            return {};
	}

	POINT* pts[4];
    std::vector<std::unique_ptr<Collision>> collisions;

	//detect point to tri collision
	for (int k = 0; k < 2; ++k)
	for (int i = 0; i < 3; ++i)
    {
	    const TRI* tmp_tri1 = (k == 0) ? a : b;
	    const TRI* tmp_tri2 = (k == 0) ? b : a;

	    for (int j = 0; j < 3; ++j)
            pts[j] = Point_of_tri(tmp_tri1)[j];
        pts[3] = Point_of_tri(tmp_tri2)[i];

        //TODO: ensure AABB::getCandidates() eliminates
        //      this possibility, and remove.
        //
	    //Don't consider point against the triangle it belongs to
	    if (pts[0] == pts[3] ||
            pts[1] == pts[3] || pts[2] == pts[3])
            continue; 
        
        std::unique_ptr<Collision> collsn = MovingPointToTri(pts,h);
        if (collsn)
            collisions.push_back(std::move(collsn));
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
		
            //TODO: ensure AABB::getCandidates() eliminates
            //      this possibility, and remove.
            //
		    //Don't consider edges with a shared enpoint
            if (pts[0] == pts[2] || pts[0] == pts[3] ||
                pts[1] == pts[2] || pts[1] == pts[3])
                continue;
    
            std::unique_ptr<Collision> collsn = MovingEdgeToEdge(pts,h);
            if (collsn)
                collisions.push_back(std::move(collsn));
	    }
    }

    double min_dist = HUGE;
    std::unique_ptr<Collision> closest;

    std::vector<unique_ptr<Collision>>::iterator it;
    for (it = collisions.begin(); it < collisions.end(); ++it)
    {
        if ((*it)->dist < min_dist)
        {
            min_dist = (*it)->dist;
            closest = std::move(*it);
        }
    }

    return closest;
}

std::unique_ptr<Collision> MovingTriToBond(const TRI* tri,const BOND* bd, double h)
{
	/* do not consider bond point that is a tri vertex */
	for (int i = 0; i < 3; ++i)
	{
        //TODO: ensure AABB::getCandidates() eliminates
        //      this possibility, and remove.
        //
        //do not consider bond point that is a tri vertex
	    if (Point_of_tri(tri)[i] == bd->start ||
            Point_of_tri(tri)[i] == bd->end)
            return {};
	}

	POINT* pts[4];
    std::vector<std::unique_ptr<Collision>> collisions;

	for (int i = 0; i < 3; ++i)
	    pts[i] = Point_of_tri(tri)[i];

    //detect proximity of triangle to bond end points.
    for (int i = 0; i < 2; ++i)
    {
        pts[3] = (i == 0) ? bd->start : bd->end;

         std::unique_ptr<Collision> collsn = MovingPointToTri(pts,h);
         if (collsn)
             collisions.push_back(std::move(collsn));
    }

    /* detect collision of each of tri edge w.r.t to bond */
	pts[2] = bd->start;
	pts[3] = bd->end;
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri)[i];
	    pts[1] = Point_of_tri(tri)[(i+1)%3];
        
        std::unique_ptr<Collision> collsn = MovingEdgeToEdge(pts,h);
        if (collsn)
            collisions.push_back(std::move(collsn));
	}

    double min_dist = HUGE;
    std::unique_ptr<Collision> closest;

    std::vector<unique_ptr<Collision>>::iterator it;
    for (it = collisions.begin(); it < collisions.end(); ++it)
    {
        if ((*it)->dist < min_dist)
        {
            min_dist = (*it)->dist;
            closest = std::move(*it);
        }
    }

    return closest;
}

std::unique_ptr<Collision> MovingBondToBond(const BOND* b1, const BOND* b2, double tol)
{
	POINT* pts[4];

	pts[0] = b1->start; pts[1] = b1->end;
	pts[2] = b2->start; pts[3] = b2->end;

    //TODO: ensure AABB::getCandidates() eliminates
    //      this possibility, and remove.
    //
	/* do not consider two bonds that share a common point */
	if (pts[0] == pts[2] || pts[0] == pts[3]
        || pts[1] == pts[2] || pts[1] == pts[3])
    {
        return {};
    }

    return MovingEdgeToEdge(pts,tol);
}

static std::unique_ptr<Collision> MovingPointToTri(POINT** pts, double h)
{
	double maxdt = CollisionSolver3d::getTimeStepSize();
	double dt[4] = {-1,-1,-1,maxdt};

	if (isCoplanar(pts,maxdt,dt))
    {
        for (int i = 0; i < 4; ++i)
        {
            if (dt[i] < 0.0)
                continue;

            std::unique_ptr<Collision> collsn = KineticPointToTri(pts,h,dt[i]);
            if (collsn)
                return collsn;
        }
	}

    return {};
}

static std::unique_ptr<Collision> MovingEdgeToEdge(POINT** pts, double h)
{
    double maxdt = CollisionSolver3d::getTimeStepSize();
	double dt[4] = {-1,-1,-1,maxdt};

	if (isCoplanar(pts,maxdt,dt))
    {
        for (int i = 0; i < 4; ++i)
        {
            if (dt[i] < 0.0)
                continue;

            std::unique_ptr<Collision> collsn = KineticEdgeToEdge(pts,h,dt[i],maxdt);
            if (collsn)
                return collsn;
        }
    }

    return {};
}

static void isCoplanarHelper(double* s[], double v[][3])
{
    v[0][0] = s[0][0];         v[0][1] = s[0][1];         v[0][2] = s[0][2];
	v[1][0] = s[1][0]-s[0][0]; v[1][1] = s[1][1]-s[0][1]; v[1][2] = s[1][2]-s[0][2];
	v[2][0] = s[2][0]-s[0][0]; v[2][1] = s[2][1]-s[0][1]; v[2][2] = s[2][2]-s[0][2];
	v[3][0] = s[3][0]-s[0][0]; v[3][1] = s[3][1]-s[0][1]; v[3][2] = s[3][2]-s[0][2];
}

static bool isCoplanar(POINT* pts[], double maxdt, double* roots)
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
	for (int i = 0; i < 3; ++i)
    {
        roots[i] = roots[i]-MACH_EPS;
        if (roots[i] < 0.0 || roots[i] > maxdt)
            roots[i] = -1;
	}
	//sort the roots
	if (roots[0] > roots[1])
	    std::swap(roots[0], roots[1]);
	if (roots[0] > roots[2])
	    std::swap(roots[0], roots[2]);
	if (roots[1] > roots[2])
	    std::swap(roots[1], roots[2]);

	if (roots[0] > MACH_EPS ||
        roots[1] > MACH_EPS || roots[2] > MACH_EPS)
    {
        return true;
    }

    return false;
}

static std::unique_ptr<Collision> KineticPointToTri(
        POINT** pts,
        double h,
        double dt,
        double maxdt)
{
    double nor[3];
    double w[3];
    double dist;

    STATE* sl;
    for (int j = 0; j < 4; ++j)
    {
        sl = (STATE*) left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k] + dt*sl->avgVel[k];
    }

    std::unique_ptr<Collision> collision;

    if (PointToTri(pts,h,nor,w,&dist))
    {
        collision =
            std::unique_ptr<Collision>(new PointTriCollision(pts,nor,w,dist,dt,maxdt));

    }

    //restore coordinates of points
    for (int j = 0; j < 4; ++j)
    {
        sl = (STATE*)left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k];
    }

    return collision;
}

static std::unique_ptr<Collision> KineticEdgeToEdge(
        POINT** pts,
        double h,
        double dt,
        double maxdt)
{
    double nor[3];
    double a, b, dist;

    STATE* sl;
    for (int j = 0; j < 4; ++j)
    {
        sl = (STATE*) left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k] + dt*sl->avgVel[k];
    }

    std::unique_ptr<Collision> collision;

    if (EdgeToEdge(pts,h,nor,&a,&b,&dist))
    {
        collision =
            std::unique_ptr<Collision>(new EdgeEdgeCollision(pts,nor,a,b,dist,dt,maxdt));

    }

    //restore coordinates of points
    for (int j = 0; j < 4; ++j)
    {
        sl = (STATE*)left_state(pts[j]);
        for (int k = 0; k < 3; ++k)
            Coords(pts[j])[k] = sl->x_old[k];
    }

    return collision;
}

static bool PointToTri(POINT** pts, double h, double* nor, double* w, double* dist)
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

    double distance = Dot3d(x34,tri_nor);
    if (distance < 0.0)
    {
        distance = fabs(dist);
        scalarMult(-1.0,tri_nor,tri_nor);
    }

    if (distance > h)
        return false;

	double W[3];
    W[0] = (Dot3d(x23,x23)*Dot3d(x13,x43)-Dot3d(x13,x23)*Dot3d(x23,x43))/det;
    W[1] = (Dot3d(x13,x13)*Dot3d(x23,x43)-Dot3d(x13,x23)*Dot3d(x13,x43))/det;
    W[2] = 1.0 - W[0] - W[1];
    
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
        if (W[i] < -1.0*eps || W[i] > 1.0 + eps)
            return false;
    }

    for (int i = 0; i < 3; ++i)
    {
        w[i] = W[i];
        nor[i] = tri_nor[i];
    }
    *dist = distance;

    return true;

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
}

//For details of this implementation see:
//http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
static bool EdgeToEdge(
        POINT** pts,
        double h,
        double* nor,
        double* A,
        double* B,
        double* dist)
{
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

    //Compute the solution and obtain the closest pair of points.
    if (sN == sD)
        sC = 1.0;
    else
        sC = fabs(sN) < ROUND_EPS ? 0.0 : sN/sD;

    if (tN == tD)
        tC = 1.0;
    else
        tC = fabs(tN) < ROUND_EPS ? 0.0 : tN/tD;
    
    
    //TODO: Make sure nor is pointing in the correct direction.
    
    //"normal vector" always points from
    //the x12 edge to the x34 edge
    
	for (int j = 0; j < 3; ++j)
	{
	    nor[j]  = (1.0 - tC)*Coords(pts[2])[j] + tC*Coords(pts[3])[j];
	    nor[j] -= (1.0 - sC)*Coords(pts[0])[j] + sC*Coords(pts[1])[j];
    }

    double distance = Mag3d(nor);
    if (distance > h)
        return false;

    if (distance > ROUND_EPS)
    {
        scalarMult(1.0/(distance),nor,nor);
        *A = sC;
        *B = tC;
        *dist = distance;
    }
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

    return true;
}

/*
static void PointToLine(POINT** pts, double* a)
static void PointToLine(POINT** pts, double& a)
{
//	x1 -----projP---- x2
//		  |
//		  |  dist
//		  *x3
	
    double x12[3], x13[3];
	Pts2Vec(pts[0],pts[1],x12);
    Pts2Vec(pts[0],pts[2],x13);
    // *a = Dot3d(x13,x12)/Dot3d(x12,x12);
    a = Dot3d(x13,x12)/Dot3d(x12,x12);
}
*/


PointTriProximity::PointTriProximity(POINT** Pts, double* Nor,
                                    double* W, double Dist)
    : pts{Pts}, dist{Dist}
{
    for (int i = 0; i < 3; ++i)
    {
        w[i] = W[i];
        nor[i] = Nor[i];
    }
}

void PointTriProximity::computeImpulse()
{
    PointToTriProximityImpulse(pts,nor,w,dist);
}

void PointTriProximity::updateAverageVelocity()
{
    UpdateAverageVelocity(pts);
}

void PointTriProximity::computePostCollisionImpulse(double dt)
{
    PointToTriPostCollisionProximityImpulse(pts,nor,w,dist,dt);
}

/*
void PointTriProximity::updatePostCollisionState(double dt)
{
    UpdateState(pts,dt);
}
*/

EdgeEdgeProximity::EdgeEdgeProximity(POINT** Pts, double* Nor,
                                    double A, double B, double Dist)
    : pts{Pts}, a{A}, b{B}, dist{Dist}
{
    for (int i = 0; i < 3; ++i)
        nor[i] = Nor[i];
}

void EdgeEdgeProximity::computeImpulse()
{
    EdgeToEdgeProximityImpulse(pts,nor,a,b,dist);
}

void EdgeEdgeProximity::updateAverageVelocity()
{
    UpdateAverageVelocity(pts);
}

void EdgeEdgeProximity::computePostCollisionImpulse(double dt)
{
    EdgeToEdgePostCollisionProximityImpulse(pts,nor,a,b,dist,dt);
}

/*
void EdgeEdgeProximity::updatePostCollisionState(double dt)
{
    UpdateState(pts,dt);
}
*/

PointTriCollision::PointTriCollision(POINT** Pts, double* Nor, double* W,
                                    double Dist, double Dt, double MaxDt)
    : pts{Pts}, dist{Dist}, dt{Dt}, maxdt{MaxDt}
{
    for (int i = 0; i < 3; ++i)
    {
        w[i] = W[i];
        nor[i] = Nor[i];
    }
}

void PointTriCollision::computeImpulse()
{
    PointToTriCollisionImpulse(pts,nor,w,dist,dt);
}

void PointTriCollision::updateState()
{
    UpdateState(pts,dt);
}

void PointTriCollision::checkNewStateProximity(double tol)
{
    std::unique_ptr<Proximity> proximity = StaticPointToTri(pts,tol);
    if (proximity)
    {
        double avg_dt = 0.5*(dt + maxdt);
        proximity->computePostCollisionImpulse(avg_dt);
        UpdateNewState(pts,avg_dt);
        //proximity->updatePostCollisionState(avg_dt);
    }
}

void PointTriCollision::mergeImpactZones()
{
    CreateImpZone(pts,4);
}

void PointTriCollision::restorePrevState()
{
    RestorePrevState(pts);
}

EdgeEdgeCollision::EdgeEdgeCollision(POINT** Pts, double* Nor, double A,
                                double B, double Dist, double Dt, double MaxDt)
    : pts{Pts}, a{A}, b{B}, dist{Dist}, dt{DT}, maxdt{MaxDt}
{
    for (int i = 0; i < 3; ++i)
        nor[i] = Nor[i];
}

void EdgeEdgeCollision::computeImpulse()
{
    EdgeToEdgeCollisionImpulse(pts,nor,a,b,dist,dt);
}

void EdgeEdgeCollision::updateState()
{
    UpdateState(pts,dt);
}

void EdgeEdgeCollision::checkNewStateProximity(double tol)
{
    std::unique_ptr<Proximity> proximity = StaticEdgeToEdge(pts,tol);
    if (proximity)
    {
        double avg_dt = 0.5*(dt + maxdt);
        proximity->computePostCollisionImpulse(avg_dt);
        UpdateNewState(pts,avg_dt);
        //proximity->updatePostCollisionState(avg_dt);
    }
}

void EdgeEdgeCollision::mergeImpactZones()
{
    CreateImpZone(pts,4);
}

void EdgeEdgeCollision::restorePrevState()
{
    RestorePrevState(pts);
}

//TODO: make method of Proximity class?
//
//Call after all impulses have been computed;
//see CollisionSolver3d::processProximityCandidates()
static void UpdateAverageVelocity(POINT** pts)
{
	POINT *p;
	STATE *sl;
	
    for (int i = 0; i < 4; ++i)
    {
        p = pts[i];
        if (isStaticRigidBody(p))
            continue;

        sl = (STATE*) left_state(p);

        //sl->has_collsn = true;
    
        //TODO: divide only the inelastic collision impulse by sl->collsn_num
        //      or multiply by some relaxation parameter (e.g. 0.25)?
        
        //TODO: need to divide by sl->collsn_num at all?

        //TODO: does friction impulse also need to be divided
        //      by the number of collisions?
        
        for (int k = 0; k < 3; ++k)
        {
            sl->avgVel[k] += sl->collsnImpulse[k]/sl->collsn_num;
            //sl->avgVel[k] += sl->friction[k]/sl->collsn_num;
            //sl->avgVel[k] += sl->collsnImpulse[k];
            sl->avgVel[k] += sl->friction[k];
            
            sl->collsnImpulse[k] = 0.0;
            sl->friction[k] = 0.0;
            sl->collsn_num = 0;

            //save for Collision RestoreState()
            sl->avgVel_old[k] = sl->avgVel[k];
            
            if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
            {
                printf("inf/nan vel[%d]: impulse = %f, friction = %f, collsn_num = %d\n",
                        k,sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
                clean_up(ERROR);
            }
        }
    
        /*
        //TODO: do we really need to differentiate now?
        // test for RG
        if (sl->collsn_num_RG > 0)
        {
            for (int k = 0; k < 3; ++k)
                sl->avgVel[k] += sl->collsnImpulse_RG[k]/sl->collsn_num_RG;
            sl->collsn_num_RG = 0;
        }
        */

    }
	
    /*
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects
    */
}

/*
static void SaveState(POINT** pts)
{
	POINT *p;
	STATE *sl;

    for (int i = 0; i < 4; ++i)
    {
        p = pts[i];
        sl = (STATE*) left_state(p);

        for (int k = 0; k < 3; ++k)
        {
            sl->x_old[k] = Coords(p)[k];
            sl->avgVel_old[k] = sl->avgVel[k];
        }
    }
}
*/

//TODO: make method of Collision class?
//
//void UpdateAverageVelocityCollision(POINT** pts, double dt)
static void UpdateState(POINT** pts, double dt)
{
	POINT *p;
	STATE *sl;
	
    for (int i = 0; i < 4; ++i)
    {
        p = pts[i];
        if (isStaticRigidBody(p))
            continue;

        sl = (STATE*) left_state(p);

        //sl->has_collsn = true;
    
        //TODO: divide only the inelastic collision impulse by sl->collsn_num
        //      or multiply by some relaxation parameter (e.g. 0.25)?
        
        //TODO: need to divide by sl->collsn_num at all?
        for (int k = 0; k < 3; ++k)
        {
            sl->avgVel[k] += sl->collsnImpulse[k];
            //sl->avgVel[k] += sl->collsnImpulse[k]/sl->collsn_num;
            
            sl->collsnImpulse[k] = 0.0;

            //TODO: is this what we want?
            sl->collsn_num = 0;
            
            if (std::isinf(sl->avgVel[k]) || std::isnan(sl->avgVel[k])) 
            {
                printf("inf/nan vel[%d]: impulse = %f, collsn_num = %d\n",
                        k,sl->collsnImpulse[k],sl->collsn_num);
                clean_up(ERROR);
            }
            
            //compute new position and new effective velocity
            Coords(p)[k] = sl->x_old[k] + dt*sl->avgVel[k];
            sl->avgVel[k] = (Coords(p)[j] - sl->x_old[j])/dt;
        }

    
        /*
        //TODO: do we really need to differentiate now?
        // test for RG
        if (sl->collsn_num_RG > 0)
        {
            for (int k = 0; k < 3; ++k)
                sl->avgVel[k] += sl->collsnImpulse_RG[k]/sl->collsn_num_RG;
            sl->collsn_num_RG = 0;
        }
        */

    }
	
    /*
    if (getTimeStepSize() > 0.0)
	    updateImpactZoneVelocityForRG(); // test for moving objects
    */
}

//TODO: make method of Collision class?
//void RestoreAverageVelocity(POINT** pts)
static void RestorePrevState(POINT** pts)
{
    POINT *p;
	STATE *sl;
	
    for (int i = 0; i < 4; ++i)
    {
        p = pts[i];
        sl = (STATE*) left_state(p);

        for (int k = 0; k < 3; ++k)
        {
            Coords(p)[k] = sl->x_old[k];
            sl->avgVel[k] = sl->avgVel_old[k];
        }
    }

}


