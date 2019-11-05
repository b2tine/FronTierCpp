#include "collid.h"

static void PointToTriImpulse(POINT**,double*,double*,double,MotionState,double);

static void PointToTriInelasticImpulse(double,POINT**,double*,double*,double*,double*);

static void PointToTriElasticImpulse(double,double,POINT**,double*,double*,
                                     double,double,double,double);

static void EdgeToEdgeImpulse(POINT**,double*,double,double,double,MotionState,double);

static void EdgeToEdgeInelasticImpulse(double,POINT**,double*,double*,double*);

static void EdgeToEdgeElasticImpulse(double,double,POINT**,double*,double*,
                                     double,double,double,double);

static void SpreadImpactZoneImpulse(POINT*, double, double*);


void PointToTriProximityImpulse(
        POINT** pts,
        double* nor,
        double* w,
        double dist)
{
    MotionState mstate = MotionState::STATIC;
    PointToTriImpulse(pts,nor,w,dist,mstate,-1.0);
}

void PointToTriCollisionImpulse(
        POINT** pts,
        double* nor,
        double* w,
        double dist,
        double dt)
{
    MotionState mstate = MotionState::MOVING;
    PointToTriImpulse(pts,nor,w,dist,mstate,dt);
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
    {
        if (fabs(m_impulse) > 0.0)
        {
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
            //TODO: should probably move this somewhere else
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

void EdgeToEdgeCollisionImpulse(
        POINT** pts,
        double* nor,
        double a,
        double b,
        double dist,
        double dt)
{
    MotionState mstate = MotionState::MOVING;
    EdgeToEdgeImpulse(pts,nor,a,b,dist,mstate,dt);
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
    {
        if (fabs(m_impulse) > 0.0)
        {
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
            //TODO: should probably move this somewhere else
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

static void SpreadImpactZoneImpulse(POINT* p, double impulse, double* nor)
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







