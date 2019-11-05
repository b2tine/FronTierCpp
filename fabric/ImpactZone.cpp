#include <armadillo>
#include "collid.h"

POINT*& UF_Root(POINT* p)
{
	STATE* sl = (STATE*) left_state(p);
	return sl->impZone.root;
}

POINT*& UF_NextPoint(POINT* p)
{
	STATE* sl = (STATE*) left_state(p);
    return sl->impZone.next_pt;
}

POINT*& UF_Tail(POINT* p)
{
	STATE* sl = (STATE*) left_state(p);
	return sl->impZone.tail;
}

int& UF_Weight(POINT* p)
{
	STATE* sl = (STATE*) left_state(p);
	return sl->impZone.num_pts;
}

void UF_MakeDisjointSets(std::vector<CD_HSE*>& hselist)
{
	STATE* sl;
	POINT* pt;

    for (std::vector<CD_HSE*>::iterator it = hselist.begin();
            it < hselist.end(); ++it)
    {
        for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            sorted(pt) = NO;

            sl = (STATE*) left_state(pt);

            sl->impZone.root = pt;
            sl->impZone.next_pt = nullptr;
            sl->impZone.tail = pt;
            sl->impZone.num_pts = 1;
        }
    }
}

POINT* UF_FindSet(POINT* p)
{
	if (UF_Root(p) != p)
		UF_Root(p) = UF_FindSet(UF_Root(p));
	return UF_Root(p);
}

void UF_MergePoints(POINT* X, POINT* Y)
{
	POINT* PX = UF_FindSet(X);
	POINT* PY = UF_FindSet(Y);
	
    if (PX == PY)
        return;
	
    if (UF_Weight(PX) > UF_Weight(PY))
    {
	    //update root after merge
	    UF_Weight(PX) += UF_Weight(PY);
	    UF_Root(PY) = PX;

	    //link lists, update tail
	    UF_NextPoint(UF_Tail(PX)) = PY;
	    UF_Tail(PX) = UF_Tail(PY); 
	}
	else
    {
	    //update root after merge
	    UF_Weight(PY) += UF_Weight(PX);
	    UF_Root(PX) = PY;

	    //link lists, update tail
	    UF_NextPoint(UF_Tail(PY)) = PX;
	    UF_Tail(PY) = UF_Tail(PX); 
	}
}

//Note: num has default value of 4,
//and first has default value of false
void CreateImpactZone(POINT** pts, int num, bool first)
{
	for (int i = 0; i < num; ++i)
	{
	    for (int j = 0; j < i; ++j)
	    {
	        if ((!first) && (isMovableRigidBody(pts[i])
                || isMovableRigidBody(pts[j])))
            {
                continue;
            }
            mergePoint(pts[i],pts[j]); 
	    }
	}
    //TODO: Check that all points of the interface have not
    //      been merged into a single impact zone.
    //      If this occurs, then we should have a collision
    //      free state (double check this claim) and the collision
    //      handling iterations can be terminated.
}

// test function for creating impact zone for each movable RG
void CollisionSolver3d::createImpZoneForRG(const INTERFACE* intfc)
{
	SURFACE** s;
	TRI* tri;

	intfc_surface_loop(intfc, s)
	{
	    if (is_bdry(*s))
            continue;
	    
        if (!isMovableRigidBody(Point_of_tri(first_tri(*s))[0]))
            continue;

        surf_tri_loop(*s, tri)
	    {
    		CreateImpactZone(Point_of_tri(tri), 3, YES);
	    }
	}
}


void CollisionSolver3d::computeImpactZones()
{
    if (debugging("collision"))
        std::cout<<"Starting fail-safe (Impact Zone) method:\n";
	
	turnOnImpZone();

    int numImpactZones = 0;
    bool collision_free = false; 

    while (!collision_free)
    {
        aabbCollision();
        collisionCandidates.clear();
        collisionCandidates = abt_collision->getCandidates();

        /*
	    if (debugging("collision"))
        {
            std::cout<< "    #" << niter++ << ": " << collisionCandidates.size()
                     << " pair of collision tris\n";
            
            std::cout<< "     " << numImpactZones  << " zones of impact\n";
        }
        */

        processCollisionCandidates();
        
        updateImpactZoneVelocity(numImpactZones);

        if (Collisions.empty())
            collision_free = true;
    }
    
    turnOffImpZone();
}

void CollisionSolver3d::updateImpactZoneVelocity(int &nZones)
{
	POINT* pt;
	int numZones = 0;

	unsortHseList(hseList);

    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);

            //skip traversed or isolated pts
            if (sorted(pt) || weight(findSet(pt)) == 1)
                continue;
            else
            {
                updateImpactListVelocity(findSet(pt));
                numZones++;
            }
	    }
    }
    
    nZones = numZones;
}

void CollisionSolver3d::updateImpactZoneVelocityForRG()
{
	POINT* pt;
	unsortHseList(hseList);

    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            pt = (*it)->Point_of_hse(i);
            
            //skip traversed or isolated pts
            if (sorted(pt) || weight(findSet(pt)) == 1)
                continue;
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

void CollisionSolver3d::updateImpactListVelocity(POINT* head)
{
	STATE* sl = nullptr;
	POINT* p = head;

	double m = getPointMass();
    double x_cm[3] = {0.0};
    double v_cm[3] = {0.0};
	
    double L[3] = {0.0};    //angular momentum
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


