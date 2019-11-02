#ifndef COLLID_H
#define COLLID_H

#include "AABB.h"

#if defined(isnan)
#undef isnan
#endif

#define DEBUGGING false

//TODO: How where these determined?
const double ROUND_EPS = 1.0e-10;
const double EPS = 1.0e-6;
const double DT = 0.001;


class CollisionSolver3d {
public:
	
    int m_dim {3};
	std::vector<CD_HSE*> hseList;
	std::map< int, std::vector<double> > mrg_com;
	
    int build_count_pre = 1;
    int build_count_col = 1;

	CollisionSolver3d() = default;
    virtual ~CollisionSolver3d();

    CollisionSolver3d(const CollisionSolver3d&) = delete;
    CollisionSolver3d& operator=(const CollisionSolver3d&) = delete;
    CollisionSolver3d(CollisionSolver3d&&) = delete;
    CollisionSolver3d& operator=(CollisionSolver3d&&) = delete;

	static void setRoundingTolerance(double);
	static double getRoundingTolerance();
	static void setFabricThickness(double);
	static double getFabricThickness();
	static void setTimeStepSize(double);
	static double getTimeStepSize();
	static void setSpringConstant(double);
	static double getSpringConstant();
	static void setFrictionConstant(double);
	static double getFrictionConstant();
	static void setPointMass(double);
	static double getPointMass();
	static void setRestitutionCoef(double);
	static double getRestitutionCoef();
	static bool getImpZoneStatus();	

    double setVolumeDiff(double);

	void clearHseList();
	void assembleFromInterface(const INTERFACE*,double dt);
	void createImpZoneForRG(const INTERFACE*);
	
    void resolveCollision();
	void recordOriginalPosition();	
	void setDomainBoundary(double* L,double *U);

	double getDomainBoundary(int dir,int side) {return Boundary[dir][side];}

    POINT **gpoints;
    TRI **gtris;

	//for debugging
	static void printDebugVariable();
	static int moving_edg_to_edg;
	static int moving_pt_to_tri;
	static int is_coplanar;
	static int edg_to_edg;
	static int pt_to_tri;

    TRI *res_tris[100];
    int num_res_tris;

private:

	std::unique_ptr<AABBTree> abt_proximity {nullptr};
    std::unique_ptr<AABBTree> abt_collision {nullptr};

    //TODO: put this in AABBTree class
    double vol_diff {0.0};
    double collision_vol;
    double proximity_vol;

    std::vector<NodePair> candidates;

	static double s_eps;
	static double s_thickness;
	static double s_dt;
	static double s_m;
	static double s_k;
	static double s_mu;
	static double s_cr;
    static bool s_detImpZone;

	double Boundary[3][2]; //domain boundary[dir][side]

	static void turnOffImpZone();
	static void turnOnImpZone();
	
    bool reduceSuperelastOnce(int&);
	void computeAverageVelocity();
    void resetPositionCoordinates();
	void updateFinalPosition();
	void reduceSuperelast();
	void updateFinalVelocity();
	void updateAverageVelocity();
	void computeImpactZone();
	void updateImpactZoneVelocity(int&);
	void updateImpactZoneVelocityForRG();
    void aabbProximity();
    void aabbCollision();
	void detectProximity();
	void detectCollision();
	void processProximityCandidates();
	void processCollisionCandidates();
	void detectDomainBoundaryCollision();
	void updateFinalForRG();
	void updateImpactListVelocity(POINT*);
};

void TriToTri(const TRI*,const TRI*,double);
void TriToBond(const TRI*,const BOND*,double);
void BondToBond(const BOND*,const BOND*,double);

bool MovingTriToTri(const TRI*,const TRI*,double);
bool MovingTriToBond(const TRI*,const BOND*,double);
bool MovingBondToBond(const BOND*,const BOND*,double);

void checkProximity(const CD_HSE*,const CD_HSE*,double);
void checkCollision(const CD_HSE*,const CD_HSE*,double);

void EdgeToEdgeProximity();
void EdgeToEdgeProximity();

void EdgeToEdgeProximityImpulse(POINT**,double*,double,double,double);
void PointToTriProximityImpulse(POINT**,double*,double*,double);

class Proximity
{
    public:

        POINT** pts {nullptr};
        double nor[3] {0.0};
        double dist {HUGE};
        
        virtual void applyImpulse() = 0;
        
        virtual ~Proximity() = default;
};

class EdgeEdgeProximity : public Proximity
{
    public:

        double a {-1.0};
        double b {-1.0};
        
        EdgeEdgeProximity(POINT** Pts, double* Nor,
                double A, double B, double Dist)
            : pts{Pts}, a{A}, b{B}, dist{Dist}
        {
            for (int i = 0; i < 3; ++i)
                nor[i] = Nor[i];
        }
        
        virtual void applyImpulse() override
        {
            EdgeToEdgeProximityImpulse(pts,nor,a,b,dist);
        }
};

class PointTriProximity : public Proximity
{
    public:

        double w[3] {-1.0};

        PointTriProximity(POINT** Pts,
                double* Nor, double* W, double Dist)
            : pts{Pts}, dist{Dist}
        {
            for (int i = 0; i < 3; ++i)
            {
                w[i] = W[i];
                nor[i] = Nor[i];
            }
        }

        void applyImpulse() override
        {
            PointToTriProximityImpulse(pts,nor,w,dist);
        }
};

class EdgeEdgeCollision : public EdgeEdgeProximity
{
    public:

        double dt {-1.0};

        EdgeEdgeCollision(POINT** Pts, double* Nor,
                double A, double B, double Dist, double Dt)
            : EdgeEdgeProximity(Pts,Nor,A,B,Dist), dt{DT}
        {}

        //TODO: applyImpulse()
};

class PointTriCollision : public PointTriProximity
{
    public:

        double dt {-1.0};

        PointTriCollision(POINT** Pts, double* Nor,
                double* W, double Dist, double Dt)
            : PointTriProximity(Pts,Nor,W,Dist), dt{Dt}
        {}

        //TODO: applyImpulse()
};


void initSurfaceState(SURFACE*,const double*);
void initCurveState(CURVE*,const double*);
void initTestModule(Front&, char*);
void Pts2Vec(const POINT*, const POINT*, double*); 
void scalarMult(double a,double* v, double* ans); 
void addVec(double* v1, double* v2, double* ans); 
void minusVec(double* v1, double* v2, double* ans); 
double myDet3d(double[][3]);
double distBetweenCoords(double* v1, double* v2);
extern void printPointList(POINT**, const int);
extern void createImpZone(POINT*[],int num = 4,bool first = NO);
extern void makeSet(std::vector<CD_HSE*>&);
void unsortHseList(std::vector<CD_HSE*>&);
POINT*& next_pt(POINT*);
int& weight(POINT*);
bool isStaticRigidBody(const POINT*);
bool isStaticRigidBody(const CD_HSE*);
bool isMovableRigidBody(const POINT*);
bool isMovableRigidBody(const CD_HSE*);
bool isRigidBody(const POINT*);
bool isRigidBody(const CD_HSE*);
extern void SpreadImpactZoneImpulse(POINT*, double, double*);

void vtkplotVectorSurface(std::vector<CD_HSE*>&,const char*);



#endif
