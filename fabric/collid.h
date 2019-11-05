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


class Proximity;
class Collision;


class CollisionSolver3d
{
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

	void recordOriginalPosition();	

	void clearHseList();
	void assembleFromInterface(const INTERFACE*,double dt);
	void createImpZoneForRG(const INTERFACE*);
	void setDomainBoundary(double* L,double *U);
	
    void resolveCollision();

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

    std::vector<NodePair> proximityCandidates;
    std::vector<std::unique_ptr<Proximity>> Proximities;
    
    std::vector<NodePair> collisionCandidates;
    std::vector<std::unique_ptr<Collision>> Collisions;

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
	
	void computeAverageVelocity();
	//void updateAverageVelocity();
    void resetPositionCoordinates();
	
    void aabbProximity();
	void detectProximity();
	void processProximityCandidates();

    void aabbCollision();
	void detectCollision();
	void processCollisionCandidates();
	
    void detectDomainBoundaryCollision();
	void updateFinalForRG();

    void updateFinalPosition();
	void updateFinalVelocity();
	
    void reduceSuperelast();
    bool reduceSuperelastOnce(int&);
	
    void computeImpactZone();
	void updateImpactZoneVelocity(int&);
	void updateImpactZoneVelocityForRG();
	void updateImpactListVelocity(POINT*);
};


//dcollid.cpp
void vtkplotVectorSurface(std::vector<CD_HSE*>&,const char*);

void unsort_surface_point(SURFACE *surf);
void unsortHseList(std::vector<CD_HSE*>&);
void printPointList(POINT**, const int);

bool isStaticRigidBody(const POINT*);
bool isStaticRigidBody(const CD_HSE*);
bool isMovableRigidBody(const POINT*);
bool isMovableRigidBody(const CD_HSE*);
bool isRigidBody(const POINT*);
bool isRigidBody(const CD_HSE*);

void scalarMult(double a,double* v, double* ans); 
void addVec(double* v1, double* v2, double* ans); 
void minusVec(double* v1, double* v2, double* ans); 
void Pts2Vec(const POINT*, const POINT*, double*); 
double distBetweenCoords(double* v1, double* v2);
double myDet3d(double[][3]);

//Proximity.cpp
std::unique_ptr<Proximity> checkProximity(const CD_HSE*,const CD_HSE*,double);
std::unique_ptr<Collision> checkCollision(const CD_HSE*,const CD_HSE*,double);

//Impulse.cpp
void PointToTriProximityImpulse(POINT**,double*,double*,double);
void PointToTriPostCollisionProximityImpulse(POINT**,double*,double*,double);
void EdgeToEdgeCollisionImpulse(POINT**,double*,double,double,double,double);

void EdgeToEdgeProximityImpulse(POINT**,double*,double,double,double);
void EdgeToEdgePostCollisionProximityImpulse(POINT**,double*,double,double,double);
void PointToTriCollisionImpulse(POINT**,double*,double*,double,double);

//ImpactZone.cpp
void createImpZone(POINT**, int num = 4, bool first = NO);
int& weight(POINT* p);
POINT*& root(POINT* p);
POINT*& next_pt(POINT* p);
POINT*& tail(POINT* p);
void makeSet(std::vector<CD_HSE*>& hselist);
POINT* findSet(POINT* p);
void mergePoint(POINT* X, POINT* Y);



#endif
