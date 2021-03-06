#ifndef COLLID_H
#define COLLID_H

#include "AABB.h"

#include <numeric>

#if defined(isnan)
#undef isnan
#endif

#define DEBUGGING false

const double ROUND_EPS = DBL_EPSILON;
//const double ROUND_EPS = 1.0e-10;
const double EPS = 1.0e-06;
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

	static void setTimeStepSize(double);
	static double getTimeStepSize();
	static void setRestitutionCoef(double);
	static double getRestitutionCoef();
	static bool getImpZoneStatus();	
	
    static void setFabricThickness(double);
	static double getFabricThickness();
	static void setFabricSpringConstant(double);
	static double getFabricSpringConstant();
	static void setFabricFrictionConstant(double);
	static double getFabricFrictionConstant();
	static void setFabricPointMass(double);
	static double getFabricPointMass();
	static void setFabricRoundingTolerance(double);
	static double getFabricRoundingTolerance();

	static void setStringThickness(double);
	static double getStringThickness();
	static void setStringSpringConstant(double);
	static double getStringSpringConstant();
	static void setStringFrictionConstant(double);
	static double getStringFrictionConstant();
	static void setStringPointMass(double);
	static double getStringPointMass();
	static void setStringRoundingTolerance(double);
	static double getStringRoundingTolerance();

    
    static void clearCollisionTimes()
    {
        CollisionTimes.clear();
    }

    static void setSizeCollisionTimes(unsigned int size)
    {
        CollisionTimes.reserve(size);
    }

    static void addCollisionTime(double collsn_dt)
    {
        CollisionTimes.push_back(collsn_dt);
    }

    static double getAverageCollisionTime()
    {
        double avg_dt =
            std::accumulate(CollisionTimes.begin(),CollisionTimes.end(),0.0);
        avg_dt /= CollisionTimes.size();
        return avg_dt;
    }

	
    void setStrainLimit(double);
	//double getStrainLimit();
	void setStrainRateLimit(double);
	//double getStrainRateLimit();


    double setVolumeDiff(double);

	void clearHseList();
	void assembleFromInterface(const INTERFACE*,double dt);
	void createImpZoneForRG(const INTERFACE*);
	
    void resolveCollision();
	void recordOriginalPosition();	
	void setDomainBoundary(double* L,double *U);

	double getDomainBoundary(int dir,int side) {return Boundary[dir][side];}
	bool hasCollision() {return has_collision;}

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

    bool has_collision;
	double Boundary[3][2]; //domain boundary[dir][side]

    double volume;
    double vol_diff {0.0};
    
    static std::vector<double> CollisionTimes;

	static double s_dt;
	static double s_cr;

	static double s_thickness; //fabric thickness
	static double s_eps;
	static double s_m;
	static double s_k;
	static double s_mu;
    
	static double l_thickness; //string thickness
	static double l_eps;
	static double l_m;
	static double l_k;
	static double l_mu;

    double strain_limit {0.1};
    double strainrate_limit {0.1};

    static bool s_detImpZone;

    int numImpactZones {0};
    int numImpactZonePoints {0};
	
    static void turnOffImpZone();
	static void turnOnImpZone();

    int numStrainRateEdges {0};
    int numStrainEdges {0};

    void limitStrainRate();
    void modifyStrainRate();
    void limitStrain();
    void modifyStrain();

	void computeAverageVelocity();
    void resetPositionCoordinates();
	void updateFinalPosition();
	void updateFinalVelocity();
	void updateAverageVelocity();
	void updateExternalImpulse();
	void computeImpactZone();
	void infoImpactZones();
	void markImpactZonePoints(POINT* head);
	void updateImpactZoneVelocity();
	void updateImpactZoneVelocityForRG();
	void detectProximity();
	void detectProximityEndStep();
	void detectCollision();
    void aabbProximity();
    void aabbCollision();
	void detectDomainBoundaryCollision();
	void updateFinalForRG();
	void setHasCollision(bool judge) {has_collision = judge;}
};


bool BondToBond(const BOND*,const BOND*);
bool TriToBond(const TRI*,const BOND*);
bool TriToTri(const TRI*,const TRI*);
bool MovingBondToBond(const BOND*,const BOND*);
bool MovingTriToBond(const TRI*,const BOND*);
bool MovingTriToTri(const TRI*,const TRI*);

void updateImpactListVelocity(POINT* head);

void initSurfaceState(SURFACE*,const double*);
void initCurveState(CURVE*,const double*);
void initTestModule(Front&, char*);
void Pts2Vec(const POINT*, const POINT*, double*); 
void scalarMult(double a,double* v, double* ans); 
void addVec(double* v1, double* v2, double* ans); 
void minusVec(double* v1, double* v2, double* ans); 
double myDet3d(double[][3]);
double distBetweenCoords(double* x1, double* x2);
extern void printPointList(POINT**, const int);
extern void createImpZone(POINT*[],int num = 4,bool first = NO);
extern void makeSet(std::vector<CD_HSE*>&);
void unsortHseList(std::vector<CD_HSE*>&);
void unsort_surface_point(SURFACE *surf);

POINT* findSet(POINT* p);
POINT*& next_pt(POINT* p);
void mergePoint(POINT* X, POINT* Y);
int& weight(POINT* p);

bool isStaticRigidBody(const POINT*);
bool isStaticRigidBody(const CD_HSE*);
bool isMovableRigidBody(const POINT*);
bool isMovableRigidBody(const CD_HSE*);
bool isRigidBody(const POINT*);
bool isRigidBody(const CD_HSE*);

extern void SpreadImpactZoneImpulse(POINT*, double, double*);

void vtkplotVectorSurface(std::vector<CD_HSE*>&,const char*);

#endif
