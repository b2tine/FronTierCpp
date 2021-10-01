#ifndef COLLID_H
#define COLLID_H

#include "AABB.h"

#if defined(isnan)
#undef isnan
#endif

#define DEBUGGING false

//TODO: Consistent use of ROUND_EPS AND MACH_EPS.
//      Currently they are equal, and we can change all to MACH_EPS. 
const double ROUND_EPS = DBL_EPSILON;


struct StrainStats
{
    int n_edges;
    double total_edge_length;
};


struct CollisionTimeStats
{
    double avg_dt;
    double min_dt;
    double max_dt;
};

struct FABRIC_COLLISION_PARAMS
{
    double fabric_eps {1.0e-05};
    double fabric_thickness {0.001};
    double mu_s;
    double k_s;
    double m_s;

    double string_eps {4.0e-05};
    double string_thickness {0.004};
    double mu_l;
    double k_l;
    double m_l;

    double overlap_coefficient {0.5};

    double strain_limit {0.01};
    double compressive_strain_limit {0.0};
    double strainrate_limit {0.05};
    double strain_vel_tol {0.15};

    double coefRestitution {1.0};

    bool collision_off {false};
};


class CollisionSolver3d 
{
private:

    static Front* front; //static so we can call FT_Save() and FT_Draw() for debugging 
    
    static bool collision_off;

public:
	
    int m_dim {3};
	std::vector<CD_HSE*> hseList;
    std::vector<CD_HSE*> fabricTriList;
    std::vector<CD_HSE*> staticRigidTriList;
    std::vector<CD_HSE*> movableRigidTriList;
    std::vector<CD_HSE*> stringBondList;
    std::vector<CD_HSE*> elasticHseList;
	
    std::map<int,std::vector<double>> mrg_com;
	
    int build_count_pre = 1;
    int build_count_col = 1;

	CollisionSolver3d() = default;
        //CollisionSolver3d(Front* fr); //Won't work right now since we need front to be static
    virtual ~CollisionSolver3d();

    CollisionSolver3d(const CollisionSolver3d&) = delete;
    CollisionSolver3d& operator=(const CollisionSolver3d&) = delete;
    CollisionSolver3d(CollisionSolver3d&&) = delete;
    CollisionSolver3d& operator=(CollisionSolver3d&&) = delete;

    static void turnCollision_ON();
    static void turnCollision_OFF();
    static bool collisionEnabled();

	static void setTimeStepSize(double);
	static double getTimeStepSize();
	
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

    static void setOverlapCoefficient(double);
    static double getOverlapCoefficient();
    
    void setStrainLimit(double);
	double getStrainLimit();
    void setCompressiveStrainLimit(double);
	double getCompressiveStrainLimit();
	void setStrainRateLimit(double);
	double getStrainRateLimit();
    void setStrainVelocityTol(double);
	double getStrainVeloctiyTol();
	
    static void setRestitutionCoef(double);
	static double getRestitutionCoef();
	
    static bool getImpZoneStatus();	
    static bool getGsUpdateStatus();	


	void clearHseList();
    const std::vector<CD_HSE*>& getHseList() const;

    void initFront(Front*fr);
    void initializeSystem(std::vector<CD_HSE*>& list);

    void initializeSystem(Front* fr);
	void assembleFromInterface(INTERFACE*);
    void assembleFromSurf(SURFACE* surf); //TODO: Need this for parallel runs?
    void assembleFromCurve(CURVE* curve); //TODO: Need this for parallel runs?
	void recordOriginalPosition();	
    void setHseTypeLists();
    void initializeImpactZones();
	void initRigidBodyImpactZones();
	
    void resolveCollision();
    void resolveCollisionSubstep();

	void setDomainBoundary(double* L,double *U);
	double getDomainBoundary(int dir,int side) {return Boundary[dir][side];}
	
    double setVolumeDiff(double);
    
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

    static void clearCollisionTimes();
    static void setSizeCollisionTimes(unsigned int size);
    static void addCollisionTime(double collsn_dt);
    static double getAverageCollisionTime();
    static double getMaxCollisionTime();
    static double getMinCollisionTime();
    static CollisionTimeStats getCollisionTimeStats();

    static int tstep;
    static int getStep() {return tstep;}
    static void setStep(int step) {tstep = step;}
	static void setFrameTimeStepSize(double);
	static double getFrameTimeStepSize();

    static std::string outdir;
    static std::string getOutputDirectory() {return outdir;}
    static void setOutputDirectory(std::string dir) {outdir = dir;}
    static void saveFront() {FT_Save(front);}
    static void drawFront() {FT_Draw(front);}

private:
    
    double max_fabric_speed {0.0};
    double prev_max_fabric_speed {0.0};

	std::unique_ptr<AABBTree> abt_proximity {nullptr};
    std::unique_ptr<AABBTree> abt_collision {nullptr};

	double Boundary[3][2]; //domain boundary[dir][side]

    double volume;
    double vol_diff {0.0};
    
    static std::vector<double> CollisionTimes;

	static double s_dt;
	static double frame_dt;
	
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

    static double overlap_coefficient;

    static bool gs_update;
	static void turnOnGsUpdate();
    static void turnOffGsUpdate();

    double strain_limit {0.1};
    double compressive_strain_limit {0.01};
    double strainrate_limit {0.05};
    double strain_vel_tol {0.15};

    static bool s_detImpZone;
	static void turnOnImpZone();
    static void turnOffImpZone();

    int numImpactZones {0};
    int numImpactZonePoints {0};

    std::vector<CD_HSE*> getHseTypeList(CD_HSE_TYPE type);
    std::vector<CD_HSE*> shuffleHseList(const std::vector<CD_HSE*>& list) const;
    
    void limitStrainPosnJac(MotionState mstate);
    void limitStrainPosnGS(MotionState mstate);
    StrainStats computeStrainImpulsesPosn(std::vector<CD_HSE*>& list, MotionState mstate);
    
    void limitStrainRatePosnJac(MotionState mstate);
    void limitStrainRatePosnGS(MotionState mstate);
    StrainStats computeStrainRateImpulsesPosn(std::vector<CD_HSE*>& list, MotionState mstate);

    void limitStrainVelJAC();
    void limitStrainVelGS();
    int computeStrainImpulsesVel(std::vector<CD_HSE*>& list);
    
    void applyStrainImpulses(MotionState mstate);

	void computeMaxSpeed();
	void computeAverageVelocity();
    void resetPositionCoordinates();
	void updateFinalPosition();
	void updateFinalVelocity();
    void updateFinalStates();
	void updateAverageVelocity(MotionState mstate);
	void saveAverageVelocity();
	void revertAverageVelocity();

	void computeImpactZoneGS(std::vector<CD_HSE*>& list);
	void computeImpactZoneJac(std::vector<CD_HSE*>& list);
    void connectNearbyImpactZones(std::vector<CD_HSE*>& list);
	void markImpactZonePoints(POINT* head);
	void updateImpactZoneVelocity();
	void updateImpactZoneVelocityForRG();
	void infoImpactZones();
	void debugImpactZones();
	
    void detectProximityRGB();
	void detectProximity(std::vector<CD_HSE*>& list);
	void detectCollision(std::vector<CD_HSE*>& list);
    void aabbProximity(std::vector<CD_HSE*>& list);
    void aabbCollision(std::vector<CD_HSE*>& list);
	void detectDomainBoundaryCollision();
	void updateFinalForRG();

    void writeCollisionPoints();
    std::vector<POINT*> getCollisionPoints();
    //TODO: write these functions
        //void writeProximityPoints();
        //std::vector<POINT*> getProximityPoints();
};

void unsortHseList(std::vector<CD_HSE*>&);

bool BondToBond(const BOND*,const BOND*);
bool TriToBond(const TRI*,const BOND*);
bool TriToTri(const TRI*,const TRI*);
bool MovingBondToBondGS(const BOND*,const BOND*);
bool MovingTriToBondGS(const TRI*,const BOND*);
bool MovingTriToTriGS(const TRI*,const TRI*);
bool MovingBondToBondJac(const BOND*,const BOND*);
bool MovingTriToBondJac(const TRI*,const BOND*);
bool MovingTriToTriJac(const TRI*,const TRI*);

void makeSet(std::vector<CD_HSE*>&);
void createImpZone(POINT*[],int num = 4,bool first = NO);
void createImpactZone(POINT*[],int num);
void createImpactZoneRigidBody(POINT*[],int num);
void updateImpactListVelocity(POINT* head);
void SpreadImpactZoneImpulse(POINT*, double, double*);
void printPointList(POINT**, const int);

void vtkplotVectorSurface(std::vector<CD_HSE*>&,const char*);

POINT* findSet(POINT* p);
POINT*& next_pt(POINT* p);
void mergePoint(POINT* X, POINT* Y);
int& weight(POINT* p);

bool isStaticRigidBody(const POINT*);
bool isStaticRigidBody(const STATE*);
bool isStaticRigidBody(const CD_HSE*);
bool isMovableRigidBody(const POINT*);
bool isMovableRigidBody(const STATE*);
bool isMovableRigidBody(const CD_HSE*);
bool isRigidBody(const POINT*);
bool isRigidBody(const STATE*);
bool isRigidBody(const CD_HSE*);
bool isRegisteredPoint(const POINT*);
bool isRegisteredPoint(const STATE*);
bool isConstrainedPoint(const POINT*);
bool isConstrainedPoint(const STATE*);
bool isImpactZonePoint(POINT*);

void initSurfaceState(SURFACE*,const double*);
void initCurveState(CURVE*,const double*);
void initTestModule(Front&, char*);

//TODO: Move into new file and make available globally
void Pts2Vec(const POINT*, const POINT*, double*);
void scalarMult(double a,double* v, double* ans);
void addVec(double* v1, double* v2, double* ans);
void minusVec(double* v1, double* v2, double* ans);
double myDet3d(double[][3]);
double distBetweenCoords(double* x1, double* x2);

#endif
