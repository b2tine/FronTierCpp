#ifndef COLLID_H
#define COLLID_H

#include "AABB.h"

#if defined(isnan)
#undef isnan
#endif

#define DEBUGGING false

//TODO: Determine ROUND_EPS optimal value; currently it seems
//      to be too small to be meaningful in the context it is used.
const double ROUND_EPS = DBL_EPSILON;
//const double ROUND_EPS = 1.0e-10;
const double EPS = 1.0e-06;
const double DT = 0.001;


class CollisionSolver3d {
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
	
    double setVolumeDiff(double);
    
    void setStrainLimit(double);
	double getStrainLimit();
    void setCompressiveStrainLimit(double);
	double getCompressiveStrainLimit();
	void setStrainRateLimit(double);
	double getStrainRateLimit();
	static bool getGsUpdateStatus();	

	void clearHseList();
    const std::vector<CD_HSE*>& getHseList() const;

    void initializeSystem(Front* front);
	void assembleFromInterface(INTERFACE*);
    void assembleFromSurf(SURFACE* surf);//TODO: Need this for parallel runs?
    void assembleFromCurve(CURVE* curve);//TODO: Need this for parallel runs?
    void assembleFromNode(NODE* node);//TODO: Need this for parallel runs?
	void recordOriginalPosition();	
    void setHseTypeLists();
    void initializeImpactZones();
    void initializeImpactZones(const INTERFACE* intfc);
	void initRigidBodyImpactZones();
	    //void initRigidBodyImpactZones(const INTERFACE* intfc);//TODO: remove?
	
    void resolveCollision();

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

    static void clearCollisionTimes();
    static void setSizeCollisionTimes(unsigned int size);
    static void addCollisionTime(double collsn_dt);
    static double getAverageCollisionTime();

    static int tstep;
    static int getStep() {return tstep;}
    static void setStep(int step) {tstep = step;}

    static std::string outdir;
    static std::string getOutputDirectory() {return outdir;}
    static void setOutputDirectory(std::string dir) {outdir = dir;}
    static void saveFront() {FT_Save(ft);}
    static void drawFront() {FT_Draw(ft);}

private:

    static Front* ft; //static so we can call FT_Save() and FT_Draw() for debugging 
    
    double max_fabric_speed {0.0};
    double prev_max_fabric_speed {0.0};

	std::unique_ptr<AABBTree> abt_proximity {nullptr};
    std::unique_ptr<AABBTree> abt_collision {nullptr};

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

    static bool gs_update;
	static void turnOnGsUpdate();
    static void turnOffGsUpdate();

    double strain_limit {0.1};
    double compressive_strain_limit {0.1};
    double strainrate_limit {0.1};

    static bool s_detImpZone;
	static void turnOnImpZone();
    static void turnOffImpZone();

    int numImpactZones {0};
    int numImpactZonePoints {0};

    std::vector<CD_HSE*> getHseTypeList(CD_HSE_TYPE type);
    
    void limitStrainPosnJac(MotionState mstate);
    void limitStrainPosnGS(MotionState mstate);
    int computeStrainImpulsesPosn(std::vector<CD_HSE*>& list, MotionState mstate);
    void limitStrainRatePosnJac(MotionState mstate);
    void limitStrainRatePosnGS(MotionState mstate);
    int computeStrainRateImpulsesPosn(std::vector<CD_HSE*>& list, MotionState mstate);
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
	void infoImpactZones();
	void debugImpactZones();
	void markImpactZonePoints(POINT* head);
	void updateImpactZoneVelocity();
	void updateImpactZoneVelocityForRG();
	void detectProximityRGB();
	void detectProximity(std::vector<CD_HSE*>& list);
	void detectCollision(std::vector<CD_HSE*>& list);
    void aabbProximity(std::vector<CD_HSE*>& list);
    void aabbCollision(std::vector<CD_HSE*>& list);
	void detectDomainBoundaryCollision();
	void updateFinalForRG();

    void writeCollisionPoints();
    std::vector<POINT*> getCollisionPoints();
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

bool isImpactZonePoint(POINT*);
bool isStaticRigidBody(const POINT*);
bool isStaticRigidBody(const STATE*);
bool isStaticRigidBody(const CD_HSE*);
bool isMovableRigidBody(const POINT*);
bool isMovableRigidBody(const STATE*);
bool isMovableRigidBody(const CD_HSE*);
bool isRigidBody(const POINT*);
bool isRigidBody(const STATE*);
bool isRigidBody(const CD_HSE*);

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
