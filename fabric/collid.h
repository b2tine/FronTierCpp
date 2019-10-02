#ifndef COLLID_H
#define COLLID_H

#include <FronTier.h>
#include "state.h"

#include <functional>
#include <map>
#include <fstream>
#include <memory>

#if defined(isnan)
#undef isnan
#endif

//TODO: How where these determined?
#define DEBUGGING false
const double ROUND_EPS = 1e-10;
const double EPS = 1e-6;
const double DT = 0.001;

/*
user-defined state should include the following
struct UF{
	POINT* next_pt;
	POINT* root;
	POINT* tail;
	int num_pts;
};

struct STATE{
	double vel[3];
	double collsnImpulse[3];
	double friction[3];
	double avgVel[3];
	double x_old[3];
	int    collsn_num;
	bool   has_collsn;
	bool   is_fixed;
	UF     impZone;
};*/

//abstract base class for hypersurface element(HSE)
//can be a point or a bond or a triangle
class CD_HSE{
public:
        std::string name;
	virtual double max_static_coord(int) = 0;
	virtual double min_static_coord(int) = 0;
	virtual double max_moving_coord(int,double) = 0;
	virtual double min_moving_coord(int,double) = 0;
	virtual POINT* Point_of_hse(int) const  = 0;
	virtual int num_pts() const= 0;
	virtual ~CD_HSE(){};
};

//wrap class for triangle
class CD_TRI: public CD_HSE{
public:
	TRI* m_tri;
	CD_TRI(TRI* tri, const char* n):m_tri(tri){
            name = n;
        }
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts() const {return 3;}
};

//wrap class for bond
class CD_BOND: public CD_HSE{
public:
	BOND* m_bond;
	int m_dim;
	CD_BOND(BOND* bond, int dim, const char* n):m_bond(bond), m_dim(dim){
            name = n;
        }
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const{return 2;}
};

//wrap class for point
class CD_POINT: public CD_HSE{
public:
	POINT* m_point;
	CD_POINT(POINT* point):m_point(point){}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const {return 1;}
};

//Forward declaration of AABBTree
class AABBTree;

//abstract base class for collision detection and handling
class CollisionSolver {
private:
	std::unique_ptr<AABBTree> abt_proximity;
    std::unique_ptr<AABBTree> abt_collision;
    double volume;
    double vol_diff {0.0};
	static double s_eps;
	static double s_thickness;
	static double s_dt;
	static double s_m;
	static double s_k;
	static double s_lambda;
	static double s_cr;
	bool has_collision;
	double Boundary[3][2]; //domain boundary[dir][side]
	static void turnOffImpZone();
	static void turnOnImpZone();
	bool reduceSuperelastOnce(int&);
	void computeAverageVelocity();
	void updateFinalPosition();
	void reduceSuperelast();
	void updateFinalVelocity();
	void updateAverageVelocity();
	void computeImpactZone();
	void updateImpactZoneVelocity(int&);
	void updateImpactZoneVelocityForRG();
	void detectProximity();
	void detectCollision();
    void aabbProximity();
    void aabbCollision();
	void detectDomainBoundaryCollision();
	void updateFinalForRG();
	void setHasCollision(bool judge) {has_collision = judge;}

	virtual void updateImpactListVelocity(POINT*) = 0;
	virtual bool BondToBond(const BOND*,const BOND*,double) = 0;
	virtual bool TriToTri(const TRI*,const TRI*,double) = 0;
	virtual bool TriToBond(const TRI*,const BOND*,double)=0;
	virtual bool MovingBondToBond(const BOND*,const BOND*,double) = 0;
	virtual bool MovingTriToTri(const TRI*,const TRI*,double) = 0;
	virtual bool MovingTriToBond(const TRI*,const BOND*,double)=0;
protected:
	int m_dim;
	std::vector<CD_HSE*> hseList;
	std::map< int, std::vector<double> > mrg_com;
	static bool s_detImpZone;
	void clearHseList();
public:
    int build_count_pre = 1;
    int build_count_col = 1;
    enum {STATIC, MOVING};

	CollisionSolver(int);
	CollisionSolver();
    CollisionSolver(CollisionSolver&&);
    CollisionSolver& operator=(CollisionSolver&&);
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
	virtual ~CollisionSolver(); //virtual destructor
	//pure virtual functions
	virtual void assembleFromInterface(const INTERFACE*,double dt) = 0;
	virtual void createImpZoneForRG(const INTERFACE*) = 0;
	bool getProximity(const CD_HSE*,const CD_HSE*);	
	bool getCollision(const CD_HSE*,const CD_HSE*);
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
};

//derived 3D-class for collision detection and handling
class CollisionSolver3d : public CollisionSolver {
private:
	void updateImpactListVelocity(POINT*);
	bool BondToBond(const BOND*,const BOND*,double);
	bool TriToTri(const TRI*,const TRI*,double);
	bool TriToBond(const TRI*,const BOND*,double);
	bool MovingBondToBond(const BOND*,const BOND*,double);
	bool MovingTriToTri(const TRI*,const TRI*,double);
	bool MovingTriToBond(const TRI*,const BOND*,double);
public:
	CollisionSolver3d():CollisionSolver(3){}
	void assembleFromInterface(const INTERFACE*,double dt);
	void createImpZoneForRG(const INTERFACE*);
};


void initSurfaceState(SURFACE*,const double*);
void initCurveState(CURVE*,const double*);
void initTestModule(Front&, char*);
void Pts2Vec(const POINT*, const POINT*, double*); 
void scalarMult(double a,double* v, double* ans); 
void addVec(double* v1, double* v2, double* ans); 
void minusVec(double* v1, double* v2, double* ans); 
bool LeftTurn(std::vector<double>&a,
        std::vector<double>&b, std::vector<double>&c);
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
