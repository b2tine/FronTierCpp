#ifndef AABB_H
#define AABB_H

#include <FronTier.h>

#include <fstream>
#include <memory>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <stack>
#include <queue>
#include <set>
#include <map>

using CPoint = std::vector<double>;

enum class MotionState {STATIC, MOVING};

//abstract base class for hypersurface element(HSE)
//can be a point, bond, or triangle
struct CD_HSE
{
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
struct CD_TRI: public CD_HSE
{
    TRI* m_tri;
	
    CD_TRI(TRI* tri, const char* n)
        : m_tri(tri)
    {
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
struct CD_BOND: public CD_HSE
{
	int m_dim;
    BOND* m_bond;
	
    CD_BOND(BOND* bond, int dim, const char* n)
        : m_bond(bond), m_dim(dim)
    {
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
struct CD_POINT: public CD_HSE
{
    POINT* m_point;

    CD_POINT(POINT* point)
        : m_point(point)
    {}

	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const {return 1;}
};


class Node;
class AABBTree;

class AABB {
    friend class Node;
    friend class AABBTree;
    CPoint lowerbound;
    CPoint upperbound;
    // indices will store the index of points on the  
    // corresponding triangle or bond.
    std::vector<long> indices;
    void updateAABBInfo(double);
    //void updateAABBInfo(const std::unordered_map<long, POINT*>&);
    bool contain(const AABB*);
    CD_HSE* hse = nullptr;
    double dt;
    MotionState abType;
    double tol;
public:
    // constructor
    AABB() {}
    AABB(double, CD_HSE*, MotionState);
    AABB(double, CD_HSE*, MotionState, double);
    AABB(const CPoint&, const CPoint&);
    // explicit saying that we need a default version of 
    // copy and move operations 
    AABB(const AABB&) = default;
    AABB& operator=(const AABB&) = default;
    AABB(AABB&&) = default;
    AABB& operator=(AABB&&) = default;
    ~AABB() = default;
    // merge this with anther AABB to get a
    // merged AABB and construct the corresponding AABB tree
    AABB merge(const AABB&) const;
    // get the volume of the AABB
    double volume();
    bool isCollid(const AABB&);
};

class Node {
public:
    friend class AABBTree;
    // AABB stored in node. May store information for branch AABB
    // and may be adjusted for dynamic AABB 
    AABB box;
    // if leaf, point to the corresponding AABB
    // empty for branch
    std::unique_ptr<AABB> data;
    // parent node
    std::weak_ptr<Node> parent;
    // left and right children node
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
    void updateBranch();
    // make this node to be brance from two Node parameter
    void setBranch(std::shared_ptr<Node>, std::shared_ptr<Node>, std::shared_ptr<Node>);
    // judge if this node is a leaf
    bool isLeaf();
    // set an AABB element to be a leaf
    void setLeaf(AABB*);
    bool isCollid(Node*);
    void updateAABB();
    Node* getSibling() const;
    ~Node();
};

class AABBTree {
public:
    std::shared_ptr<Node> root;
    // node needed to be removed and reinsert to the tree
    std::unordered_map<long, POINT*> ump;
    // map from object's indices (2 or 3 points' global indices)
    // to corresponding CD_HSE* in collision library 
    std::map<std::vector<long>, CD_HSE*> vhMap;
    std::unordered_set<Node*> nodeSet;
    std::vector<std::shared_ptr<Node>> nodeArray;
    int count;
    int numLeaf = 0;
    double treeHeight(Node*); 
    double dt;
    bool isProximity;
    bool isCollsn;
    
    // query all collid pairs
    void query(double tol);

    // insert a node into the subtree with parent 
    // as the root
    void insertNode(std::shared_ptr<Node>, std::shared_ptr<Node>&);
    MotionState type;
    double tolerance;
    
    AABBTree(MotionState mstate);
    ~AABBTree();
    void deleteTree();

    // don't want tree to be copied or moved
    AABBTree(const AABBTree&) = delete;
    AABBTree& operator=(const AABBTree&) = delete;
    AABBTree(AABBTree&&) = delete;
    AABBTree& operator=(AABBTree&&) = delete;
    
    // add an AABB element into a tree
    void addAABB(AABB*);
    
    int getCount() { return count; }
    double getVolume() { return root->box.volume(); } 
    void updateTreeStructure();
    void updatePointMap(const std::vector<CD_HSE*>&);
    void setTimeStep(double t) { dt = t; }
    bool getCollsnState() { return isCollsn; }
    void updateAABBTree(const std::vector<CD_HSE*>&);
    MotionState getType() { return type; }

private:

    bool queryProximity(Node* n,double tol);
    bool queryCollision(Node* n,double tol);
};


//dcollid.cpp
bool getProximity(const CD_HSE*,const CD_HSE*,double);
bool getCollision(const CD_HSE*,const CD_HSE*,double);


#endif
