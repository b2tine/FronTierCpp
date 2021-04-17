#ifndef AABB_H
#define AABB_H

#include "CD_HSE.h"

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
//using collision_pair = std::pair<CD_HSE*,CD_HSE*>;


enum class MotionState
{
    STATIC,
    MOVING,
    POSTCOLLISION
};


class Node;
class AABBTree;

class AABB
{
    //Node and AABBTree classes may access private members of AABB
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
    AABB(double, CD_HSE*);
    AABB(double, CD_HSE*, double);
    AABB(const CPoint&, const CPoint&);
    
    AABB() = default;
    AABB(const AABB&) = default;
    AABB& operator=(const AABB&) = default;
    AABB(AABB&&) = default;
    AABB& operator=(AABB&&) = default;
    ~AABB() = default;
    // merge this with anther AABB to get a
    // merged AABB and construct the corresponding AABB tree
    AABB merge(const AABB&) const;
    // get the volume of the AABB
    double volume() const;
    bool isCollid(const AABB&) const;
};

class Node : public std::enable_shared_from_this<Node>
{
public:
     
    //Allows the AABBTree class to access private members of Node
    friend class AABBTree; //TODO: which it has none ...

    // AABB stored in node. May store information for branch AABB
    // and may be adjusted for dynamic AABB 
    AABB box;

    // if leaf, point to the corresponding AABB
    // empty for branch
    std::unique_ptr<AABB> data;
    
    std::weak_ptr<Node> parent;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;

    void updateBranch();
    
    // make this node to be branch from two Node parameter
    void setBranch(std::shared_ptr<Node> n1,
                   std::shared_ptr<Node> n2,
                   std::shared_ptr<Node> parent);
    
    // check if this node is a leaf
    bool isLeaf() const;
    
    // set an AABB element to be a leaf
    void setLeaf(AABB*);
    bool isCollid(std::shared_ptr<Node> node);
    void updateAABB();
    std::shared_ptr<Node> getSibling() const;
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
    std::unordered_set<std::shared_ptr<Node>> nodeSet;
        //std::unordered_set<Node*> nodeSet;
    std::vector<std::shared_ptr<Node>> nodeArray;

    int count {0};
    std::vector<std::pair<CD_HSE*,CD_HSE*>> interference_pairs;
    std::vector<std::pair<CD_HSE*,CD_HSE*>> getInterferencePairs() const;

    int numLeaf {0};
    double treeHeight(std::shared_ptr<Node> node); 
    double dt;
    bool isProximity {false};
    bool isCollsn {false};
    
    // query all collid pairs
    void query();

    // insert a node into the subtree with parent 
    // as the root
    void insertNode(std::shared_ptr<Node>, std::shared_ptr<Node>&);
    MotionState type;
    double tolerance;
    
    explicit AABBTree(MotionState mstate);
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
    MotionState getType() {return type;}

    bool turn_on_GS_update() {gauss_seidel = true;}
    bool turn_off_GS_update() {gauss_seidel = false;}
    bool get_GS_update_status() const {return gauss_seidel;}

private:

    bool gauss_seidel {false};

    bool queryProximity(std::shared_ptr<Node> n);
    bool queryCollision(std::shared_ptr<Node> n);
    bool getCollision(const CD_HSE* a, const CD_HSE* b);
};


//dcollid.cpp
bool getProximity(const CD_HSE* a, const CD_HSE* b);
bool getCollisionGS(const CD_HSE* a, const CD_HSE* b);
bool getCollisionJac(const CD_HSE* a, const CD_HSE* b);


#endif
