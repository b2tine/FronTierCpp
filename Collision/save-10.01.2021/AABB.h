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
    
    // indices store the index of points on the  
    // corresponding triangle or bond.
    std::vector<long> indices;
    
    void updateAABBInfo(double);
    //void updateAABBInfo(const std::unordered_map<long, POINT*>&);
    
    bool contain(const AABB*);
    CD_HSE* hse {nullptr};
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

class Node
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
    
    Node* parent {nullptr}; //std::weak_ptr<Node> parent;
    Node* left {nullptr};
    Node* right {nullptr};

    void updateBranch();
    
    // make this node to be branch from two Node parameter
    void setBranch(Node* n1, Node* n2);
    
    // check if this node is a leaf
    bool isLeaf() const;
    
    // set an AABB element to be a leaf
    void setLeaf(AABB*);
    bool isCollid(Node* node) const;
    void updateAABB();
    Node* getSibling() const;

    ~Node();
};

class AABBTree {
public:
    Node* root {nullptr};
    
    // node needed to be removed and reinsert to the tree
    std::unordered_map<long, POINT*> ump;
    
    // map from object's indices (2 or 3 points' global indices)
    // to corresponding CD_HSE* in collision library 
    std::map<std::vector<long>, CD_HSE*> vhMap;
    std::unordered_set<Node*> nodeSet;

    //only used for rebuilding tree in call to updateTreeStructure()
    std::vector<Node*> nodeArray;

    int count {0};
    std::vector<std::pair<CD_HSE*,CD_HSE*>> interference_pairs;
    std::vector<std::pair<CD_HSE*,CD_HSE*>> getInterferencePairs() const;

    int numLeaf {0};
    double treeHeight(Node* node); 
    double dt;
    bool isProximity {false};
    bool isCollsn {false};
    
    // query all collid pairs
    void query();

    // insert a node into the subtree with parent as the root
    Node* insertNode(Node* n, Node* parentNode);
    
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

    void turn_on_GS_update() {gauss_seidel = true;}
    void turn_off_GS_update() {gauss_seidel = false;}
    bool get_GS_update_status() const {return gauss_seidel;}

private:

    bool gauss_seidel {false};

    bool queryProximity(Node* n);
    bool queryCollision(Node* n);
    bool getCollision(const CD_HSE* a, const CD_HSE* b);
};


//dcollid.cpp
bool getProximity(const CD_HSE* a, const CD_HSE* b);
bool getCollisionGS(const CD_HSE* a, const CD_HSE* b);
bool getCollisionJac(const CD_HSE* a, const CD_HSE* b);


#endif
