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


enum class MotionState {STATIC, MOVING};

class Node;
class AABBTree;

using CPoint = std::vector<double>;
using NodePair = std::pair<std::shared_ptr<Node>,std::shared_ptr<Node>>;


class AABB {
    
    friend class Node;
    friend class AABBTree;

public:
    
    CPoint lowerbound;
    CPoint upperbound;
    // indices will store the index of points on the  
    // corresponding triangle or bond.
    std::vector<long> indices;
    void updateAABBInfo(double);
    //void updateAABBInfo(const std::unordered_map<long, POINT*>&);
    CD_HSE* hse = nullptr;
    double dt;
    MotionState abType;
    double tol;

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
    const double volume() const;
    const bool overlaps(const AABB&) const;
    const bool contains(const AABB&) const;
};

class Node
{
public:

    friend class AABBTree;
    
    // AABB stored in node. May store information for branch AABB
    // and may be adjusted for dynamic AABB 
    AABB box;

    const double volume() const;
    const bool overlaps(Node*) const;
    const bool overlaps(std::shared_ptr<Node> node) const;
    
    // if leaf, point to the corresponding AABB (empty for branch)
    //std::unique_ptr<AABB> data;
    AABB* data;
    const CD_HSE* const getHSE() const;
    
    const bool hasAdjacentHSE(std::shared_ptr<Node> node) const;

    
    std::weak_ptr<Node> parent;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
    
    Node* getLeftChild() const;
    Node* getRightChild() const;
    Node* getSibling() const;

    void updateBranch();
    
    // make this node to be branch from two Node parameter
    void setBranch(std::shared_ptr<Node>, std::shared_ptr<Node>, std::shared_ptr<Node>);
    
    void setLeaf(AABB*);
    const bool isLeaf() const;

    void updateAABB();

    ~Node();
};

class AABBTree
{

public:

    std::shared_ptr<Node> root;
    // node needed to be removed and reinsert to the tree
    std::unordered_map<long, POINT*> ump;
    // map from object's indices (2 or 3 points' global indices)
    // to corresponding CD_HSE* in collision library 
    std::map<std::vector<long>, CD_HSE*> vhMap;
    std::unordered_set<Node*> nodeSet;
    std::vector<std::shared_ptr<Node>> nodeArray;

    int count {0};
    int numLeaf {0};
    double treeHeight(Node*); 
    double dt;
    
    // query all collid pairs
    void query(double tol);
    std::vector<NodePair> getCandidates();

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
    void updateAABBTree(const std::vector<CD_HSE*>&);
    MotionState getType() { return type; }

private:

    //bool queryProximity(Node* n,double tol);
    //bool queryCollision(Node* n,double tol);
};


//dcollid.cpp
//bool getProximity(const CD_HSE*,const CD_HSE*,double);
//bool getCollision(const CD_HSE*,const CD_HSE*,double);


#endif
