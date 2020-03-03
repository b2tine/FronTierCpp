#ifndef AABB_NODE_H
#define AABB_NODE_H

#include "newAABB.h"
#include "CD_HSE.h"

#include <utility>
#include <stack>


class AABB_Node
{
    public:

        //default destructor, standard and move ctors
        virtual ~AABB_Node() = default;
        AABB_Node() = default;
        AABB_Node(AABB_Node&&) = default;
        AABB_Node& operator=(AABB_Node&&) = default;

        //disable copying for now
        AABB_Node(const AABB_Node&) = delete;
        AABB_Node& operator=(const AABB_Node&) = delete;

        virtual const bool isLeaf() const noexcept = 0;

        void setParent(AABB_Node* const node) noexcept;
        const AABB& getAABB() const noexcept;




    /* 
    // if leaf, point to the corresponding AABB
    // empty for branch
    
    std::unique_ptr<AABB> data;
    

    // set an AABB element to be a leaf
    
    void setLeaf(AABB*);
    bool isCollid(Node*);
    void updateAABB();
    Node* getSibling() const;
    ~Node();

    void updateBranch();
    // make this node to be branch from two Node parameter
    
    void setBranch(std::shared_ptr<Node>, std::shared_ptr<Node>, std::shared_ptr<Node>);
    
    */


    //private:
    protected:

        AABB bv;
        AABB_Node* parent {nullptr};
        
};

class InternalNode : public AABB_Node
{
    public:

        InternalNode(AABB_Node* const lc, AABB_Node* const rc);

        InternalNode(InternalNode&&) = default;
        InternalNode& operator=(InternalNode&&) = default;

        InternalNode() = delete;
        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;

        const bool isLeaf() const noexcept;


    private:

        AABB_Node* left {nullptr};
        AABB_Node* right {nullptr};

};

class LeafNode : public AABB_Node
{
    public:

        LeafNode(CD_HSE* const hse, double pad);

        const bool isLeaf() const noexcept;

    private:

        double pad;
        CD_HSE* hse {nullptr};

};




/*
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

    int count {0};
    int numLeaf {0};
    double treeHeight(Node*); 
    double dt;
    bool isProximity;
    bool isCollsn;
    
    // query all collid pairs
    void query();

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

    bool queryProximity(Node* n);
    bool queryCollision(Node* n);
};
*/


/*
//dcollid.cpp
bool getProximity(const CD_HSE*,const CD_HSE*);
bool getCollision(const CD_HSE*,const CD_HSE*);
*/

#endif
