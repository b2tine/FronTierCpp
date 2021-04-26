#include "AABB.h"

//for proximity detection
AABB::AABB(double t, CD_HSE* h)
    : tol{t}, hse{h}, abType{MotionState::STATIC},
    lowerbound(3), upperbound(3)
{
    for (int i = 0; i < 3; i++)
    {
         lowerbound[i] = h->min_static_coord(i) - tol;
         upperbound[i] = h->max_static_coord(i) + tol;
    }
    for (int i = 0; i < h->num_pts(); i++)
         indices.push_back(h->Point_of_hse(i)->global_index);
}

//for collision detection
AABB::AABB(double t, CD_HSE* h, double Dt)
    : tol{t}, hse{h}, abType{MotionState::MOVING},
    dt{Dt}, lowerbound(3), upperbound(3)
{   
    for (int i = 0; i < 3; i++)
    {
         lowerbound[i] = h->min_moving_coord(i, dt) - tol;
         upperbound[i] = h->max_moving_coord(i, dt) + tol;
    }
    for (int i = 0; i < h->num_pts(); i++) 
         indices.push_back(h->Point_of_hse(i)->global_index);
}

AABB::AABB(const CPoint& pl, const CPoint& pu)
    : lowerbound(pl), upperbound(pu)
{}

AABB AABB::merge(const AABB& ab) const
{
    CPoint pl(3), pu(3);
    for (int i = 0; i < 3; i++)
    {
         pl[i] = std::min(lowerbound[i], ab.lowerbound[i]);
         pu[i] = std::max(upperbound[i], ab.upperbound[i]);
    }
    return AABB(pl, pu);
} 

double AABB::volume() const
{
    double val = 1.0;
    for (int i = 0; i < 3; ++i)
        val *= upperbound[i] - lowerbound[i];
    return val;
}

//This is the intersection test for AABB's.
bool AABB::isCollid(const AABB& ab) const
{
    for (int i = 0; i < 3; ++i)
    {
        if (ab.upperbound[i] < lowerbound[i]) return false;
        if (ab.lowerbound[i] > upperbound[i]) return false;
    }
    return true;
}

void AABB::updateAABBInfo(double dt)
{
    if (abType == MotionState::STATIC)
    {
        for (int i = 0; i < 3; i++)
        {
             lowerbound[i] = hse->min_static_coord(i) - tol;
             upperbound[i] = hse->max_static_coord(i) + tol;
        }
    }
    else
    {
        for (int i = 0; i < 3; i++)
        {
             lowerbound[i] = hse->min_moving_coord(i, dt) - tol;
             upperbound[i] = hse->max_moving_coord(i, dt) + tol;
        }
    }
}

bool AABB::contain(const AABB* ab)
{
    return lowerbound[0] <= ab->lowerbound[0] && lowerbound[1] <= ab->lowerbound[1] &&
        lowerbound[2] <= ab->lowerbound[2] && upperbound[0] >= ab->upperbound[0] &&
        upperbound[1] >= ab->upperbound[1] && upperbound[2] >= ab->upperbound[2];
}

void Node::setBranch(Node* n1, Node* n2)
{
    this->left = n1;
    this->right = n2;
    n1->parent = this;
    n2->parent = this;
}

bool Node::isLeaf() const
{
    return left == nullptr && right == nullptr;
}

void Node::setLeaf(AABB* ab)
{
    data.reset(ab);
}

void Node::updateBranch()
{
    if (this->isLeaf()) return;
    
    for (int i = 0; i < 3; i++)
    {
         box.lowerbound[i] = std::min(left->box.lowerbound[i], right->box.lowerbound[i]);
         box.upperbound[i] = std::max(left->box.upperbound[i], right->box.upperbound[i]);
    }
}

void Node::updateAABB()
{
    if (this->isLeaf())
    {
        box.lowerbound = data->lowerbound;
        box.upperbound = data->upperbound;
    }
    else
    {
        // branch node has no AABB yet
        if (box.lowerbound.size() == 0)
            box = left->box.merge(right->box);
        else 
            updateBranch();
    }
}

bool Node::isCollid(Node* node) const
{
    return box.isCollid(node->box);
}

Node* Node::getSibling() const
{
    if (!parent) return nullptr;
    return (this == parent->left) ? parent->right : parent->left;
}

Node::~Node()
{
    data.reset();
}

AABBTree::AABBTree(MotionState mstate)
    : type{mstate}
{}

AABBTree::~AABBTree()
{
    deleteTree();
}

//Iterative Breadth First (Level Order) Traversal
void AABBTree::deleteTree()
{
    std::queue<Node*> q; //std::queue<std::shared_ptr<Node>> q;
    q.push(this->root);

    while (!q.empty())
    {
        auto node = q.front();
        q.pop();

        if (node->left != nullptr)
            q.push(node->left);
        
        if (node->right != nullptr)
            q.push(node->right);

        delete node;
    }
}

void AABBTree::addAABB(AABB* ab)
{
    if (root)
    {
        Node* node = new Node;
        node->setLeaf(ab);
        node->updateAABB();
        root = insertNode(node,root);
        nodeArray.push_back(node);
        numLeaf++;
    }
    else
    {
        root = new Node;
        root->setLeaf(ab);
        root->updateAABB();
        nodeArray.push_back(root);
        numLeaf++;
    }
}

// reorganize the tree structure
void AABBTree::updateTreeStructure()
{
    root = nullptr;
    for (auto node : nodeArray)
    {
         if (root) 
             root = insertNode(node,root);
         else 
             root = node;
    }
}

Node* AABBTree::insertNode(Node* n, Node* parentNode)
{
    Node* p = parentNode;
    
    // if parent is a leaf node, then create a branch
    // with n and parent to be two children
    if (p->isLeaf())
    {
        Node* newParentNode = new Node;
        
        Node* gp = p->parent;
        newParentNode->parent = gp;

        if (gp)
        {
            if (gp->left == parentNode)
                gp->left = newParentNode;
            else
                gp->right = newParentNode;
        }

        n->parent = newParentNode;
        p->parent = newParentNode;
        newParentNode->left = n;
        newParentNode->right = p;
            
        newParentNode->updateAABB();
        return newParentNode;
    }
    else
    {
        // we have to decide which subtree to insert to
        // the rule is insert to the subtree with smaller volume 
        AABB& abl = p->left->box;
        AABB& abr = p->right->box;
        // get volume after inserting current node to 
        // left or right subtree
            //double vdiff1 = (abl.merge(n->box).volume()-abl.volume())/abl.volume();
            //double vdiff2 = (abr.merge(n->box).volume()-abr.volume())/abr.volume();
            //double vdiff1 = abl.merge(n->box).volume()-abl.volume();
            //double vdiff2 = abr.merge(n->box).volume()-abr.volume();
        double vdiff1 = abl.merge(n->box).volume();
        double vdiff2 = abr.merge(n->box).volume();
            // int vdiff1 = treeHeight(p->left);
            // int vdiff2 = treeHeight(p->right);
        
        // insert to left subtree
        if (vdiff1 < vdiff2)
            p->left = insertNode(n, p->left);
        else
            p->right = insertNode(n, p->right);
    
        parentNode->updateAABB();
        return parentNode;
    }
}

void AABBTree::updatePointMap(const std::vector<CD_HSE*>& hseList)
{
    vhMap.clear();
    nodeSet.clear();
    count = 0; 

    for (auto it : hseList)
    {
         std::vector<long> ids;
         for (int i = 0; i < it->num_pts(); i++) 
              ids.push_back(it->Point_of_hse(i)->global_index);
         vhMap.insert({ids, it});
    }
}

void AABBTree::updateAABBTree(const std::vector<CD_HSE*>& hseList)
{
    updatePointMap(hseList);

    if (!root) return;

    std::stack<Node*> sn;
    Node* cur = root;
    
    //iterative postorder traversal of tree
    do {
        while (cur)
        {
            if (cur->right)
                sn.push(cur->right);
            sn.push(cur);
            cur = cur->left;
        }

        cur = sn.top();
        sn.pop();

        if (!sn.empty() && cur->right && cur->right == sn.top())
        {
            sn.pop();
            sn.push(cur);
            cur = cur->right;
        }
        else
        {
            if (cur->isLeaf())
            {
                cur->data->hse = vhMap[cur->data->indices];
                cur->data->updateAABBInfo(dt);
                cur->updateAABB();    
            }

            if (!cur->isLeaf())
                cur->updateBranch();

            cur = nullptr;
        }

    } while (!sn.empty());
}

double AABBTree::treeHeight(Node* node)
{
    if (!node) return 0;
    return std::max(treeHeight(node->left), treeHeight(node->right)) + 1;
}

// inorder traversal of the tree and whenever come up with a leaf node, 
// find collided pairs corresponding to it.
void AABBTree::query()
{
    Node* cur = root;
    std::stack<Node*> sn;

    while (cur || !sn.empty())
    {
        while (cur)
        {
            sn.push(cur);
            cur = cur->left;
        }

        cur = sn.top();
        sn.pop();
        
        if (cur->isLeaf())
        {
            if (type == MotionState::STATIC)
                isProximity = queryProximity(cur);
            else
                isCollsn = queryCollision(cur);
            
            nodeSet.insert(cur);
        }
        
        cur = cur->right;
    }
}

// For AABB inside Node n, find all intersecting AABBs.
// Preorder traverse the tree and if find a collided node to be 
// (1) leaf, find a pair and add to the list
// (2) branch, push two children into the stack
bool AABBTree::queryProximity(Node* n)
{
    std::stack<Node*> sn;
    Node* cur = root;

    while (cur || !sn.empty())
    {
        while (cur)
        {
            if (cur->isCollid(n))
            {
                if (cur->isLeaf() && n != cur)
                {
                    if (nodeSet.find(cur) == nodeSet.end())
                    {
                        CD_HSE* a = cur->data->hse;
                        CD_HSE* b = n->data->hse;
                        if (!adjacentHSE(a,b))
                        {
                            if (getProximity(a,b))
                            {
                                interference_pairs.push_back(std::make_pair(a,b));
                                count++; 
                            }
                        }
                    }
                }

                sn.push(cur);
                cur = cur->left;
            }   
            else
            { 
                //if the AABB of the subtree does not collid with 
                //node n, we ignore the whole subtree
                break;
            }
        }

        if (sn.empty())
            break;

        cur = sn.top();
        sn.pop();
        cur = cur->right;
    }

    return count > 0;
}

bool AABBTree::queryCollision(Node* n)
{
    std::stack<Node*> sn;
    Node* cur = root;

    while (cur || !sn.empty())
    {
        while (cur)
        {
            if (cur->isCollid(n))
            {
                if (cur->isLeaf() && n != cur)
                {
                    if (nodeSet.find(cur) == nodeSet.end())
                    {
                        CD_HSE* a = cur->data->hse;
                        CD_HSE* b = n->data->hse;
                        if (!adjacentHSE(a,b))
                        {
                            if (getCollision(a,b)) 
                            {
                                interference_pairs.push_back(std::make_pair(a,b));
                                count++;
                            }
                        }
                    }
                }

                sn.push(cur);
                cur = cur->left;
            }
            else 
            {
                //if the AABB of the subtree does not collid with 
                //node n, we ignore the whole subtree
                break;
            }
        }
        
        if (sn.empty()) break;

        cur = sn.top();
        sn.pop();
        cur = cur->right;
    }

    return count > 0;
}

bool AABBTree::getCollision(const CD_HSE* a, const CD_HSE* b)
{
    if (gauss_seidel)
        return getCollisionGS(a,b);
    else
        return getCollisionJac(a,b);
}

std::vector<std::pair<CD_HSE*,CD_HSE*>> AABBTree::getInterferencePairs() const
{
    return interference_pairs;
}

