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

AABB AABB::merge(const AABB& ab) const {
    CPoint pl(3), pu(3);

    for (int i = 0; i < 3; i++) {
         pl[i] = std::min(lowerbound[i], ab.lowerbound[i]);
         pu[i] = std::max(upperbound[i], ab.upperbound[i]);
    }
    return AABB(pl, pu);
} 

double AABB::volume() {
    return (upperbound[0]-lowerbound[0])*(upperbound[1]-lowerbound[1])*
            (upperbound[2]-lowerbound[2]);
}

//This is the intersection test for AABB's.
//Not a collision or geometric primitive check.
bool AABB::isCollid(const AABB& ab)
{
    for (int i = 0; i < 3; ++i)
    {
        if (ab.upperbound[i] < lowerbound[i]) return false;
        if (ab.lowerbound[i] > upperbound[i]) return false;
    }
    return true;
    /*
    return (lowerbound[0] <= ab.upperbound[0] && upperbound[0] >= ab.lowerbound[0]) && 
           (lowerbound[1] <= ab.upperbound[1] && upperbound[1] >= ab.lowerbound[1]) && 
           (lowerbound[2] <= ab.upperbound[2] && upperbound[2] >= ab.lowerbound[2]); 
    */
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

bool AABB::contain(const AABB* ab) {
    return lowerbound[0] <= ab->lowerbound[0] && lowerbound[1] <= ab->lowerbound[1] &&
        lowerbound[2] <= ab->lowerbound[2] && upperbound[0] >= ab->upperbound[0] &&
        upperbound[1] >= ab->upperbound[1] && upperbound[2] >= ab->upperbound[2];
}

void Node::setBranch(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2,
            std::shared_ptr<Node> parent) {
    n1->parent = parent; 
    n2->parent = parent;
    left = n1;
    right = n2;
}

bool Node::isLeaf() {
    return left == nullptr && right == nullptr;
}

void Node::setLeaf(AABB* ab) {
    data.reset(ab);
}

void Node::updateBranch() {
    if (isLeaf())
        return;
    for (int i = 0; i < 3; i++) {
         box.lowerbound[i] = std::min(left->box.lowerbound[i], right->box.lowerbound[i]);
         box.upperbound[i] = std::max(left->box.upperbound[i], right->box.upperbound[i]);
    }
}

void Node::updateAABB() {
    if (isLeaf()) {
        box.lowerbound = data->lowerbound;
        box.upperbound = data->upperbound;
    }
    else {
        // branch node has no AABB yet
        if (box.lowerbound.size() == 0)
            box = left->box.merge(right->box);
        else 
            updateBranch();
    }
}

bool Node::isCollid(Node* n) {
    return box.isCollid(n->box);
}

Node* Node::getSibling() const {
    auto parent = this->parent.lock();
    if (!parent.get())
        return nullptr;
    return this == parent->left.get() ? parent->right.get() :
        parent->left.get();
}

Node::~Node()
{
    left.reset();
    right.reset();
    parent.reset();
    data.reset();
}

AABBTree::AABBTree(MotionState mstate)
    : type{mstate}
{}

AABBTree::~AABBTree()
{
    deleteTree();
}

void AABBTree::deleteTree()
{
    std::queue<std::shared_ptr<Node> > q;

    q.push(this->root);
    while (!q.empty())
    {
        auto node = q.front();
        q.pop();

        if (node->left != nullptr)
            q.push(node->left);
        
        if (node->right != nullptr)
            q.push(node->right);

        node.reset();
    }
}

void AABBTree::addAABB(AABB* ab) {
    if (root.get()) {
        auto node = std::make_shared<Node>();
        node->setLeaf(ab);
        node->updateAABB();
        insertNode(node, root);
        nodeArray.push_back(node);
        numLeaf++;
    }
    else {
        root = std::make_shared<Node>();
        root->setLeaf(ab);
        root->updateAABB();
        nodeArray.push_back(root);
        numLeaf++;
    }
}

// reorganize the tree structure
void AABBTree::updateTreeStructure() {
    root.reset();
    for (auto node : nodeArray) {
         if (root.get()) 
             insertNode(node, root);
         else 
             root = node;
    }
}

void AABBTree::insertNode(std::shared_ptr<Node> n, std::shared_ptr<Node>& parentNode) {
    std::shared_ptr<Node> p = parentNode;
    // if parent is a leaf node, then create a branch
    // with n and parent to be two children
    if (p->isLeaf()) {
        auto newParentNode = std::make_shared<Node>();
        
        newParentNode->parent = p->parent;
        auto par = p->parent.lock();

        if (par.get())
            par->left.get() == parentNode.get() ? par->left = newParentNode :
              par->right = newParentNode;
        newParentNode->setBranch(n, p, newParentNode);
        parentNode = newParentNode;
    }
    // we have to decide which subtree to insert to
    // the rule is insert to the subtree with smaller volume 
    else {
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
        if (vdiff1 < vdiff2) {
            insertNode(n, p->left);
        }
        else {
            insertNode(n, p->right);
        }
    }
    // this will guarantee all relavent ancestor will be 
    // updated
    parentNode->updateAABB();
}

//sets AABBTree::count = 0
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

    std::stack<Node*> sn;
    Node* cur = root.get();

    if (!root.get())
        return;
    
    //iterative postorder traversal of tree
    do {
        while (cur)
        {
            if (cur->right)
                sn.push(cur->right.get());
            sn.push(cur);
            cur = cur->left.get();
        }

        cur = sn.top();
        sn.pop();

        if (cur->right && !sn.empty() && cur->right.get() == sn.top())
        {
            sn.pop();
            sn.push(cur);
            cur = cur->right.get();
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


double AABBTree::treeHeight(Node* root) {
    if (!root)
        return 0;
    return std::max(treeHeight(root->left.get()),
                    treeHeight(root->right.get())) + 1;
}

// inorder traverse the tree and whenever come up with a leaf node, 
// find collided pairs correspond to it.
void AABBTree::query()
{
    Node* cur = root.get();
    std::stack<Node*> sn;

    while (cur || !sn.empty())
    {
        while (cur)
        {
            sn.push(cur);
            cur = cur->left.get();
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
        
        cur = cur->right.get();
    }
}

// For AABB inside Node n, find all intersecting AABBs.
// Preorder traverse the tree and if find a collided node to be 
// (1) leaf, find a pair and add to the list
// (2) branch, push two children into the stack
bool AABBTree::queryProximity(Node* n)
{
    std::stack<Node*> sn;
    Node* cur = root.get();

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
                                count++; 
                        }
                    }
                }

                sn.push(cur);
                cur = cur->left.get();
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
        cur = cur->right.get();
    }

    return count > 0;
}

bool AABBTree::queryCollision(Node* n)
{
    std::stack<Node*> sn;
    Node* cur = root.get();

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
                                count++;
                        }
                    }
                }

                sn.push(cur);
                cur = cur->left.get();
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
        cur = cur->right.get();
    }

    return count > 0;
}
