#include "AABB_Tree.h"

AABB_Tree::AABB_Tree(MotionState ms)
    : mstate{ms}
{}

//AABB_Node* AABB_Tree::createLeafNode(CD_HSE* hse, double pad)
AABB_Node* AABB_Tree::createLeafNode(CD_HSE* const hse, double pad)
{
    assert(hse);
    return new LeafNode(hse,pad);
}

AABB_Node* AABB_Tree::createInternalNode(
        AABB_Node* const lc, AABB_Node* const rc)
{
    assert(lc && rc);
    return new InternalNode(lc,rc);
}

const bool AABB_Tree::isEmpty() const noexcept
{
    return (!root) ? true : false;
}

const AABB_Node* const AABB_Tree::getRoot() const noexcept
{
    return root;
}

void AABB_Tree::setFabricPad(double pad) noexcept
{
    fabricPad = pad;
}

void AABB_Tree::setStringPad(double pad) noexcept
{
    stringPad = pad;
}

void AABB_Tree::buildTree(const std::vector<CD_HSE*>& hseList)
{
    constructLeafNodes(hseList);

    initChildren();
    while( children.size() != 1 )
    {
        //alternate sorting direction at each level
        if( sort_iter % 2 == 0 )
        {
            std::reverse(children.begin(),children.end());
        }
        //drawHeirarchyLevel();//TODO
        constructParentNodes();
    }
   
    assert(root != nullptr);
    children.clear(); //Point_Node_Vector().swap(children);
    //drawbool = false;
}

void AABB_Tree::constructLeafNodes(const std::vector<CD_HSE*>& hseList)
{
    leaves.clear();
    leaves.reserve(hseList.size());

    std::vector<CD_HSE*>::const_iterator it;
    for (it = hseList.cbegin(); it != hseList.cend(); ++it)
    {
        double pad = fabricPad;
        if ((*it)->type == CD_HSE_TYPE::STRING_BOND)
            pad = stringPad;

        AABB_Node* leaf = AABB_Tree::createLeafNode(*it,pad);
        leaves.push_back(leaf);
    }
}

void AABB_Tree::initChildren()
{
    sort_iter = 0;
    children = getLeafSortingData();
    sortChildren();
}

Point_Node_Vector AABB_Tree::getLeafSortingData() const
{
    Point_Node_Vector leafdata;
    leafdata.reserve(leaves.size());

    std::vector<AABB_Node*>::const_iterator it;
    for (it = leaves.cbegin(); it != leaves.cend(); ++it)
    {
        AABB_Node* node = *it;
        Point_Node_Pair ctr_node_pair(node->getAABB().centroid(),node);
        leafdata.push_back(ctr_node_pair);
    }
    return leafdata;
}

void AABB_Tree::sortChildren()
{
    sort_iter++;
    if (children.size() == 1)
    {
        root = children[0].second;
        //drawHeirarchyLevel();//TODO
    }
    else
    {
        CGAL::hilbert_sort(children.begin(),children.end(),hst);
    }
}

void AABB_Tree::constructParentNodes()
{
    Point_Node_Vector parents;
    parents.reserve(children.size()/2 + 1);
    
    //greedily pair off sorted children
    for (int i = 0; i < children.size()-1; i += 2)
    {
        auto lc = children[i].second;
        auto rc = children[i+1].second;
        auto p = AABB_Tree::createInternalNode(lc,rc);
        Point_Node_Pair ctr_node_pair(p->getAABB().centroid(),p);
        parents.push_back(ctr_node_pair);
    }

    //if an odd number of children, move the unpaired
    //node up to parent level unchanged
    if (children.size() % 2 != 0)
    {
        auto orphan = children[children.size()-1].second;
        Point_Node_Pair ctr_node_pair(orphan->getAABB().centroid(),orphan);
        parents.push_back(ctr_node_pair);
    }

    std::swap(parents,children);
    children.shrink_to_fit();
    sortChildren();
}




