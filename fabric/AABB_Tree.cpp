#include "AABB_Tree.h"

AABB_Tree::AABB_Tree(MotionState ms)
    : mstate{ms}
{}

const bool AABB_Tree::isEmpty() const noexcept
{
    return (!root) ? true : false;
}

const AABB_Node* const AABB_Tree::getRoot() const noexcept
{
    return root;
}

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
