#include "AABB_Tree.h"


const bool AABB_Tree::isEmpty() const noexcept
{
    return (!root) ? true : false;
}

const AABB_Node* const AABB_Tree::getRoot() const noexcept
{
    return root;
}

//AABB_Node* AABB_Tree::createLeafNode(CD_HSE* hse, double pad)
AABB_Node* AABB_Tree::createLeafNode(const CD_HSE* const hse, double pad)
{
    assert(hse);
    return new LeafNode(hse,pad);
}

AABB_Node* AABB_Tree::createInternalNode(
        const AABB_Node* const lc, const AABB_Node* const rc)
{
    assert(lc && rc);
    return new InternalNode(lc,rc);
}
