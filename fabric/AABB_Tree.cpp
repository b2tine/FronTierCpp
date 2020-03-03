#include "AABB_Tree.h"


const bool AABB_Tree::isEmpty() const noexcept
{
    return (!root) ? true : false;
}

const AABB_Node* const getRoot() const noexcept
{
    return root;
}
