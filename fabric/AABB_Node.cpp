#include "AABB_Node.h"


/////////////////////////////////////////////////////////
////////           AABB_Node Methods            ////////
///////////////////////////////////////////////////////

void AABB_Node::setParent(AABB_Node* node) noexcept
{
    parent = node;
}

const AABB& AABB_Node::getAABB() const noexcept
{
    return bv;
}



/////////////////////////////////////////////////////////
////////          InternalNode Methods          ////////
///////////////////////////////////////////////////////

InternalNode::InternalNode(
        AABB_Node* const lc, AABB_Node* const rc)
{
    bv = AABB(lc->getAABB(),rc->getAABB());
    lc->setParent(this);
    left = lc;
    rc->setParent(this);
    right = rc;
}

const bool InternalNode::isLeaf() const noexcept
{
    return false;
}








/////////////////////////////////////////////////////////
////////            LeafNode Methods            ////////
///////////////////////////////////////////////////////

LeafNode::LeafNode(CD_HSE* const Hse, double Pad)
    : hse{Hse}, pad{Pad}
{
    bv = AABB(hse,pad);
}

LeafNode::LeafNode(CD_HSE* const Hse, double Dt, double Pad)
    : hse{Hse}, dt{Dt}, pad{Pad}
{
    bv = AABB(hse,dt,pad);
}

const bool LeafNode::isLeaf() const noexcept
{
    return true;
}

