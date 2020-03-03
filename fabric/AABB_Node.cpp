#include "AABB_Node.h"


/////////////////////////////////////////////////////////
////////           AABB_Node Methods            ////////
///////////////////////////////////////////////////////

void AABB_Node::setParent(AABB_Node* const node) noexcept
{
    parent = p;
}

const AABB& AABB_Node::getAABB() const noexcept
{
    return bv;
}



/////////////////////////////////////////////////////////
////////          InternalNode Methods          ////////
///////////////////////////////////////////////////////

InternalNode::InternalNode(
        const AABB_Node* const lc, const AABB_Node* const rc)
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

LeafNode::LeafNode(const CD_HSE* const Hse, double Pad)
    : hse{Hse}, pad{Pad}
{
    bv = AABB(hse,pad);
}

const bool LeafNode::isLeaf() const noexcept
{
    return true;
}

