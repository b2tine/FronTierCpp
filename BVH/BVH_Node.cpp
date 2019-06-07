#include "BVH_Node.h"

/////////////////////////////////////////////////////////
////////           BVH_Node Methods             ////////
///////////////////////////////////////////////////////

void BVH_Node::setBV(BoundingVolume BV)
{
    bv = BV;
}

const BoundingVolume BVH_Node::getBV() const noexcept
{
    return bv;
}

const bool BVH_Node::overlaps(const BVH_Node* const node) const
{
    BoundingVolume bv = this->getBV();
    return bv.overlaps(node->getBV());
}

const double BVH_Node::volume() const noexcept
{
    return bv.volume();
}
        
BVH_Node* BVH_Node::getParent() const noexcept
{
    return parent;
}

void BVH_Node::setParent(BVH_Node* const p) noexcept
{
    parent = p;
}

void BVH_Node::expandBV(double pad)
{
    bv.expand(pad);
}

//InternalNode uses the following default implementations
const Hse* const BVH_Node::getHse() const noexcept
{
    return {};
}

const bool BVH_Node::hasAdjacentHse(BVH_Node* node) const noexcept
{
    return false;
}


//LeafNode uses the following default implementations
BVH_Node* BVH_Node::getLeftChild() const noexcept
{
    return {};
}

BVH_Node* BVH_Node::getRightChild() const noexcept
{
    return {};
}

void BVH_Node::setChildren(
        BVH_Node* lc, BVH_Node* rc) noexcept
{
    return;
}

/////////////////////////////////////////////////////////
////////         InternalNode Methods           ////////
///////////////////////////////////////////////////////

InternalNode::InternalNode(BVH_Node* lc, BVH_Node* rc)
{
    setBV(BoundingVolume(lc->getBV(),rc->getBV()));
    setChildren(lc,rc);
}

const bool InternalNode::isLeaf() const noexcept
{
    return false;
}

void InternalNode::setChildren(
        BVH_Node* lc, BVH_Node* rc) noexcept
{
    setLeftChild(lc);
    setRightChild(rc);
}

void InternalNode::setLeftChild(BVH_Node* lc) noexcept
{
    lc->setParent(this);
    left = lc;
}

void InternalNode::setRightChild(BVH_Node* rc) noexcept
{
    rc->setParent(this);
    right = rc;
}

BVH_Node* InternalNode::getLeftChild() const noexcept
{
    return left;
}

BVH_Node* InternalNode::getRightChild() const noexcept
{
    return right;
}

void InternalNode::refitBV()
{
    bv.encloseBVs(getLeftChild()->bv,getRightChild()->bv);
}

/////////////////////////////////////////////////////////
///////            LeafNode Methods             ////////
///////////////////////////////////////////////////////

LeafNode::LeafNode(Hse* h)
    : hse{h} 
{
    bv = BoundingVolume(hse);
}

const bool LeafNode::isLeaf() const noexcept
{
    return true;
}

const Hse* const LeafNode::getHse() const noexcept
{
    return hse;
}

const bool LeafNode::hasAdjacentHse(BVH_Node* node) const noexcept
{
    if( !node->isLeaf() )
        return false;
    return areAdjacentHse(this->getHse(),node->getHse());
}
        
void LeafNode::refitBV()
{
    bv.encloseHse(this->getHse());
}


