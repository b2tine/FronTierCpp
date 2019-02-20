#include "BVH_Node.h"

/////////////////////////////////////////////////////////
////////           BVH_Node Methods             ////////
///////////////////////////////////////////////////////
        
void BVH_Node::setBV(BoundingVolume BV)
{
    bv = std::move(BV);
}
const BoundingVolume& BVH_Node::getBV() const
{
    return bv;
}

void BVH_Node::setParent(std::shared_ptr<BVH_Node> P)
{
    parent = std::move(P);
}
        
const std::weak_ptr<BVH_Node> BVH_Node::getParent() const
{
    return std::weak_ptr<InternalNode>(parent);
}

const bool BVH_Node::overlaps(const std::shared_ptr<BVH_Node>& node) const
{
    auto bv = getBV();
    return bv.overlaps(node->getBV());
}


/////////////////////////////////////////////////////////
////////         InternalNode Methods           ////////
///////////////////////////////////////////////////////

InternalNode::InternalNode(std::shared_ptr<BVH_Node> lc,
        std::shared_ptr<BVH_Node> rc)
{
    assert(lc && rc);
    setBV(BoundingVolume(lc->getBV(),rc->getBV()));
}

//SetChildren() is a temporary solution for testing.
//It is not possible to call shared_from_this() in the
//constructor of InternalNode since an existing shared_ptr
//managing "this" must already exist.
//The above constructor and SetChildren() are consolidated
//inside a static factory function of the BVH class,
//but setChildren(), unfortunately, remains exposed to the
//public interface.
//
//Explicitly construct the shared_ptr instead?
void InternalNode::setChildren(std::shared_ptr<BVH_Node> lc,
        std::shared_ptr<BVH_Node> rc)
{
    setLeftChild(std::move(lc));
    setRightChild(std::move(rc));
}

void InternalNode::setLeftChild(std::shared_ptr<BVH_Node> lc)
{
    lc->setParent(shared_from_this());
    left = std::move(lc);
}

void InternalNode::setRightChild(std::shared_ptr<BVH_Node> rc)
{
    rc->setParent(shared_from_this());
    right = std::move(rc);
}

const std::weak_ptr<BVH_Node> InternalNode::getLeftChild() const
{
    return std::weak_ptr<BVH_Node>(left);
}

const std::weak_ptr<BVH_Node> InternalNode::getRightChild() const
{
    return std::weak_ptr<BVH_Node>(right);
}

const bool InternalNode::isLeaf() const
{
    return false;
}


/////////////////////////////////////////////////////////
///////            LeafNode Methods             ////////
///////////////////////////////////////////////////////

LeafNode::LeafNode(Hse* h)
    : hse{h} 
{
    setBV(BoundingVolume(h));
}

const Hse* const LeafNode::getHse() const
{
    return hse;
}

const bool LeafNode::isLeaf() const
{
    return true;
}



