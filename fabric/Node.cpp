#include "Node.h"

/////////////////////////////////////////////////////////
////////           Node Methods             ////////
///////////////////////////////////////////////////////

void Node::setBV(AABB BV)
{
    bv = BV;
}
const AABB Node::getBV() const noexcept
{
    return bv;
}

const bool Node::overlaps(const Node* const node) const
{
    AABB bv = getBV();
    return bv.overlaps(node->getBV());
}

const double Node::volume() const noexcept
{
    return bv.volume();
}
        
Node* Node::getParent() const noexcept
{
    return parent;
}

void Node::setParent(Node* const p) noexcept
{
    parent = p;
}

//InternalNode uses these default implementations
const Hse* const Node::getHse() const noexcept
{
    return {};
}

const bool Node::hasAdjacentHse(Node* node) const noexcept
{
    return false;
}

//LeafNode uses these default implementations
void Node::expandBV(double pad)
{
    bv.expand(pad);
}

Node* Node::getLeftChild() const noexcept
{
    return {};
}

Node* Node::getRightChild() const noexcept
{
    return {};
}

void Node::setChildren(
        Node* lc, Node* rc) noexcept
{
    return;
}

/////////////////////////////////////////////////////////
////////         InternalNode Methods           ////////
///////////////////////////////////////////////////////

InternalNode::InternalNode(Node* lc, Node* rc)
{
    setBV(BoundingVolume(lc->getBV(),rc->getBV()));
    setChildren(lc,rc);
}

const bool InternalNode::isLeaf() const noexcept
{
    return false;
}

void InternalNode::expandBV(double pad) noexcept
{
     return;
}

void InternalNode::setChildren(
        Node* lc, Node* rc) noexcept
{
    setLeftChild(lc);
    setRightChild(rc);
}

void InternalNode::setLeftChild(Node* lc) noexcept
{
    lc->setParent(this);
    left = lc;
}

void InternalNode::setRightChild(Node* rc) noexcept
{
    rc->setParent(this);
    right = rc;
}

Node* InternalNode::getLeftChild() const noexcept
{
    return left;
}

Node* InternalNode::getRightChild() const noexcept
{
    return right;
}

/////////////////////////////////////////////////////////
///////            LeafNode Methods             ////////
///////////////////////////////////////////////////////

LeafNode::LeafNode(Hse* h)
    : hse{h} 
{
    setBV(BoundingVolume(h));
}

const bool LeafNode::isLeaf() const noexcept
{
    return true;
}

const Hse* const LeafNode::getHse() const noexcept
{
    return hse;
}

const bool LeafNode::hasAdjacentHse(Node* node) const noexcept
{
    if( !node->isLeaf() )
        return false;
    return areAdjacentHse(this->getHse(),node->getHse());
}

