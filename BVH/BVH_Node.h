#ifndef BVH_NODE_H
#define BVH_NODE_H

#include "BoundingVolume.h"

#include <memory>
#include <utility>
#include <stack>


//Design Considerations:
//
//  1. Do we really need smart_ptrs for Node linkage?
//     No deletion or rotations are going to be performed,
//     and the objects will persist for the entirety of
//     the program.
//  2. How much of an impact will the use of smart_ptrs
//     have on performance?

using BoundingVolume = AABB;


class BVH_Node :
    public std::enable_shared_from_this<BVH_Node>
{
    private:

        BoundingVolume bv;
        std::weak_ptr<BVH_Node> parent;
        
    public:

        BVH_Node() = default;
        BVH_Node(BVH_Node&&) = default;
        BVH_Node& operator=(BVH_Node&&) = default;
        virtual ~BVH_Node() = default;

        BVH_Node(const BVH_Node&) = delete;
        BVH_Node& operator=(const BVH_Node&) = delete;

        virtual const bool isLeaf() const = 0;

        void setBV(BoundingVolume);
        const BoundingVolume& getBV() const;
        
        void setParent(std::shared_ptr<BVH_Node>);
        const std::weak_ptr<BVH_Node> getParent() const;

        const bool overlaps(const std::shared_ptr<BVH_Node>&) const;
        const double volume() const;

        virtual const Hse* const getHse() const;
        virtual void expandBV(double);

        virtual const std::weak_ptr<BVH_Node> getLeftChild() const;
        virtual const std::weak_ptr<BVH_Node> getRightChild() const;

        virtual void setChildren(std::shared_ptr<BVH_Node> lc,
                std::shared_ptr<BVH_Node> rc);
       
};


class InternalNode : public BVH_Node
{
    private:

        std::shared_ptr<BVH_Node> left;
        std::shared_ptr<BVH_Node> right;
        
        void setLeftChild(std::shared_ptr<BVH_Node> lc);
        void setRightChild(std::shared_ptr<BVH_Node> rc);

    public:

        InternalNode(std::shared_ptr<BVH_Node> lc,
                std::shared_ptr<BVH_Node> rc);

        InternalNode(InternalNode&&) = default;
        InternalNode& operator=(InternalNode&&) = default;
        ~InternalNode() = default;

        InternalNode() = delete;
        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;
        
        const bool isLeaf() const noexcept override;
        void expandBV(double) override;

        //Would be better if this could be made private,
        //or coupled to the constructor. However, this does
        //not seem possible given the use of smart_ptrs for
        //node linkage. See InternalNode constructor in
        //BVH_Node.cpp for details.
        void setChildren(std::shared_ptr<BVH_Node> lc,
                std::shared_ptr<BVH_Node> rc) override;
       
        const std::weak_ptr<BVH_Node> getLeftChild() const override;
        const std::weak_ptr<BVH_Node> getRightChild() const override;
};


class LeafNode : public BVH_Node
{
    private:

        Hse* hse{nullptr};

    public:

        explicit LeafNode(Hse* h);
        LeafNode(LeafNode&&) = default;
        LeafNode& operator=(LeafNode&&) = default;
        ~LeafNode() = default;

        LeafNode() = delete;
        LeafNode(const LeafNode&) = delete;
        LeafNode& operator=(const LeafNode&) = delete;

        const bool isLeaf() const noexcept override;
        const Hse* const getHse() const override;
};



#endif


