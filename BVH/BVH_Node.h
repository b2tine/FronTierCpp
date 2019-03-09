#ifndef BVH_NODE_H
#define BVH_NODE_H

#include "BoundingVolume.h"

#include <utility>
#include <stack>


using BoundingVolume = AABB;



class BVH_Node
{
    private:

        BoundingVolume bv;
        BVH_Node* parent{nullptr};
        
    public:

        BVH_Node() = default;
        BVH_Node(BVH_Node&&) = default;
        BVH_Node& operator=(BVH_Node&&) = default;
        virtual ~BVH_Node() = default;

        BVH_Node(const BVH_Node&) = delete;
        BVH_Node& operator=(const BVH_Node&) = delete;

        virtual const bool isLeaf() const = 0;

        void setBV(BoundingVolume);
        const BoundingVolume getBV() const noexcept;
        const bool overlaps(const BVH_Node* const) const;
        const double volume() const noexcept;

        virtual const Hse* const getHse() const noexcept;
        virtual void expandBV(double);

        void setParent(BVH_Node* const p) noexcept;
        virtual void setChildren(BVH_Node* const lc,
                                 BVH_Node* const rc) noexcept;

        const BVH_Node* getParent() const noexcept;
        virtual const BVH_Node* const getLeftChild() const noexcept;
        virtual const BVH_Node* const getRightChild() const noexcept;
       
};


class InternalNode : public BVH_Node
{
    private:

        BVH_Node* left{nullptr};
        BVH_Node* right{nullptr};
        
        void setLeftChild(BVH_Node* lc) noexcept;
        void setRightChild(BVH_Node* rc) noexcept;

    public:

        InternalNode(BVH_Node* const lc, BVH_Node* const rc);

        InternalNode(InternalNode&&) = default;
        InternalNode& operator=(InternalNode&&) = default;
        ~InternalNode() = default;

        InternalNode() = delete;
        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;
        
        const bool isLeaf() const noexcept override;
        void expandBV(double) noexcept override;

        void setChildren(BVH_Node* const lc, BVH_Node* const rc) noexcept override;
        const BVH_Node* const getLeftChild() const noexcept override;
        const BVH_Node* const getRightChild() const noexcept override;
};


class LeafNode : public BVH_Node
{
    private:

        //TODO: consider unique_ptr<Hse> hse
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
        const Hse* const getHse() const noexcept override;
};


#endif


