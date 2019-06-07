#ifndef BVH_NODE_H
#define BVH_NODE_H

#include "BoundingVolume.h"

#include <utility>
#include <stack>


using BoundingVolume = AABB;



class BVH_Node
{
    private:

        BVH_Node* parent{nullptr};
    
    protected:

        //BoundingVolume bv;

    public:

        BoundingVolume bv;

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
        void expandBV(double);
        virtual void refitBV() = 0;

        virtual const Hse* const getHse() const noexcept;
        virtual const bool hasAdjacentHse(BVH_Node*) const noexcept;

        void setParent(BVH_Node* const p) noexcept;
        virtual void setChildren(BVH_Node* lc, BVH_Node* rc) noexcept;

        BVH_Node* getParent() const noexcept;
        virtual BVH_Node* getLeftChild() const noexcept;
        virtual BVH_Node* getRightChild() const noexcept;
       
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

        void setChildren(BVH_Node* lc, BVH_Node* rc) noexcept override;
        BVH_Node* getLeftChild() const noexcept override;
        BVH_Node* getRightChild() const noexcept override;

        void refitBV() override;
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
        const bool hasAdjacentHse(BVH_Node*) const noexcept override;

        void refitBV() override;
};


#endif


