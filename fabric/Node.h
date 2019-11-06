#ifndef NODE_H
#define NODE_H

#include "NewAABB.h"

#include <utility>
#include <stack>



class Node
{
    private:

        AABB bv;
        Node* parent{nullptr};
        
    public:

        Node() = default;
        Node(Node&&) = default;
        Node& operator=(Node&&) = default;
        virtual ~Node() = default;

        Node(const Node&) = delete;
        Node& operator=(const Node&) = delete;

        virtual const bool isLeaf() const = 0;

        void setBV(AABB);
        const AABB getBV() const noexcept;
        const bool overlaps(const Node* const) const;
        const double volume() const noexcept;

        //virtual const CD_HSE* const getHSE() const noexcept;
        //virtual const bool hasAdjacentHSE(Node*) const noexcept;
        //virtual void expandBV(double);

        void setParent(Node* const p) noexcept;
        //virtual void setChildren(Node* lc, Node* rc) noexcept;

        Node* getParent() const noexcept;
        //virtual Node* getLeftChild() const noexcept;
        //virtual Node* getRightChild() const noexcept;
       
};


class InternalNode : public Node
{
    private:

        Node* left{nullptr};
        Node* right{nullptr};
        
        void setLeftChild(Node* lc) noexcept;
        void setRightChild(Node* rc) noexcept;

    public:

        InternalNode(Node* const lc, Node* const rc);

        InternalNode(InternalNode&&) = default;
        InternalNode& operator=(InternalNode&&) = default;
        ~InternalNode() = default;

        InternalNode() = delete;
        InternalNode(const InternalNode&) = delete;
        InternalNode& operator=(const InternalNode&) = delete;
        
        const bool isLeaf() const noexcept override;
        //void expandBV(double) noexcept override;

        //void setChildren(Node* lc, Node* rc) noexcept override;
        //Node* getLeftChild() const noexcept override;
        //Node* getRightChild() const noexcept override;
};


class LeafNode : public Node
{
    private:

        CD_HSE* hse{nullptr};

    public:

        explicit LeafNode(CD_HSE* h);
        LeafNode(LeafNode&&) = default;
        LeafNode& operator=(LeafNode&&) = default;
        ~LeafNode() = default;

        LeafNode() = delete;
        LeafNode(const LeafNode&) = delete;
        LeafNode& operator=(const LeafNode&) = delete;

        const bool isLeaf() const noexcept override;
        
        //const Hse* const getHse() const noexcept override;
        //const bool hasAdjacentHse(Node*) const noexcept override;
};


#endif


