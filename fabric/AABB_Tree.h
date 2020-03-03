#ifndef AABB_TREE_H
#define AABB_TREE_H

#include "AABB_Node.h"


class AABB_Tree
{
    public:

        AABB_Tree() = default;
        ~AABB_Tree() = default;

        //disable copy and move ctors
        AABB_Tree(const AABB_Tree&) = delete;
        AABB_Tree& operator=(const AABB_Tree&) = delete;
        AABB_Tree(AABB_Tree&&) = delete;
        AABB_Tree& operator=(AABB_Tree&&) = delete;

        const bool isEmpty() const noexcept;
        const AABB_Node* const getRoot() const noexcept;

        static AABB_Node* createLeafNode(const CD_HSE* const hse);
        static AABB_Node* createInternalNode(
                const AABB_Node* const lc, const AABB_Node* const rc);


    private:

        AABB_Node* root {nullptr};

};




#endif



