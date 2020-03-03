#ifndef AABB_TREE_H
#define AABB_TREE_H

#include "AABB_Node.h"

#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>


using Point_Node_Pair = std::pair<CGAL_Point,AABB_Node*>;
using Point_Node_Vector = std::vector<Point_Node_Pair>;

using pMap = CGAL::First_of_pair_property_map<Point_Node_Pair>;
using AABB_HilbertSortingTraits = CGAL::Spatial_sort_traits_adapter_3<Kernel,pMap>;



class AABB_Tree
{
    public:

        MotionState mstate;

        explicit AABB_Tree(MotionState ms);

        AABB_Tree() = default;
        ~AABB_Tree() = default;

        //disable copy and move ctors
        AABB_Tree(const AABB_Tree&) = delete;
        AABB_Tree& operator=(const AABB_Tree&) = delete;
        AABB_Tree(AABB_Tree&&) = delete;
        AABB_Tree& operator=(AABB_Tree&&) = delete;

        static AABB_Node* createLeafNode(CD_HSE* const hse, double pad);
        static AABB_Node* createInternalNode(
                AABB_Node* const lc, AABB_Node* const rc);

        const bool isEmpty() const noexcept;
        const AABB_Node* const getRoot() const noexcept;

        void setFabricPad(double pad) noexcept;
        void setStringPad(double pad) noexcept;

        void buildTree(const std::vector<CD_HSE*>& hseList);


    private:

        AABB_Node* root {nullptr};

        double fabricPad {0.0};
        double stringPad {0.0};

        std::vector<AABB_Node*> leaves;

        Point_Node_Vector children;
        AABB_HilbertSortingTraits hst;
        int sort_iter {0};

        void constructLeafNodes(const std::vector<CD_HSE*>& hseList);

        void initChildren();
        Point_Node_Vector getLeafSortingData() const;
        void sortChildren();
            //void sortChildren(Point_Node_Vector& cvec);
        
        void constructParentNodes();
};




#endif



