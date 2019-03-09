#ifndef BVH_H
#define BVH_H

/*
#include "BVH_Node.h"

#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>


using Point_with_Node = std::pair<CGAL_Point,const BVH_Node* const>;
using Point_Node_Vector = std::vector<Point_with_Node>;

using pMap = CGAL::First_of_pair_property_map<Point_with_Node>;
using BV_HilbertSortingTraits = CGAL::Spatial_sort_traits_adapter_3<K,pMap>;


class BVH
{
    private:
       
        //TODO: consider unique_ptr for root
        BVH_Node* root{nullptr};
        std::vector<const BVH_Node*> leaves;
        
        int sort_iter{0};
        BV_HilbertSortingTraits hst;
        Point_Node_Vector children;

        void buildHeirarchy();
        void constructLeafNodes(const INTERFACE* const intfc);
        void constructParentNodes();
        void constructRootNode();
        //void clearVectors();
        
        void initChildren();
        void sortChildren();

        const Point_Node_Vector getLeafSortingData() const;
        const Point_Node_Vector getSortedLeafData() const;

        //hard coded in BVH.cpp, where it must be initialized, for now;
        //too early to tell where/how this should be set.
        //1.0e-03 is default for static proximity boxes, and
        //1.0e-06 is default for kinetic collision boxes.
        static double proximityPad;
        
        static BVH_Node* createLeafNode(Hse* h);
        static BVH_Node* createInternalNode(
                const BVH_Node* lc, const BVH_Node* rc);

    public:

        explicit BVH(const Front* const);
              
        BVH() = default;
        ~BVH() = default;

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        const BVH_Node* const getRoot() const noexcept;
        const bool isEmpty() const noexcept; 

        //temp functions for testing/debugging
        void buildTester(std::vector<Hse*>);
        void writeHilbertCurveFile(std::string,int) const;
};



const bool checkProximity(const BVH* const , const BVH* const);
*/


#endif
