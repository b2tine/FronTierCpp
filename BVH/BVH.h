#ifndef BVH_H
#define BVH_H

#include "BVH_Node.h"

#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>


using Point_with_Node = std::pair<CGAL_Point,std::shared_ptr<BVH_Node>>;
using Point_Node_Vector = std::vector<Point_with_Node>;

using pMap = CGAL::First_of_pair_property_map<Point_with_Node>;
using BV_HilbertSortingTraits = CGAL::Spatial_sort_traits_adapter_3<K,pMap>;


class BVH
{
    private:
        
        std::shared_ptr<BVH_Node> root{nullptr};
        
        int sort_iter{0};
        BV_HilbertSortingTraits hst;
        std::vector<std::shared_ptr<BVH_Node>> leaves;
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
        static std::shared_ptr<BVH_Node> createLeafNode(Hse* h);
        static std::shared_ptr<BVH_Node> createInternalNode(
                std::shared_ptr<BVH_Node> lc, std::shared_ptr<BVH_Node> rc);

    public:

        BVH(const Front* const);
              
        BVH() = default;
        ~BVH() = default;

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        const std::weak_ptr<BVH_Node> getRoot() const noexcept;
        const bool isEmpty() const noexcept; 

        //temp functions for testing/debugging
        void buildTester(std::vector<Hse*>);
        void writeHilbertCurveFile(std::string) const;
};


const bool checkProximity(BVH*, BVH*);

std::stack< std::pair<std::shared_ptr<BVH_Node>,
    std::shared_ptr<BVH_Node>> >
    queryProximity(std::shared_ptr<BVH_Node>,std::shared_ptr<BVH_Node>);



#endif
