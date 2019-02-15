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
        
        std::shared_ptr<InternalNode> root{nullptr};
        
        int numLeaves{0};
        std::vector<std::shared_ptr<BVH_Node>> leaves;

        Point_Node_Vector children;
        BV_HilbertSortingTraits hst;

        int sort_iter{0};
        void sortChildNodes();
        void constructLeafNodes(const INTERFACE* const intfc);
        void constructParentNodes();
        void constructRootNode();
        //void clearVectors();

        static std::shared_ptr<LeafNode> createLeafNode(Hse* h);
        static std::shared_ptr<InternalNode> createInternalNode(
                std::shared_ptr<BVH_Node> lc, std::shared_ptr<BVH_Node> rc);

    public:

        BVH(const Front* const);
              
        BVH() = default;
        ~BVH() = default;

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        const std::weak_ptr<InternalNode> getRoot() const;

        //temp function for testing/debugging
        void writeHilbertCurveFile(std::string,std::string);
        void buildTester(std::vector<Hse*>);
};



#endif
