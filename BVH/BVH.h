#ifndef BVH_H
#define BVH_H

#include "BVH_Node.h"

#include <CGAL/hilbert_sort.h>
//#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>

#include <utility>


using Point_with_Node = std::pair<CGAL_Point,std::shared_ptr<BVH_Node>>;
using Point_Node_Vector = std::vector<Point_with_Node>;

using pMap = CGAL::First_of_pair_property_map<Point_with_Node>;
using BV_HilbertSortingTraits = CGAL::Spatial_sort_traits_adapter_3<K,pMap>;


class BVH
{
    private:
        
        //std::shared_ptr<InternalNode> root{nullptr};
        //int lastHseCount{0};
        //std::vector<Hse*> hseList;
        int lastCountLeaves{0};
        //keep leaves in own vector for now
        std::vector<std::shared_ptr<BVH_Node>> leaves;
        std::vector<std::shared_ptr<BVH_Node>> children;
        std::vector<std::shared_ptr<BVH_Node>> parents;

        Point_Node_Vector centroids;

    public:

        //TODO: Will want to enforce the invariant that
        //      the root is always initialized, but this
        //      requires a build routine since we are
        //      performing a bottom up construction.
              
        BVH() = default;
        ~BVH() = default;

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        /*
        const std::weak_ptr<const InternalNode> getRoot() const
        {
            return std::weak_ptr<InternalNode>(root);
        }
        */

        static std::shared_ptr<LeafNode> createLeafNode(Hse* h);
        static std::shared_ptr<InternalNode> createInternalNode(
                std::shared_ptr<BVH_Node> lc, std::shared_ptr<BVH_Node> rc);

        //void assembleHseListFromInterface(const INTERFACE* const intfc);
        //void clearHseList();

        void constructLeafNodes(const INTERFACE* const intfc);
        void clearLeafNodes();

        void sortNodes();
        //void constructParentNodes();

        //temp function for prototype debugging
        void writeHilbertCurveFile(std::string,std::string);
        
};




#endif
