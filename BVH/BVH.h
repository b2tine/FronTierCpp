#ifndef BVH_H
#define BVH_H

#include "BVH_Node.h"


class BVH
{
    private:
        
        //std::shared_ptr<InternalNode> root{nullptr};
        //int lastHseCount{0};
        //std::vector<Hse*> hseList;
        int lastCountLeaves{0};
        std::vector<std::shared_ptr<BVH_Node>> leaves;

        std::vector<CGAL_Point> ctrVec;
        std::map<CGAL_Point,std::shared_ptr<BVH_Node>> bvMap;

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
        void writeHilbertCurve(std::string,std::string);
        
};




#endif
