#ifndef BVH_H
#define BVH_H

#include "BVH_Node.h"

#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>


using Point_with_Node = std::pair<CGAL_Point,BVH_Node*>;
using Point_Node_Vector = std::vector<Point_with_Node>;

using pMap = CGAL::First_of_pair_property_map<Point_with_Node>;
using BV_HilbertSortingTraits = CGAL::Spatial_sort_traits_adapter_3<K,pMap>;


class BVH
{
    private:
       
        BVH_Node* root{nullptr};

        int num_tris{0};
        int num_bonds{0};
        std::vector<BVH_Node*> leaves;

        Point_Node_Vector children;
        BV_HilbertSortingTraits hst;
        int sort_iter{0};

        void processSurfaces(SURFACE**);
        void processCurves(CURVE**);

        void constructParentNodes();
        void constructRootNode();

        void initChildren();
        void sortChildren();

        const Point_Node_Vector getLeafSortingData() const;
        //const Point_Node_Vector getSortedLeafData() const;

        bool drawbool{false};
        std::string drawdir;
        void drawHeirarchyLevel() const;
        void writeHilbertCurveFiles(int level) const;

    protected:
        
        void buildHeirarchy();
        void constructLeafNodes(INTERFACE*);
        void constructLeafNodes(std::vector<Hse*>);
        
        static BVH_Node* createLeafNode(Hse* h);
        static BVH_Node* createInternalNode(BVH_Node* lc, BVH_Node* rc);

        static double proximityPad;
        //proximityPad is hardcoded in BVH.cpp,
        //where it must be initialized for now (static member data).
        //Too early to tell where/how this should be set.
        //1.0e-03 is default for static proximity boxes, and
        //1.0e-06 is default for kinetic collision boxes.

    public:

        explicit BVH(Front* front, bool draw = false);

        BVH() = default;
        ~BVH() = default;

        BVH(const BVH&) = delete;
        BVH& operator=(const BVH&) = delete;
        BVH(BVH&&) = delete;
        BVH& operator=(BVH&&) = delete;

        BVH_Node* const getRoot() const noexcept;
        const bool isEmpty() const noexcept; 

        void setDrawBool(bool draw);
        void setDrawDirectory(std::string dir);
};



const bool checkProximity(const BVH* A, const BVH* B);


#endif
