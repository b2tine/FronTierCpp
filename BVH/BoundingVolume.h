#ifndef BOUNDING_VOLUME_H
#define BOUNDING_VOLUME_H

#include "HyperSurfElement.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGAL_Point = K::Point_3;



enum class BV_Type
{
    AABB,
    OBB,
    KDOP,
    SPHERE
};


using BV_Point = std::vector<double>;


//Axis Aligned Bounding Box (AABB)
class AABB
{
    public:

        //TODO: would like to make these private,
        //      but it makes some testing difficult..
        BV_Point lower;
        BV_Point upper;

        AABB();
        explicit AABB(Hse*);
        AABB(const AABB&,const AABB&);

        AABB(const AABB&) = default;
        AABB& operator=(const AABB&) = default;
        AABB(AABB&&) = default;
        AABB& operator=(AABB&&) = default;
        ~AABB() = default;

        //const BV_Type getBvType() const noexcept;
        const CGAL_Point Centroid() const;
        const double volume() const noexcept;
        void expand(double pad);

        //TODO: need to carefully test this when
        //      considering self intersection checks
        const bool overlaps(const AABB&) const;
        const bool contains(const AABB&) const;

        void print() const;
};


#endif
