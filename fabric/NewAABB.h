#ifndef AABB_H
#define AABB_H

#include "CD_HSE.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGAL_Point = K::Point_3;

using BV_Point = std::vector<double>;


class AABB
{
    public:

        BV_Point lower;
        BV_Point upper;

        AABB();
        explicit AABB(CD_HSE*);
        AABB(const AABB&,const AABB&);

        AABB(const AABB&) = default;
        AABB& operator=(const AABB&) = default;
        AABB(AABB&&) = default;
        AABB& operator=(AABB&&) = default;
        ~AABB() = default;

        const CGAL_Point Centroid() const;
        const double volume() const noexcept;
        void expand(double pad);

        const bool overlaps(const AABB&) const;
        const bool contains(const AABB&) const;

        void print() const;
};


#endif
