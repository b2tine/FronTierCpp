#include "BoundingVolume.h"

///////////////////////////////////
//////     AABB methods     //////
/////////////////////////////////

AABB::AABB()
    : lower(3,HUGE), upper(3,-HUGE)
{}

AABB::AABB(Hse* h)
{
    for( int i = 0; i < 3; ++i )
    {
        lower.push_back(h->min_coord(i));
        upper.push_back(h->max_coord(i));
    }
}

AABB::AABB(const BV_Point& L, const BV_Point& U)
    : lower{L}, upper{U}
{}

AABB::AABB(const AABB& A, const AABB& B)
{
    for( int i = 0; i < 3; ++i )
    {
        lower.push_back(std::min(A.lower[i],B.lower[i]));
        upper.push_back(std::max(A.upper[i],B.upper[i]));
    }
}

const BV_Type AABB::getBvType() const
{
    return BV_Type::AABB;
}

const CGAL_Point AABB::Centroid() const
{
    double ctr[3];
    for( int i = 0; i < 3; ++i )
        ctr[i] = 0.5*(lower[i]+upper[i]);
    return CGAL_Point(ctr[0],ctr[1],ctr[2]);
}

bool AABB::contains(const AABB& BB) const
{
    for( int i = 0; i < 3; ++i )
    {
        if( lower[i] >= BB.lower[i] ) return false;
        if( upper[i] <= BB.upper[i] ) return false;
    }
    return true;
}

bool AABB::overlaps(const AABB& BB) const
{
    for( int i = 0; i < 3; ++i )
    {
        if( lower[i] > BB.upper[i] ) return false;
        if( upper[i] < BB.lower[i] ) return false;
    }
    return true;
}

const double AABB::volume() const
{
    double volume = 1.0;
    for( int i = 0; i < 3; ++i )
        volume *= upper[i] - lower[i];
    return volume;

}

void AABB::print() const
{
    CGAL_Point ctr = this->Centroid();
    printf("\n   upper: (%3g,%3g,%3g) \n", upper[0], upper[1], upper[2]);
    printf("centroid: (%3g,%3g,%3g) \n", ctr.x(), ctr.y(), ctr.z());
    printf("   lower: (%3g,%3g,%3g) \n\n", lower[0], lower[1], lower[2]);
}


