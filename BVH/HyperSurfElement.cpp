#include "HyperSurfElement.h"

//////////////////////////////////////////////////
//              Hse Methods                     //
//////////////////////////////////////////////////

Hse::Hse(HseTag Tag)
    : tag{Tag}
{}

/*
void Hse::setTag(HseTag Tag)
{
    tag = Tag;
}
*/

const HseTag Hse::getTag() const noexcept
{
    return tag;
}

//////////////////////////////////////////////////
//              HsPoint Methods                 //
//////////////////////////////////////////////////

HsPoint::HsPoint(POINT* p)
    : point{p}
{}

HsPoint::HsPoint(POINT* p, HseTag tag)
    : Hse(tag), point{p}
{}

POINT* HsPoint::Point_of_hse(int i) const
{
    assert( i == 0 );
    assert( point != nullptr );
    return point;
}

double HsPoint::min_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return Coords(this->Point_of_hse(0))[dim];
}

double HsPoint::max_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return Coords(this->Point_of_hse(0))[dim];
}

//////////////////////////////////////////////////
//              HsBond Methods                  //
//////////////////////////////////////////////////

HsBond::HsBond(BOND* b)
    : bond{b}
{}

HsBond::HsBond(BOND* b, HseTag tag)
    : Hse(tag), bond{b}
{}

POINT* HsBond::Point_of_hse(int i) const
{
    assert( i >= 0 && i < this->num_pts() );
    assert( bond != nullptr );
    return i == 0 ? bond->start : bond->end;
}

double HsBond::min_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return std::min(Coords(this->Point_of_hse(0))[dim],
            Coords(this->Point_of_hse(1))[dim]);
}

double HsBond::max_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    return std::max(Coords(this->Point_of_hse(0))[dim],
            Coords(this->Point_of_hse(1))[dim]);
}

//////////////////////////////////////////////////
//              HsTri Methods                   //
//////////////////////////////////////////////////

HsTri::HsTri(TRI* t)
    : tri{t}
{}

HsTri::HsTri(TRI* t, HseTag tag)
    : Hse(tag), tri{t}
{}

POINT* HsTri::Point_of_hse(int i) const
{
    assert( i >= 0 && i < this->num_pts() );
    assert( tri != nullptr );
    return Point_of_tri(tri)[i];
}

double HsTri::min_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    double val = HUGE;
    for( int i = 0; i < 3; ++i )
    {
        val = std::min(val,
                Coords(this->Point_of_hse(i))[dim]);
    }
    return val;
}

double HsTri::max_coord(int dim) const
{
    assert( dim >= 0 && dim < 3 );
    double val = -HUGE;
    for( int i = 0; i < 3; ++i )
    {
        val = std::max(val,
                Coords(this->Point_of_hse(i))[dim]);
    }
    return val;
}

//////////////////////////////////////////////////////////

//TODO: should probably move this into BVH class
const bool areAdjacentHse(Hse* A, Hse* B)
{
    assert( A && B )

    //NOTE: These first 2 conditionals evaluate true
    //      for the purpose of avoiding false positive
    //      proximity queries. 
    if( A == B)
        return true;

    int n1 = A->num_points();
    int n2 = A->num_points();
    
    //POINTS don't have a BV
    if( n1 == 1 || n2 == 1 )
        return true;

    POINT* ptsA[2];
    POINT* ptsB[2];
    
    for( int i = 0; i < 2; ++i )
    {
        ptsA[i] = A->Point_of_hse(i);
        ptsB[i] = B->Point_of_hse(i);
    }
    
    int count = 0;
    for( int i = 0; i < 2; ++i )
    {
        for( int j = 0; j < 2; ++j )
        {
            if( ptsA[i] == ptsB[j] )
                count++;
        }
    }
    
    if( count == 0 )
        return false;

    if( n1 == 2 || n2 == 2 )
    {
        //bond-bond or bond-tri: 1 common point
        if( count == 1 )
        {
            return true;
        }
        return false;
    }
    else
    {
        //tri-tri: 2 common points
        if( count == 2 )
            return true;

        POINT* Ap3 = A->Point_of_hse(3);
        POINT* Bp3 = B->Point_of_hse(3);
    
        for( int i = 0; i < 2; ++i )
        {
            if( ptsA[i] == Bp3 || ptsB[i] == Ap3)
                return true;
        }
        return false;
    }

}

