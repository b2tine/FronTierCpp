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

/*
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
*/

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


//adjacent elements have at least 1 common point
const bool areAdjacentHse(Hse* A, Hse* B)
{
    assert( A && B );
    if( A == B) return true;

    int nA = A->num_pts();
    int nB = A->num_pts();
    
    for( int i = 0; i < nA; ++i )
    {
        POINT* pA = A->Point_of_hse(i);
        for( int j = 0; j < nB; ++j )
        {
            if( pA == B->Point_of_hse(j) )
                return true;
        }
    }
}

