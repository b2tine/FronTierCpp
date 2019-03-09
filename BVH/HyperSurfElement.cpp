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
//              HsBond Methods                  //
//////////////////////////////////////////////////

HsBond::HsBond(BOND* b)
    : bond{b}
{}

HsBond::HsBond(BOND* b, HseTag tag)
    : Hse(tag), bond{b}
{}

const POINT* const HsBond::Point_of_hse(int i) const
{
    assert(this->bond);
    assert( i >= 0 && i < 2 );
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

const POINT* const HsTri::Point_of_hse(int i) const
{
    assert(this->tri);
    assert( i >= 0 && i < 3 );
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
const bool areAdjacentHse(const Hse* const A, const Hse* const B)
{
    assert( A && B );

    int nA = A->num_pts();
    int nB = B->num_pts();
    
    for( int i = 0; i < nA; ++i )
    {
        const POINT* pA = A->Point_of_hse(i);
        for( int j = 0; j < nB; ++j )
        {
            if( pA == B->Point_of_hse(j) )
                return true;
        }
    }
    return false;
}

