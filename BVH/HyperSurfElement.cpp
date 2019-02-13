#include "HyperSurfElement.h"

//////////////////////////////////////////////////
//              Hse Methods                     //
//////////////////////////////////////////////////

void Hse::setTag(HseTag Tag)
{
    tag = Tag;
}

HseTag Hse::getTag()
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
    : point{p}
{
    this->setTag(tag);
}

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
    : bond{b}
{
    this->setTag(tag);
}

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
    : tri{t}
{
    this->setTag(tag);
}

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

