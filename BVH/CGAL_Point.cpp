#include "CGAL_Point.h"

CGAL_Point::CGAL_Point()
    : point(3,0.0)
{}

CGAL_Point::CGAL_Point(double x, double y, double z)
    : point{x,y,z}
{}

const double& CGAL_Point::operator[](int i) const
{
    assert(i >=0 && i <= 2);
    return point[i];
}

double& CGAL_Point::operator[](int i)
{
    const_cast<double&>(
            static_cast<const CGAL_Point&>(*this)[i]);
}

bool CGAL_Point::operator < (const CGAL_Point& rhs) const
{
    if( point[0] == rhs[0] )
    {
        if( point[1] == rhs[1] )
        {
            return point[2] < rhs[2];
        }
        return point[1] < rhs[1];
    }
    return point[0] < rhs[0];
}
