#ifndef CGAL_POINT_H
#define CGAL_POINT_H

#include <CGAL/hilbert_sort.h>
#include <vector>


class CGAL_Point
{
    private:

        std::vector<double> point;
    
    public:

        CGAL_Point();
        CGAL_Point(double, double, double);

        ~CGAL_Point() = default;
        CGAL_Point(const CGAL_Point&) = default;
        CGAL_Point(CGAL_Point&&) = default;

        CGAL_Point& operator=(const CGAL_Point&) = delete;
        CGAL_Point& operator=(CGAL_Point&&) = delete;

        const double& operator[](int i) const;
        double& operator[](int i);

        bool operator < (const CGAL_Point& rhs) const;

};

//additional structures allowing CGAL_Point to
//be used in the CGAL function hilbert_sort();
struct BV_LessX
{
    bool operator()(const CGAL_Point& p, const CGAL_Point& q) const
    {
        return p[0] < q[0];
    }
};

struct BV_LessY
{
    bool operator()(const CGAL_Point& p, const CGAL_Point& q) const
    {
        return p[1] < q[1];
    }
};

struct BV_LessZ
{
    bool operator()(const CGAL_Point& p, const CGAL_Point& q) const
    {
        return p[2] < q[2];
    }
};

struct BV_HilbertSortingTraits
{
    using Point_3 = CGAL_Point;
    using Less_x_3 = BV_LessX;
    using Less_y_3 = BV_LessY;
    using Less_z_3 = BV_LessZ;

    Less_x_3 less_x_3_object() const
    {
        return Less_x_3();
    }

    Less_y_3 less_y_3_object() const
    {
        return Less_y_3();
    }

    Less_z_3 less_z_3_object() const
    {
        return Less_z_3();
    }
};

#endif
