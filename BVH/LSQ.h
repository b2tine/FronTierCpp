#ifndef LSQ_H
#define LSQ_H

#include <vector>
#include <cassert>



class LeastSquares2d
{
    public:

        void setA(double xx, double xy,
                  double yx, double yy)
        {
            A[0][0] = xx; A[0][1] = xy;
            A[1][0] = yx; A[0][1] = yy;
        }

        void setRHS(double b0, double b1)
        {
            b[0] = b0;
            b[1] = b1;
        }

        std::vector<double> solve()
        {
            assert(A[0][1] ==  A[1][0]);
            assert(A[0][0] > 0 && A[1][1] > 0);

            double D = Det2d(A[0][0], A[0][1],
                             A[1][0], A[1][1]);
            
            double D1 = Det2d(b[0], A[0][1],
                              b[1], A[1][1]);
            
            double D2 = Det2d(A[0][0], b[0],
                              A[1][0], b[1]);
            
            return std::vector<double> {D1/D, D2/D};
        }

    private:

        double A[2][2] = {{0,0},{0,0}};
        double b[2] = {0,0};

        double Det2d(double a, double b,
                     double c, double d)
        {
            return a*d - b*c;
        }

};

#endif
