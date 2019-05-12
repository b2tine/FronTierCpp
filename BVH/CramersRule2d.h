#ifndef CRAMERS_RULE_2D_H
#define CRAMERS_RULE_2D_H

#include <vector>
#include <cassert>



class CramersRule2d
{
    public:

        void setA(double a00, double a01,
                  double a10, double a11) noexcept
        {
            A[0][0] = a00; A[0][1] = a01;
            A[1][0] = a10; A[1][1] = a11;
        }

        void setRHS(double b0, double b1) noexcept
        {
            b[0] = b0;
            b[1] = b1;
        }

        std::vector<double> solve() noexcept
        {
            double D = Determinant2d(A[0][0], A[0][1],
                                     A[1][0], A[1][1]);

            //Return empty vector if system is degenerate
            if( D < 1.0e-08 )
            {
                return {};
            }

            double D1 = Determinant2d(b[0], A[0][1],
                                      b[1], A[1][1]);
            
            double D2 = Determinant2d(A[0][0], b[0],
                                      A[1][0], b[1]);
            
            return std::vector<double> {D1/D, D2/D};
        }

    private:

        double A[2][2] = {{0,0},{0,0}};
        double b[2] = {0,0};

        double Determinant2d(double a00, double a01,
                             double a10, double a11) noexcept
        {
            return a00*a11 - a01*a10;
        }

};

#endif
