#ifndef IF_TURB_H
#define IF_TURB_H

#include "iFluid.h"

struct SpaldingWallLaw
{
    double u;       //tangential fluid velocity
    double y;       //distance to wall
    double nu;      //kinematic viscosity (laminar)
    double B;       //model constant ~ 5.2 for smooth walls

    SpaldingWallLaw(double v_tan, double dist, double nu_lam, double b)
        : u(v_tan), y(dist), nu(nu_lam), B(b)
    {}

    double operator() (double u_star) const
    {
        double Kup = 0.41*u/u_star;
        return exp(0.41*B)*(y/nu*u_star - u/u_star)
            + 1.0 + Kup + 0.5*Kup*kup + Kup*Kup*Kup/6.0 - exp(Kup);
    }
};


template<typename F>
double secantMethod(const F& f, double xa, double xb)
{
    const int MAXITER = 2000;
    const double TOL = 1.0e-08;

    double x0 = xa;
    double x1 = xb;
    //TODO: handle input edge cases

    double xn = x0;
    for (int i = 0; i < MAXITER; ++i)
    {
        xn = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
        x0 = x1;
        x1 = xn;

        double fval = f(xn);
        if (fabs(fval) < TOL)
            return xn;
    }

    //TODO: better debugging output
    printf("ERROR: secantMethod() could not find root\n");
    LOC(); clean_up(EXIT_FAILURE);
}

#endif
