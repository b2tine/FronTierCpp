#ifndef IF_TURB_H
#define IF_TURB_H

//TODO: Not finding this 
struct SpaldingWallLaw
{
    double u;           //tangential fluid velocity magnitude
    double y;           //distance to wall
    double nu;          //kinematic viscosity (laminar)
    double B {5.2};     //model constant ~ 5.2 for smooth walls

    SpaldingWallLaw(double u_tan, double dist, double nu_lam)
        : u(u_tan), y(dist), nu(nu_lam)
    {}

    SpaldingWallLaw(double u_tan, double dist, double nu_lam, double b)
        : u(u_tan), y(dist), nu(nu_lam), B(b)
    {}

    //u_star is the friction velocity
    double operator() (double u_star) const
    {
        //TODO: dividing by zero whe u_star = 0 is an initial guess,
        //      need a better formulation of the root finding problem.
        double Kup = 0.41*u/u_star;
        return exp(0.41*B)*(y/nu*u_star - u/u_star)
            + 1.0 + Kup + 0.5*Kup*Kup + Kup*Kup*Kup/6.0 - exp(Kup);
    }
};


template<typename F>
double secantMethod(const F& f, double xa, double xb)
{
    const int MAXITER = 5000;
    const double TOL = 1.0e-08;

    double x0 = xa;
    double x1 = xb;
    //TODO: handle input edge cases

    double xn = x0;
    for (int i = 0; i < MAXITER; ++i)
    {
        double fval = f(xn);
        if (fabs(fval) < TOL)
            return xn;

        xn = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
        x0 = x1;
        x1 = xn;
    }

    //TODO: better debugging output
    printf("ERROR: secantMethod() could not find root\n");
    LOC(); clean_up(EXIT_FAILURE);
}


double computeWallShearStress(double u_tan, double walldist, double mu, double rho);
double computeFrictionVelocity(double u_tan, double walldist, double mu, double rho);


#endif
