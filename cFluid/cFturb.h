#ifndef CF_TURB_H
#define CF_TURB_H

struct SpaldingWallLaw
{
public:

    SpaldingWallLaw(double u_tan, double dist, double nu_lam)
        : u{u_tan}, y{dist}, nu{nu_lam}
    {}

    SpaldingWallLaw(double u_tan, double dist, double nu_lam, double b)
        : u{u_tan}, y{dist}, nu{nu_lam}, B{b}
    {}

    ~SpaldingWallLaw() = default;

    SpaldingWallLaw() = delete;
    SpaldingWallLaw(const SpaldingWallLaw&) = delete;
    SpaldingWallLaw& operator=(const SpaldingWallLaw&) = delete;
    SpaldingWallLaw(SpaldingWallLaw&&) = delete;
    SpaldingWallLaw& operator=(SpaldingWallLaw&&) = delete;

        
    const int MAXITER = 5000;
    const double TOL = 1.0e-08;
        
    //Newton's method to solve for friction velocity
    double solve(double u0)
    {
        if (u < TOL) return 0.0;

        double un = u0;
        for (int i = 0; i < MAXITER; ++i)
        {
            double fval = f(un);
            if (fabs(fval) < TOL) return un;

            un = un - f(un)/fprime(un);
        }

        printf("\nERROR: SpaldingWallLaw::solve() could not find root\n");
        printf("u = %g  y = %g  nu = %g\n",u,y,nu);
        printf("un = %g   f(un) = %g\n",un,f(un));
        LOC(); clean_up(EXIT_FAILURE);
    }

private:
    
    double u;           //tangential fluid velocity magnitude
    double y;           //physical distance to wall
    double nu;          //fluid kinematic viscosity (laminar)
    double B {5.2};     //model constant ~ 5.2 for smooth walls
    double K {0.41};    //Karman constant

    double f(double u_plus)
    {
        double kup = K*u_plus;

        double val = y*u - nu*u_plus*u_plus
                    - nu*exp(-K*B)*((exp(kup) - 1.0) - kup
                            - 0.5*kup*kup - 1.0/6.0*kup*kup*kup)*u_plus;
        
        return val;
    }

    double fprime(double u_plus)
    {
        double kup = K*u_plus;

        double val = -2.0*nu*u_plus - nu*exp(-K*B)*((exp(kup) - 1.0)
                        - kup - 0.5*kup*kup - 1.0/6.0*kup*kup*kup)
                    - nu*u_plus*exp(-K*B)*(K*exp(kup) - K - K*kup - K*0.5*kup*kup);

        return val;
    }

};


double computeWallShearStress(double u_tan, double walldist,
            double mu, double rho, double u_wall_initial_guess);

double computeFrictionVelocity(double u_tan, double walldist,
            double nu, double u_wall_initial_guess);

double computeWallVelocity(double u_tan, double walldist,
            double nu, double u_wall_initial_guess);


#endif
