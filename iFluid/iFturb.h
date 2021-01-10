#ifndef IF_TURB_H
#define IF_TURB_H

//TODO: Need to be able to solve this root
struct SpaldingWallLaw
{
public:

    SpaldingWallLaw(double u_tan, double dist, double nu_lam)
        : u(u_tan), y(dist), nu(nu_lam)
    {}

    SpaldingWallLaw(double u_tan, double dist, double nu_lam, double b)
        : u(u_tan), y(dist), nu(nu_lam), B(b)
    {}

    ~SpaldingWallLaw() = default;

    SpaldingWallLaw() = delete;
    SpaldingWallLaw(const SpaldingWallLaw&) = delete;
    SpaldingWallLaw& operator=(const SpaldingWallLaw&) = delete;
    SpaldingWallLaw(SpaldingWallLaw&&) = delete;
    SpaldingWallLaw& operator=(SpaldingWallLaw&&) = delete;

        
    const int MAXITER = 5000;
    const double TOL = 1.0e-08;
        
    double solve(double u0)
    {
            //printf("u = %g  y = %g\n",u,y);//DEBUG
        if (u < TOL) return 0.0;
            //assert(u0 > 0.0);

        double un = u0;
        for (int i = 0; i < MAXITER; ++i)
        {
            double fval = f(un);
                //printf("iter = %d  un = %g  fval = %g\n",i,un,fval);//DEBUG
            if (fabs(fval) < TOL) return un;

            un = un - f(un)/fprime(un);
            
            //TODO: check if nan, inf, or less than zero
        }

        //TODO: better debugging output
        printf("ERROR: SpaldingWallLaw::solve() could not find root\n");
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
        double val = u_plus*u_plus - y/nu*u
            + exp(-K*B)*((exp(K*u_plus) - 1.0)*u_plus
                - K*u_plus*u_plus - 0.5*K*K*u_plus*u_plus*u_plus
                - K*K*K/6.0*u_plus*u_plus*u_plus*u_plus);
        return val;
    }

    double fprime(double u_plus)
    {
        double val = 2.0*u_plus
            + exp(-K*B)*(exp(K*u_plus)*(K*u_plus + 1.0) - 1.0
                    - 2.0*K*u_plus - 1.5*K*K*u_plus*u_plus
                    - 2.0/3.0*K*K*K*u_plus*u_plus*u_plus);
        return val;
    }
};


double computeWallShearStress(double u_tan, double walldist, double mu, double rho);
double computeFrictionVelocity(double u_tan, double walldist, double mu, double rho);
double computeWallVelocity(double u_tan, double walldist, double mu, double rho);


#endif
