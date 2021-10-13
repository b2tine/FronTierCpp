#include "cFluid.h"
#include "cFturb.h"

    
//TODO: Other SGS terms.
//      -- Turbulent heat flux,
//      -- Turbulent kinetic energy transport and diffusion,
//      -- etc.
//
void G_CARTESIAN::computeSGSTerms()
{
    if (!eqn_params->use_eddy_viscosity) return;

    switch (dim)
    {
        case 2:
            computeEddyViscosity2d();
            break;
        case 3:
            computeEddyViscosity3d();
            break;
        default:
            printf("\nERROR computeSGSTerms(): dim must be equal to 2 or 3\n");
            printf("\t\t dim = %d\n\n",dim);
            LOC(); clean_up(EXIT_FAILURE);
    }
}

void G_CARTESIAN::computeEddyViscosity2d()
{
    double* mu = field.mu;

    for (int j = imin[1]; j <= imax[1]; ++j)
    for (int i = imin[0]; i <= imax[0]; ++i)
    {
        int icoords[MAXD] = {i, j, 0};
        int index = d_index(icoords,top_gmax,dim);
        mu[index] = 0.0;
        
        COMPONENT comp = top_comp[index];
        if (!gas_comp(comp)) continue;

        switch (comp)
        {
            case GAS_COMP1:
                mu[index] = eqn_params->mu1;
                break;
            case GAS_COMP2:
                mu[index] = eqn_params->mu2;
                break;
            default:
                printf("\nERROR computeEddyViscosity2d(): unrecognized component!\n");
                printf("\t\tcomp = %d\n", comp);
                LOC(); clean_up(EXIT_FAILURE);
        }

        mu[index] += computeEddyViscosityVremanModel(icoords);
            //mu[index] += computeEddyViscosityVremanModel_BdryAware(icoords);
    }
}

void G_CARTESIAN::computeEddyViscosity3d()
{
    double* mu = field.mu;

    for (int k = imin[2]; k <= imax[2]; ++k)
    for (int j = imin[1]; j <= imax[1]; ++j)
    for (int i = imin[0]; i <= imax[0]; ++i)
    {
        int icoords[MAXD] = {i, j, k};
        int index = d_index(icoords,top_gmax,dim);
        mu[index] = 0.0;
        
        COMPONENT comp = top_comp[index];
        if (!gas_comp(comp)) continue;

        switch (comp)
        {
            case GAS_COMP1:
                mu[index] = eqn_params->mu1;
                break;
            case GAS_COMP2:
                mu[index] = eqn_params->mu2;
                break;
            default:
                printf("\nERROR computeEddyViscosity3d(): unrecognized component!\n");
                printf("\t\tcomp = %d\n", comp);
                LOC(); clean_up(EXIT_FAILURE);
        }

        mu[index] += computeEddyViscosityVremanModel(icoords);
            //mu[index] += computeEddyViscosityVremanModel_BdryAware(icoords);
    }
}

double G_CARTESIAN::computeEddyViscosityVremanModel(int* icoords)
{
    //printf("\nERROR computeEddyViscosityVremanModel(): function not implemented yet\n");
    //LOC(); clean_up(EXIT_FAILURE);

    double **vel = field.vel;

    double C_v = 0.025;
        //double C_v = eqn_params->C_v;
    
    int index[6], index0;
    double alpha[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};
    double beta[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};

    //Basic version without boundary awareness
    switch (dim)
    {
        case 2:
            index0 = d_index2d(icoords[0],icoords[1],top_gmax);
            index[0] = d_index2d(icoords[0]-1,icoords[1],top_gmax);
            index[1] = d_index2d(icoords[0]+1,icoords[1],top_gmax);
            index[2] = d_index2d(icoords[0],icoords[1]-1,top_gmax);
            index[3] = d_index2d(icoords[0],icoords[1]+1,top_gmax);
            break;

        case 3:
            index0 = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
            index[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
            index[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            index[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
            index[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            index[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
            index[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            break;
    }

    COMPONENT comp = top_comp[index0];
    if (!gas_comp(comp)) return 0.0;

    double sum_alpha = 0.0;
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        /*
        //TEMP
        double vel_p = vel[j][index[2*i+1]];
        if (std::isinf(vel_p) || std::isnan(vel_p))
            vel_p = 0.0;
        
        //TEMP
        double vel_m = vel[j][index[2*i]];
        if (std::isinf(vel_m) || std::isnan(vel_m))
            vel_m = 0.0;
        
        alpha[i][j] = 0.5*(vel_p - vel_m)/top_h[i];
        */

        //TEMP
        if (std::isinf(vel[j][index[2*i+1]]) || std::isnan(vel[j][index[2*i+1]]))
        {
            vel[j][index[2*i+1]] = 0.0;
        }
        
        //TEMP
        if (std::isinf(vel[j][index[2*i]]) || std::isnan(vel[j][index[2*i]]))
        {
            vel[j][index[2*i]] = 0.0;
        }

        alpha[i][j] = 0.5*(vel[j][index[2*i+1]] - vel[j][index[2*i]])/top_h[i];
        sum_alpha += alpha[i][j]*alpha[i][j];
    }

    //TODO: Implement boundary aware computeVelocityGradient(), which includes
    //      setting slip velocity etc.
    //      Will likely be a little more complicated than iFluid version, potentially
    //      requiring the state velocity and momentum values to be set in accord with
    //      computed slip velocity.... 
        
        /*
        auto alpha = computeVelocityGradient(icoords);

        double sum_alpha = 0.0;
        for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
        {
            sum_alpha += alpha[i][j]*alpha[i][j];
        }
        */
 
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        beta[i][j] = 0.0;
        for (int k = 0; k < dim; ++k)
        {
            beta[i][j] += top_h[k]*top_h[k]*alpha[k][i]*alpha[k][j];
        }
    }

    double B_beta = beta[0][0]*beta[1][1] - beta[0][1]*beta[0][1]
                  + beta[0][0]*beta[2][2] - beta[0][2]*beta[0][2]
                  + beta[1][1]*beta[2][2] - beta[1][2]*beta[1][2];

    //see Vreman's implementation, he actually uses 1.0e-12 when checking B_beta
    //  (about 10x larger than MACH_EPS)
    double nu_t;
    if (sum_alpha < MACH_EPS || B_beta < MACH_EPS)
    {
        nu_t = 0.0;
    }
    else
    {
        nu_t = C_v*sqrt(B_beta/sum_alpha);
    }

    double mu_t = nu_t*field.dens[index0];
        //double mu_t = nu_t*field.dens[index];

    if (std::isinf(mu_t) || std::isnan(mu_t))
    {
        printf("\nERROR: inf/nan eddy viscosity!\n");
        printf("nu_t = %g  dens[%d] = %g\n",nu_t,index0,field.dens[index0]);
            //printf("nu_t = %g  dens[%d] = %g\n",nu_t,index,field.dens[index]);
        printf("B_beta = %g  sum_alpha = %g\n",B_beta,sum_alpha);

        printf("\n");
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                printf("alpha[%d][%d] = %g   ",i,j,alpha[i][j]);
            }
            printf("\n\n");
        }

        for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
        {
            printf("vel[%d][%d] = %g    vel[%d][%d] = %g\n\n",
                    j,index[2*i+1],vel[j][index[2*i+1]],
                    j,index[2*i],vel[j][index[2*i]]);
        }

        auto coords0 = cell_center[index0].getCoords();
        printf("coords0 = %f %f\n\n",coords0[0],coords0[1]);
        //printf("coords0 = %f %f %f\n\n",coords0[0],coords0[1],coords0[2]);

        LOC(); clean_up(EXIT_FAILURE);
    }

    return mu_t;
}   /* end computeEddyViscosityVremanModel */

//TODO: Implement boundary aware version, which includes setting slip velocity etc.
double G_CARTESIAN::computeEddyViscosityVremanModel_BdryAware(int* icoords)
{
    //printf("\nERROR computeEddyViscosityVremanModel(): function not implemented yet\n");
    //LOC(); clean_up(EXIT_FAILURE);

    double **vel = field.vel;

    double beta[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};

    //TODO:    //double C_v = eqn_params->C_v;
    double C_v = 0.025;
    

    //TODO: Implement computeVelocityGradient()
    auto alpha = computeVelocityGradient(icoords);

    double sum_alpha = 0.0;
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        sum_alpha += alpha[i][j]*alpha[i][j];
    }
 
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        beta[i][j] = 0.0;
        for (int k = 0; k < dim; ++k)
        {
            beta[i][j] += top_h[k]*top_h[k]*alpha[k][i]*alpha[k][j];
        }
    }

    double B_beta = beta[0][0]*beta[1][1] - beta[0][1]*beta[0][1]
                  + beta[0][0]*beta[2][2] - beta[0][2]*beta[0][2]
                  + beta[1][1]*beta[2][2] - beta[1][2]*beta[1][2];

    //see Vreman's implementation, he actually uses 1.0e-12 when checking B_beta
    //  (about 10x larger than MACH_EPS)
    double nu_t;
    if (sum_alpha < MACH_EPS || B_beta < MACH_EPS)
    {
        nu_t = 0.0;
    }
    else
    {
        nu_t = C_v*sqrt(B_beta/sum_alpha);
    }

    int index = d_index(icoords,top_gmax,dim);
    double mu_t = nu_t*field.dens[index];

    if (std::isinf(mu_t) || std::isnan(mu_t))
    {
        printf("\nERROR: inf/nan eddy viscosity!\n");
        printf("nu_t = %g  dens[%d] = %g\n",nu_t,index,field.dens[index]);
        printf("B_beta = %g  sum_alpha = %g\n",B_beta,sum_alpha);

        printf("\n");
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                printf("alpha[%d][%d] = %g   ",i,j,alpha[i][j]);
            }
            printf("\n\n");
        }

        auto coords = cell_center[index].getCoords();
        printf("coords = ");
        for (int i = 0; i < dim; ++i)
        {
            printf("%f  ", coords[i]);
        }
        printf("\n\n");

        LOC(); clean_up(EXIT_FAILURE);
    }

    return mu_t;
}   /* end computeEddyViscosityVremanModel_BdryAware */

//TODO: ADAPT FOR USE WITH G_CARTESIAN
std::vector<std::vector<double>>
G_CARTESIAN::computeVelocityGradient(
	int *icoords)
{
    ///////////////////////////////////////////////////////////////////////////////
    printf("\nERROR computeVelocityGradient(): function not implemented yet\n");
    LOC(); clean_up(EXIT_FAILURE);
    ///////////////////////////////////////////////////////////////////////////////

    double** vel = field.vel;

    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    boolean fr_crx_grid_seg;
    INTERFACE *grid_intfc = front->grid_intfc;
    STATE* intfc_state;
    HYPER_SURF *hs;
    double crx_coords[MAXD];

    std::vector<std::vector<double>> J(dim,std::vector<double>(dim,0.0));

    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];
    if (!gas_comp(comp)) return J;

    double vel_nb[2];
    double d_h[2];

    /*
    for (int l = 0; l < dim; ++l)
    for (int m = 0; m < dim; ++m)
    {
        //l derivatives in m direction
        for (int nb = 0; nb < 2; ++nb)
        {
            d_h[nb] = top_h[m];

            fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    grid_intfc,icoords,dir[m][nb],comp,
                    (POINTER*)&intfc_state,&hs,crx_coords);

            if (!fr_crx_grid_seg)
            {
                int index_nb = next_index_in_dir(icoords,dir[m][nb],dim,top_gmax);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == ELASTIC_BOUNDARY)
            {
                //NOTE: zeroing/reducing the relative fluid velocity normal to
                //      canopy via the slip boundary condition causes the pressure
                //      jump imposed by the ghost fluid method porosity model
                //      to vanish -- destroying the porosity of the canopy
                
                //TODO: What is the correct approach??
                int index_nb = next_index_in_dir(icoords,dir[m][nb],dim,top_gmax);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                     wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                if (iFparams->use_no_slip)
                {
                    vel_nb[nb] = intfc_state->vel[l];
                }
                else
                {
                    double v_slip[MAXD] = {0.0};
                    setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,field->vel,v_slip);
                    vel_nb[nb] = v_slip[l];
                }
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                //TODO: How are these handled in cFluid?
                if (boundary_state_function_name(hs) &&
                        strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                {
                    //OUTLET
                    vel_nb[nb] = intfc_state->vel[l];
                }
                else
                {
                    //INLET
                    vel_nb[nb] = intfc_state->vel[l];
                }
            }
            else
            {
                printf("ERROR: Unknown Boundary Type\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
        }

        J[l][m] = (vel_nb[1] - vel_nb[0])/(d_h[1] + d_h[0]);
    }

    return J;
    */
}

void G_CARTESIAN::setSlipBoundary(
        int *icoords,
        int idir,
        int nb,
        int comp,
        HYPER_SURF *hs,
        POINTER state,
        double** vel,
        double* v_slip)
{
    //printf("\nERROR setSlipBoundary(): function not implemented yet\n");
    //LOC(); clean_up(EXIT_FAILURE);

    setSlipBoundaryNIP(icoords,idir,nb,comp,hs,state,vel,v_slip);

    //TODO: Write GNOR implementation and compare results.
    //
    //      setSlipBoundaryGNOR(icoords,idir,nb,comp,hs,state,vel,v_slip);
}

// Based on finding the nearest interface point to the ghost point
// computed using FT_FindNearestIntfcPointInRange()
void G_CARTESIAN::setSlipBoundaryNIP(
    int *icoords,
    int idir,
    int nb,
    int comp, 
    HYPER_SURF *hs,
    POINTER state, 
    double** vel,
    double* v_slip)
{ 
    printf("\nERROR setSlipBoundaryNIP(): function not implemented yet\n");
    LOC(); clean_up(EXIT_FAILURE);

    //TODO: Implemenent
}

double computeWallShearStress(
        double u_tan,
        double walldist,
        double mu,
        double rho,
        double u_wall_initial_guess)
{
    if (u_tan < MACH_EPS) return 0.0;
    double u_friction = computeFrictionVelocity(u_tan,walldist,mu/rho,
            u_wall_initial_guess);
    double tau_wall = u_friction*u_friction*rho;
    return tau_wall;
}

double computeFrictionVelocity(
        double u_tan,
        double walldist,
        double nu,
        double u_wall_initial_guess)
{
    double u_wall = computeWallVelocity(u_tan,walldist,nu,
            u_wall_initial_guess);
    if (u_wall < MACH_EPS) return 0.0;
    double u_friction = u_tan/u_wall;
    return u_friction;

    /*
    //TODO: put into data/debugging output
    //
    double y_plus = walldist*rho/mu*u_friction;
    std::cout << "y+ = " << y_plus << "\n";
    */
}

//Computes the dimensionless wall velocity u^{+} = u_tan/u_friction
double computeWallVelocity(
        double u_tan,
        double walldist,
        double nu,
        double u_wall_initial_guess)
{
    SpaldingWallLaw wallfunc(u_tan,walldist,nu);
    //TODO: method for better initial guess?
    double u_wall = wallfunc.solve(u_wall_initial_guess);
    return u_wall;
}

