#include "cFluid.h"
#include "cFturb.h"

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};

    
void G_CARTESIAN::computeSGSTerms()
{
    computeEddyViscosity();
}

void G_CARTESIAN::computeEddyViscosity()
{
    EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
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
            printf("\nERROR computeSGSTerms(): invalid dim\n");
            printf("\t\t dim = %d\n\n",dim);
            LOC(); clean_up(EXIT_FAILURE);
    }
}

void G_CARTESIAN::computeEddyViscosity2d()
{
    double* mu = field.mu;
    double* mu_turb = field.mu_turb;
    double* k_turb = field.k_turb;
    
    EDDY_VISC_MODEL model = eqn_params->eddy_viscosity_model;
    
    for (int j = imin[1]; j <= imax[1]; ++j)
    for (int i = imin[0]; i <= imax[0]; ++i)
    {
        int icoords[MAXD] = {i, j, 0};
        int index = d_index(icoords,top_gmax,dim);
        
        COMPONENT comp = top_comp[index];
        if (!gas_comp(comp))
        {
            mu_turb[index] = 0.0;
            k_turb[index] = 0.0;
            continue;
        }

        switch (model)
        {
            case EDDY_VISC_MODEL::VREMAN:
            {
                mu_turb[index] = computeEddyViscosityVremanModel_BdryAware(icoords);
                break;
            }
            case EDDY_VISC_MODEL::WALE:
            {
                mu_turb[index] = computeEddyViscosityWALE(icoords);
                break;
            }
            default:
            {
                printf("\nERROR: unrecognized eddy viscosity model\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
        }
        
        if (mu[index] + mu_turb[index] > mu_max)
        {
            mu_max = mu[index] + mu_turb[index];
        }
    }
}

void G_CARTESIAN::computeEddyViscosity3d()
{
    double* mu = field.mu;
    double* mu_turb = field.mu_turb;
    double* k_turb = field.k_turb;

    EDDY_VISC_MODEL model = eqn_params->eddy_viscosity_model;

    for (int k = imin[2]; k <= imax[2]; ++k)
    for (int j = imin[1]; j <= imax[1]; ++j)
    for (int i = imin[0]; i <= imax[0]; ++i)
    {
        int icoords[MAXD] = {i, j, k};
        int index = d_index(icoords,top_gmax,dim);
        
        COMPONENT comp = top_comp[index];
        if (!gas_comp(comp))
        {
            mu_turb[index] = 0.0;
            k_turb[index] = 0.0;
            continue;
        }

        switch (model)
        {
            case EDDY_VISC_MODEL::VREMAN:
            {
                mu_turb[index] = computeEddyViscosityVremanModel_BdryAware(icoords);
                break;
            }
            case EDDY_VISC_MODEL::WALE:
            {
                mu_turb[index] = computeEddyViscosityWALE(icoords);
                break;
            }
            default:
            {
                printf("\nERROR: unrecognized eddy viscosity model\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
        }
    
        if (mu[index] + mu_turb[index] > mu_max)
        {
            mu_max = mu[index] + mu_turb[index];
        }
    }
}

double G_CARTESIAN::computeEddyViscosityVremanModel(int* icoords)
{
    double **vel = field.vel;

    double C_v = eqn_params->C_v; //double C_v = 0.025;
    
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
        //TEMP
        double vel_p = vel[j][index[2*i+1]];
        if (std::isinf(vel_p) || std::isnan(vel_p))
        {
            vel_p = 0.0;
        }

        //TEMP
        double vel_m = vel[j][index[2*i]];
        if (std::isinf(vel_m) || std::isnan(vel_m))
        {
            vel_m = 0.0;
        }

        alpha[i][j] = 0.5*(vel_p - vel_m)/top_h[i];
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

    double mu_t = nu_t*field.dens[index0];

    if (std::isinf(mu_t) || std::isnan(mu_t))
    {
        printf("\nERROR: inf/nan eddy viscosity!\n");
        printf("nu_t = %g  dens[%d] = %g\n",nu_t,index0,field.dens[index0]);
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

        LOC(); clean_up(EXIT_FAILURE);
    }

    return mu_t;
}   /* end computeEddyViscosityVremanModel */

double G_CARTESIAN::computeEddyViscosityVremanModel_BdryAware(int* icoords)
{
    auto alpha = computeVelocityGradient(icoords);

    double sum_alpha = 0.0;
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        sum_alpha += alpha[i][j]*alpha[i][j];
    }
 
    double beta[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};

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

    double nu_t;
    double C_v = eqn_params->C_v; //double C_v = 0.025;

    //see Vreman's implementation, he actually uses 1.0e-12 when checking B_beta
    //  (about 10x larger than MACH_EPS)
    if (sum_alpha < MACH_EPS || B_beta < MACH_EPS)
    {
        nu_t = 0.0;
    }
    else
    {
        nu_t = C_v*sqrt(B_beta/sum_alpha);
    }

    if (std::isinf(nu_t) || std::isnan(nu_t))
    {
        printf("\nERROR: inf/nan eddy viscosity!\n");
        printf("nu_t = %g\n",nu_t);
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

        int index = d_index(icoords,top_gmax,dim);
        auto coords = cell_center[index].getCoords();
        printf("coords = ");
        for (int i = 0; i < dim; ++i)
        {
            printf("%f  ", coords[i]);
        }
        printf("\n\n");

        LOC(); clean_up(EXIT_FAILURE);
    }


    int index = d_index(icoords,top_gmax,dim);
    double mu_t = nu_t*field.dens[index];

    //Compute turbulent kinetic energy: 
    //                                             
    //        k_turb = 2.0*nu_t*|S|
    //
    //  where S = 0.5*(du_i/dx_j + du_j/dx_i)
    //  
    //  and |S| = sqrt(2.0*S_ij*S_ij) = |sqrt(2.0)*S|_{Frobenius}
    
    std::vector<std::vector<double>> S(dim,std::vector<double>(dim,0.0));
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        S[i][j] = 0.5*(alpha[i][j] + alpha[j][i]);
    }
    
    double NormFrobenius = 0.0;
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        NormFrobenius += S[i][j]*S[i][j];
    }
    NormFrobenius = std::sqrt(2.0)*std::sqrt(NormFrobenius);

    field.k_turb[index] = 2.0*mu_t*NormFrobenius;

    //TODO: Is this the correct handling of turbulent kinetic energy here?
    //      Or should pressure be isolated, and compute the flux of the TKE
    //      independently in the WENO convective flux computations???
    
    /*
        field.pres[index] += 2.0/3.0*field.k_turb[index];
    */

    //TODO: Return a std::pair<double,double> instead of opaquely
    //      adding 2/3*k_turb to the pressure as a side effect of
    //      this function 
    
    return mu_t;
}   /* end computeEddyViscosityVremanModel_BdryAware */

double G_CARTESIAN::computeEddyViscosityWALE(int* icoords)
{
    auto alpha = computeVelocityGradient(icoords);

    double sum_S = 0.0;
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        sum_S += 0.5*(alpha[i][j] + alpha[j][i]); 
    }

    double sum_Sd = 0.0;
    for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
    {
        for (int k = 0; k < dim; ++k)
        {
            sum_Sd += 0.5*alpha[i][k]*alpha[k][j];
        }
    }

    for (int k = 0; k < dim; ++k)
    {
        sum_Sd -= alpha[k][k]*alpha[k][k]/3.0;
    }

    double C_v = eqn_params->C_v; //double Cv = 0.35;
    double Delta = 1.0;
    for (int k = 0; k < dim; ++k)
    {
        Delta *= top_h[k];
    }
    
    double nu_t = sqr(C_v*Delta)*std::pow(sum_Sd,3.0/2.0) /
        (std::pow(sum_S,5.0/2.0) + std::pow(sum_Sd,5.0/4.0)); 

    int index = d_index(icoords,top_gmax,dim);
    double mu_t = nu_t*field.dens[index];

    return mu_t;
}

std::vector<std::vector<double>> 
G_CARTESIAN::computeVelocityGradient(int *icoords)
{
    std::vector<std::vector<double>> J(dim,std::vector<double>(dim,0.0));

    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];
    if (!gas_comp(comp)) return J;

    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    boolean fr_crx_grid_seg;
    INTERFACE *grid_intfc = front->grid_intfc;
    STATE* intfc_state;
    HYPER_SURF *hs;
    double crx_coords[MAXD];

    double vel_nb[2];
    double d_h[2];
    
    double** vel = field.vel;

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
                int index_nb = next_index_in_dir(icoords,dir[m][nb],top_gmax,dim);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == ELASTIC_BOUNDARY || 
                     wave_type(hs) == ELASTIC_BAND_BOUNDARY)
            {
                double v_poro[MAXD] = {0.0};
                if (eqn_params->porosity == 0 || !eqn_params->with_porosity)
                {
                    setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,vel,v_poro);
                }
                else
                {
                    setPoroSlipBoundaryNIP(icoords,m,nb,comp,hs,intfc_state,vel,v_poro);
                }
                vel_nb[nb] = v_poro[l];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                     wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                double v_slip[MAXD] = {0.0};
                setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,vel,v_slip);
                vel_nb[nb] = v_slip[l];
            }
            else if (wave_type(hs) == SYMMETRY_BOUNDARY)
            {
                double v_sym[MAXD] = {0.0};
                setSymmetryBoundaryNIP(icoords,m,nb,comp,hs,intfc_state,vel,v_sym);
                vel_nb[nb] = v_sym[l];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                /*if (boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),"cF_flowThroughBoundaryState") == 0)*/
                if (boundary_state_function(hs))
                {
                    //FLOW_THROUGH_PARAMS params;
                    //POINT *oldp = ft_params->oldp;
                    //COMPONENT comp = ft_params->comp;
                    //EQN_PARAMS *eqn_params = ft_params->eqn_params;

                    //OUTLET
                    vel_nb[nb] = intfc_state->vel[l];
                }
                else
                {
                    //INLET
                    vel_nb[nb] = intfc_state->vel[l];
                }

                if (std::isnan(vel_nb[nb]))
                    vel_nb[nb] = vel[l][index];
            }
            else
            {
                printf("ERROR: Unknown Boundary Type\n");
                printf("\n\t wave_type(hs) = %d\n",wave_type(hs));
                LOC(); clean_up(EXIT_FAILURE);
            }
        }

        J[l][m] = (vel_nb[1] - vel_nb[0])/(d_h[1] + d_h[0]);

        if (std::isnan(vel_nb[0]) || std::isinf(vel_nb[0]) ||
            std::isnan(vel_nb[1]) || std::isinf(vel_nb[1]))
        {
            printf("\ncomputeVelocityGradient() ERROR: nan/inf vel\n");
            printf("\nicoords = %d %d %d\n",icoords[0],icoords[1],icoords[2]);
            printf("vel component: %d , derivative direction: %d\n", l, m);
            printf("vel_nb[0] = %f , vel_nb[1] =%f\n",vel_nb[0],vel_nb[1]);
            printf("comp =%d\n",comp);
            printf("wave_type(hs) = %d\n",wave_type(hs));
            LOC(); clean_up(EXIT_FAILURE);
        }
    }

    return J;
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
    setSlipBoundaryNIP(icoords,idir,nb,comp,hs,state,vel,v_slip);
}

//TODO: Make function return the computed slip velocity, v_slip, as a std::vector<double>
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
    //TODO: dir is unused here
    /*
    GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };
    */

    int ghost_ic[MAXD];
    double coords[MAXD], crx_coords[MAXD];
    double coords_reflect[MAXD], coords_ghost[MAXD];
    double nor[MAXD];
    
    double vel_intfc_gcrx[MAXD];
    for (int i = 0; i < dim; ++i)
    {
        vel_intfc_gcrx[i] = (*getStateVel[i])(state);
        coords[i] = top_L[i] + icoords[i]*top_h[i];
        ghost_ic[i] = icoords[i];
    }
    
    ghost_ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
    int ghost_index = d_index(ghost_ic,top_gmax,dim);
    COMPONENT ghost_comp = top_comp[ghost_index];
    
    for (int j = 0; j < dim; ++j)
    {
        coords_ghost[j] = top_L[j] + ghost_ic[j]*top_h[j];
    }

    double intrp_coeffs[MAXD] = {0.0};
    HYPER_SURF_ELEMENT* hsurf_elem;
    HYPER_SURF* hsurf;

    

    bool nip_found = nearest_interface_point(coords_ghost,comp,front->grid_intfc,
            NO_SUBDOMAIN,nullptr,crx_coords,intrp_coeffs,&hsurf_elem,&hsurf);

    /*
    bool nip_found = nearest_interface_point(coords_ghost,comp,front->interf,
            NO_SUBDOMAIN,nullptr,crx_coords,intrp_coeffs,&hsurf_elem,&hsurf);
    */

    if (!nip_found)
    {
        printf("ERROR G_CARTESIAN::setSlipBoundaryNIP(): "
                "can't find nearest interface point\n");
        LOC(); clean_up(EXIT_FAILURE);
    }


    double dist_ghost = distance_between_positions(coords_ghost,crx_coords,dim);
    
    
    //compute the normal and velocity vectors at the interface point
    double vel_intfc[MAXD] = {0.0};
    switch (dim)
	{
        case 2:
            {
                double ns[MAXD] = {0.0};
                double ne[MAXD] = {0.0};
                
                normal(Bond_of_hse(hsurf_elem)->start,hsurf_elem,hsurf,ns,front);
                normal(Bond_of_hse(hsurf_elem)->end,hsurf_elem,hsurf,ne,front);

                STATE* ss;
                STATE* se;

                if (gas_comp(negative_component(hsurf)))
                {
                    ss = (STATE*)left_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)left_state(Bond_of_hse(hsurf_elem)->end);
                }
                else if (gas_comp(positive_component(hsurf)))
                {
                    ss = (STATE*)right_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)right_state(Bond_of_hse(hsurf_elem)->end);
                }
                else
                {
                    printf("ERROR setSlipBoundaryNIP(): "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                    vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
                }
            }
            break;

        case 3:
            {
                TRI* nearTri = Tri_of_hse(hsurf_elem);
                const double* tnor = Tri_normal(nearTri);
                
                STATE* st[3];

                if (gas_comp(negative_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)left_state(Point_of_tri(nearTri)[j]);
                }
                else if (gas_comp(positive_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)right_state(Point_of_tri(nearTri)[j]);
                }
                else
                {
                    printf("ERROR setSlipBoundaryNIP(): "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = tnor[i];

                    vel_intfc[i] = 0.0;
                    for (int j = 0; j < 3; ++j)
                        vel_intfc[i] += intrp_coeffs[j]*st[j]->vel[i];
                }
            }
            break;
	}

    //NOTE: Tri_normal() does not return a unit vector
    double mag_nor = Magd(nor,dim);
    for (int i = 0; i < dim; ++i)
        nor[i] /= mag_nor;
        
    if (comp == negative_component(hsurf))
	{
	    for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
	}
    
    //NOTE: must use unit-length vectors with FT_GridSizeInDir()
    double dist_reflect = FT_GridSizeInDir(nor,front);

    //TODO: need to check if dist_reflect > dist_ghost ???
    
        /*
        // Compute dist_reflect as the diagonal length of rect grid blocks
        double dist_reflect = 0.0;
        for (int j = 0; j < 3; ++j)
             dist_reflect += sqr(top_h[j]);
        dist_reflect = sqrt(dist_reflect);
        */

    
    //The desired reflected point
    for (int j = 0; j < dim; ++j)
    {
        coords_reflect[j] = crx_coords[j] + dist_reflect*nor[j];
    }

    
    if (debugging("slip_boundary"))
    {
        printf("\nsetSlipBoundaryNIP() DEBUGGING\n");
        printf("idir = %d nb = %d\n",idir,nb);
        fprint_int_vector(stdout,"icoords",icoords,dim,", ");
        fprint_int_vector(stdout,"ghost_ic",ghost_ic,dim,"\n");
        fprint_general_vector(stdout,"coords",coords,dim,"\n");
        fprint_general_vector(stdout,"coords_ghost",coords_ghost,dim,"\n");
        fprint_general_vector(stdout,"coords_nip",crx_coords,dim,"\n");
        fprint_general_vector(stdout,"normal",nor,dim,"\n");
        fprint_general_vector(stdout,"coords_reflect",coords_reflect,dim,"\n");
        printf("dist_ghost = %g , dist_reflect = %g\n",dist_ghost,dist_reflect);
        printf("dist_ghost/dist_reflect = %g  dist_reflect - dist_ghost = %g\n",
                dist_ghost/dist_reflect, dist_reflect - dist_ghost);
    }

    
    // Interpolate the velocity at the reflected point
    int index = d_index(icoords,top_gmax,dim);
    double vel_reflect[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,vel[j],
                getStateVel[j],&vel_reflect[j],&vel[j][index]);
    }
 
    double vel_rel[MAXD] = {0.0};
    double vn = 0.0;

    for (int j = 0; j < dim; ++j)
    {
        vel_rel[j] = vel_reflect[j] - vel_intfc[j];
        vn += vel_rel[j]*nor[j];
    }

    double vel_rel_tan[MAXD] = {0.0};
    double vel_rel_nor[MAXD] = {0.0};
    double vel_ghost_nor[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
	    vel_rel_tan[j] = vel_rel[j] - vn*nor[j];
	    vel_rel_nor[j] = vn*nor[j];
	    vel_ghost_nor[j] = -1.0*(dist_ghost/dist_reflect)*vn*nor[j];
    }
    double mag_vtan = Magd(vel_rel_tan,dim);

    if (debugging("slip_boundary"))
    {
        fprint_general_vector(stdout,"vel_reflect",vel_reflect,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_tan",vel_rel_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_nor",vel_rel_nor,dim,"\n");
        printf("Magd(vel_rel_tan,dim) = %g\n",mag_vtan);
    }


    EOS_PARAMS eos = eqn_params->eos[GAS_COMP2];
    double R_specific = eos.R_specific;
    double gamma = eos.gamma;
    double Pr = eos.Pr;

    double Cp = gamma/(gamma - 1.0)*R_specific;
    double sqrmag_vel_tan = Dotd(vel_rel_tan, vel_rel_tan, dim);

    double temp_wall = eqn_params->fixed_wall_temp;
    if (!eqn_params->use_fixed_wall_temp)
    {
        //Interpolate the temperature at the reflected point
        double temp_reflect;
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.temp,
                getStateTemp,&temp_reflect,&field.temp[index]);

        //Compute Wall Temperature
        temp_wall = temp_reflect + 0.5*pow(Pr,1.0/3.0)*sqrmag_vel_tan/Cp;
    }

    //Interpolate the pressure at the reflected point
    double pres_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.pres,
            getStatePres,&pres_reflect,&field.pres[index]);

    //Compute density near wall using the wall temperature and the pressure at the reflected point
    double dens_wall = pres_reflect/temp_wall/R_specific;

    //Interpolate the viscosity at the reflected point
    double mu_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu,
            getStateMu,&mu_reflect,&field.mu[index]);

    double mu_turb_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu_turb,
            getStateMuTurb,&mu_turb_reflect,&field.mu_turb[index]);
    
    double mu_total_reflect = mu_reflect + mu_turb_reflect;

    if (std::isnan(mu_total_reflect) || std::isinf(mu_total_reflect))
    {
       printf("\nsetSlipBoundaryNIP() ERROR: nan/inf mu_reflect\n");
       printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,field.mu[index]);
       printf("mu_turb = %g , mu_turb[%d] = %g\n",mu_turb_reflect,index,field.mu_turb[index]);
       LOC(); clean_up(EXIT_FAILURE);
    }

    double tau_wall[MAXD] = {0.0};
    double mag_tau_wall = computeWallShearStress(mag_vtan,dist_reflect,mu_reflect,dens_wall,100.0);

    if (mag_vtan > MACH_EPS)
    {
        for (int j = 0; j < dim; ++j)
            tau_wall[j] = mag_tau_wall*vel_rel_tan[j]/mag_vtan;
    }

    double vel_ghost_tan[MAXD] = {0.0};
    double vel_ghost_rel[MAXD] = {0.0};

    double coeff_tau = 0.0;
    if (mu_total_reflect > MACH_EPS)
    {
        coeff_tau = (dist_reflect - dist_ghost)/mu_total_reflect;
    }
    
    double slip = 1.0;
    if (eqn_params->no_slip_wall)
    {
        slip = 0.0;
    }

    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] = slip*vel_rel_tan[j] - coeff_tau*tau_wall[j];
        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        v_slip[j] = vel_ghost_rel[j] + vel_intfc[j];
    }


    if (std::isnan(v_slip[0]) || std::isinf(v_slip[0]) ||
        std::isnan(v_slip[1]) || std::isinf(v_slip[1]) ||
        std::isnan(v_slip[2]) || std::isinf(v_slip[2]))
    {
        printf("\nsetSlipBoundaryNIP() ERROR: nan/inf v_slip\n");
        printf("\nidir = %d nb = %d\n",idir,nb);
        fprint_int_vector(stdout,"icoords",icoords,dim,"\n");
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        printf("mu_reflect = %g , mu_[%d] = %g\n",mu_reflect,index,field.mu[index]);
        LOC(); clean_up(EXIT_FAILURE);
    }

    if (debugging("slip_boundary"))
    {
        printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,field.mu[index]);
        printf("mag_tau_wall = %g\n",mag_tau_wall);
        fprint_general_vector(stdout,"tau_wall",tau_wall,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_rel",vel_ghost_rel,dim,"\n");
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
    }
    
    /*
    //store data to avoid recomputing values in the fluid solver -- see iFluid handling of this
    int fid = face_index(idir,nb);
    for (int i = 0; i < dim; ++i)
    {
        ghost_data[fid][index].vel[i] = v_slip[i];
        ghost_data[fid][index].force[i] = tau_wall[i];
    }
    */
}

//TODO: THIS SHOULD BE A MEMBER FUNCTION OF CFABRIC_CARTESIAN
void G_CARTESIAN::setPoroSlipBoundaryNIP(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double** vel,
	double* v_slip)
{
    /*
    GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };
    */

    int ghost_ic[MAXD];
    double coords[MAXD], crx_coords[MAXD];
    double coords_reflect[MAXD], coords_ghost[MAXD];
    double nor[MAXD];
    
    double vel_intfc_gcrx[MAXD];
    for (int i = 0; i < dim; ++i)
    {
        vel_intfc_gcrx[i] = (*getStateVel[i])(state);
        coords[i] = top_L[i] + icoords[i]*top_h[i];
        ghost_ic[i] = icoords[i];
    }
    
    ghost_ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
    int ghost_index = d_index(ghost_ic,top_gmax,dim);
    COMPONENT ghost_comp = top_comp[ghost_index];
    
    for (int j = 0; j < dim; ++j)
    {
        coords_ghost[j] = top_L[j] + ghost_ic[j]*top_h[j];
    }

    /*
    double real_dens = field.dens[ghost_index];
    double real_pres = field.pres[ghost_index];
    double real_vel[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
    {
        real_vel[j] = field.vel[j][ghost_index];
    }
    */

    double intrp_coeffs[MAXD] = {0.0};
    HYPER_SURF_ELEMENT* hsurf_elem;
    HYPER_SURF* hsurf;
    int range = 2;

    //TODO: Why does this fail for INCLUDE_BOUNDARIES and NO_SUBDOMAIN values?
    //      Conversely, why does it work with NO_BOUNDARIES in the backward facing
    //      step scenario -- to what degree is it working?
    
    /*
    bool nip_found = nearest_interface_point(coords_ghost,comp,front->interf,
            NO_SUBDOMAIN,hs,crx_coords,intrp_coeffs,&hsurf_elem,&hsurf);
    */
    
    //bool nip_found = nearest_interface_point(coords_ghost,comp,front->interf,
      //      NO_SUBDOMAIN,nullptr,crx_coords,intrp_coeffs,&hsurf_elem,&hsurf);
    
    bool nip_found = nearest_interface_point(coords_ghost,comp,front->grid_intfc,
            NO_SUBDOMAIN,nullptr,crx_coords,intrp_coeffs,&hsurf_elem,&hsurf);

    if (!nip_found)
    {
        printf("ERROR G_CARTESIAN::setSlipBoundaryNIP(): "
                "can't find nearest interface point\n");
        printf("comp = %d\n", comp);
        printf("ghost_comp = %d\n", ghost_comp);
        printf("positive_component(hs) = %d\n",positive_component(hs));
        printf("negative_component(hs) = %d\n",negative_component(hs));
        printf("crx_coords = %f %f %f\n",crx_coords[0],crx_coords[1],crx_coords[2]);
        LOC(); clean_up(EXIT_FAILURE);
    }


    double dist_ghost = distance_between_positions(coords_ghost,crx_coords,dim);
    
    //compute the normal and velocity vectors at the interface point
    double vel_intfc[MAXD] = {0.0};
    switch (dim)
	{
        case 2:
            {
                double ns[MAXD] = {0.0};
                double ne[MAXD] = {0.0};
                
                normal(Bond_of_hse(hsurf_elem)->start,hsurf_elem,hsurf,ns,front);
                normal(Bond_of_hse(hsurf_elem)->end,hsurf_elem,hsurf,ne,front);

                STATE* ss;
                STATE* se;

                if (gas_comp(negative_component(hsurf)))
                {
                    ss = (STATE*)left_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)left_state(Bond_of_hse(hsurf_elem)->end);
                }
                else if (gas_comp(positive_component(hsurf)))
                {
                    ss = (STATE*)right_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)right_state(Bond_of_hse(hsurf_elem)->end);
                }
                else
                {
                    printf("ERROR setSlipBoundaryNIP(): "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                    vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
                }
            }
            break;

        case 3:
            {
                TRI* nearTri = Tri_of_hse(hsurf_elem);
                const double* tnor = Tri_normal(nearTri);
                
                STATE* st[3];

                if (gas_comp(negative_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)left_state(Point_of_tri(nearTri)[j]);
                }
                else if (gas_comp(positive_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)right_state(Point_of_tri(nearTri)[j]);
                }
                else
                {
                    printf("ERROR setSlipBoundaryNIP(): "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = tnor[i];

                    vel_intfc[i] = 0.0;
                    for (int j = 0; j < 3; ++j)
                        vel_intfc[i] += intrp_coeffs[j]*st[j]->vel[i];
                }
            }
            break;
	}

    //NOTE: Tri_normal() does not return a unit vector
    double mag_nor = Magd(nor,dim);
    for (int i = 0; i < dim; ++i)
        nor[i] /= mag_nor;
        
    double vec_ghost[MAXD] = {0,0};
    for (int j = 0; j < dim; ++j)
    {
        vec_ghost[j] = coords_ghost[j] - crx_coords[j];
    }

    double side = Dotd(vec_ghost,nor,dim);
    if (side < 0)
	{
	    for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
	}

    //NOTE: must use unit-length vectors with FT_GridSizeInDir()
    double dist_reflect = FT_GridSizeInDir(nor,front);
    
    //The desired reflected point
    for (int j = 0; j < dim; ++j)
    {
        coords_reflect[j] = crx_coords[j] - dist_reflect*nor[j];
    }
    //NOTE: In this function nor points in the opposite direction 
    //      of the nor vector used in setSlipBoundaryNIP()

    
    if (debugging("poro_slip_boundary"))
    {
        printf("\nsetPoroSlipBoundaryNIP() DEBUGGING\n");
        printf("idir = %d nb = %d\n",idir,nb);
        fprint_int_vector(stdout,"icoords",icoords,dim,", ");
        fprint_int_vector(stdout,"ghost_ic",ghost_ic,dim,"\n");
        fprint_general_vector(stdout,"coords",coords,dim,"\n");
        fprint_general_vector(stdout,"coords_ghost",coords_ghost,dim,"\n");
        fprint_general_vector(stdout,"coords_nip",crx_coords,dim,"\n");
        fprint_general_vector(stdout,"normal",nor,dim,"\n");
        fprint_general_vector(stdout,"coords_reflect",coords_reflect,dim,"\n");
        printf("dist_ghost = %g , dist_reflect = %g\n",dist_ghost,dist_reflect);
        printf("dist_ghost/dist_reflect = %g  dist_reflect - dist_ghost = %g\n",
                dist_ghost/dist_reflect, dist_reflect - dist_ghost);
    }

    
    // Interpolate the velocity at the reflected point
    int index = d_index(icoords,top_gmax,dim);
    double vel_reflect[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,vel[j],
                getStateVel[j],&vel_reflect[j],&vel[j][index]);
    }
 
    double vel_rel[MAXD] = {0.0};
    double vn = 0.0;
    for (int j = 0; j < dim; ++j)
    {
        vel_rel[j] = vel_reflect[j] - vel_intfc[j];
        vn += vel_rel[j]*nor[j];
    }

    double vel_rel_tan[MAXD] = {0.0};
    double vel_rel_nor[MAXD] = {0.0};
    double vel_ghost_nor[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
	    vel_rel_tan[j] = vel_rel[j] - vn*nor[j];
	    vel_rel_nor[j] = vn*nor[j];
        vel_ghost_nor[j] = -1.0*(dist_ghost/dist_reflect)*vn*nor[j];
    }
    double mag_vtan = Magd(vel_rel_tan,dim);

    if (debugging("poro_slip_boundary"))
    {
        fprint_general_vector(stdout,"vel_reflect",vel_reflect,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_tan",vel_rel_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_nor",vel_rel_nor,dim,"\n");
        printf("Magd(vel_rel_tan,dim) = %g\n",mag_vtan);
    }


    EOS_PARAMS eos = eqn_params->eos[GAS_COMP2];
    double R_specific = eos.R_specific;
    double gamma = eos.gamma;
    double Pr = eos.Pr;

    double Cp = gamma/(gamma - 1.0)*R_specific;
    double sqrmag_vel_tan = Dotd(vel_rel_tan, vel_rel_tan, dim);

    double temp_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.temp,
            getStateTemp,&temp_reflect,&field.temp[index]);
    
    //Compute Wall Temperature
    double temp_wall = temp_reflect + 0.5*pow(Pr,1.0/3.0)*sqrmag_vel_tan/Cp;

    //Interpolate the pressure at the reflected point
    double pres_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.pres,
            getStatePres,&pres_reflect,&field.pres[index]);

    double mu_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu,
            getStateMu,&mu_reflect,&field.mu[index]);

    double mu_turb_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu_turb,
            getStateMuTurb,&mu_turb_reflect,&field.mu_turb[index]);
    
    double mu_total_reflect = mu_reflect + mu_turb_reflect;

    //Compute density near wall using the wall temperature and the pressure at the reflected point
    double dens_wall = pres_reflect/temp_wall/R_specific;


    if (std::isnan(mu_total_reflect) || std::isinf(mu_total_reflect))
    {
       printf("\nsetPoroSlipBoundaryNIP() ERROR: nan/inf mu_reflect\n");
       printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,field.mu[index]);
       printf("mu_turb_reflect = %g , mu_turb[%d] = %g\n",
               mu_turb_reflect,index,field.mu_turb[index]);
       LOC(); clean_up(EXIT_FAILURE);
    }

    double poro = eqn_params->porosity;
    double perm = eqn_params->permeability;
    double Cw = 1.0; //TODO: read into eqn_params from input file

    //Use modified distance in wall law for porous interface
    double delta_y = Cw*std::sqrt(perm/poro);
    double dist_wall = dist_reflect + delta_y;

    double tau_wall[MAXD] = {0.0};
    double mag_tau_wall = computeWallShearStress(mag_vtan,dist_wall,mu_reflect,dens_wall,100.0);

    if (mag_vtan > MACH_EPS)
    {
        for (int j = 0; j < dim; ++j)
            tau_wall[j] = mag_tau_wall*vel_rel_tan[j]/mag_vtan;
    }

    double vel_ghost_tan[MAXD] = {0.0};
    double vel_ghost_rel[MAXD] = {0.0};

    double coeff_tau = 0.0;
    if (mu_total_reflect > MACH_EPS)
    {
        coeff_tau = (dist_wall - dist_ghost)/mu_total_reflect;
            //coeff_tau = (dist_reflect - dist_ghost)/mu_total_reflect;
    }
    
    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] = vel_rel_tan[j] - coeff_tau*tau_wall[j];
        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        v_slip[j] = vel_ghost_rel[j] + vel_intfc[j];
    }


    if (std::isnan(v_slip[0]) || std::isinf(v_slip[0]) ||
        std::isnan(v_slip[1]) || std::isinf(v_slip[1]) ||
        std::isnan(v_slip[2]) || std::isinf(v_slip[2]))
    {
        printf("\nsetPoroSlipBoundaryNIP() ERROR: nan/inf v_slip\n");
        printf("\nidir = %d nb = %d\n",idir,nb);
        fprint_int_vector(stdout,"icoords",icoords,dim,"\n");
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,field.mu[index]);
        LOC(); clean_up(EXIT_FAILURE);
    }

    if (debugging("poro_slip_boundary"))
    {
        printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,field.mu[index]);
        printf("mag_tau_wall = %g\n",mag_tau_wall);
        fprint_general_vector(stdout,"tau_wall",tau_wall,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_rel",vel_ghost_rel,dim,"\n");
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
    }
    
    /*
    //store data to avoid recomputing values in the fluid solver -- see iFluid handling of this
    int fid = face_index(idir,nb);
    for (int i = 0; i < dim; ++i)
    {
        ghost_data[fid][index].vel[i] = v_slip[i];
        ghost_data[fid][index].force[i] = tau_wall[i];
    }
    */
}

void G_CARTESIAN::setSymmetryBoundaryNIP(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double** vel,
	double* v_sym)
{
    int ghost_ic[MAXD];
    double coords[MAXD], crx_coords[MAXD];
    double coords_reflect[MAXD], coords_ghost[MAXD];
    double nor[MAXD];
    
    double vel_intfc_gcrx[MAXD];
    for (int i = 0; i < dim; ++i)
    {
        vel_intfc_gcrx[i] = (*getStateVel[i])(state);
        coords[i] = top_L[i] + icoords[i]*top_h[i];
        ghost_ic[i] = icoords[i];
    }
    
    ghost_ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
    int ghost_index = d_index(ghost_ic,top_gmax,dim);
    COMPONENT ghost_comp = top_comp[ghost_index];
    
    for (int j = 0; j < dim; ++j)
    {
        coords_ghost[j] = top_L[j] + ghost_ic[j]*top_h[j];
    }

    double intrp_coeffs[MAXD] = {0.0};
    HYPER_SURF_ELEMENT* hsurf_elem;
    HYPER_SURF* hsurf;
    int range = 2;

    //TODO: Why does this fail for INCLUDE_BOUNDARIES and NO_SUBDOMAIN values?
    //      Conversely, why does it work with NO_BOUNDARIES in the backward facing
    //      step scenario -- to what degree is it working?
    //
    //      nearest_intfc_point_in_range() etc. doesn't find NEUMANN_BOUNDARY
    //      when it is domain (rect) boundary, but nearest_intfc_point() does.
    //
    //TODO: Is the above todo stil valid? Or did I just not understand what was happening?
    //
    
    /*
    FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,NO_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */
    /*
    bool nip_found = 
        FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,INCLUDE_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */
    /*
    bool nip_found = 
        FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,INCLUDE_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */
    /*
    bool nip_found = 
        FT_FindNearestIntfcPointInRange(front,comp,coords_ghost,NO_SUBDOMAIN,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */
    /*
    bool nip_found = 
        FT_FindNearestIntfcPointInRange(front, comp, coords_ghost,INCLUDE_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */

    //FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,NO_SUBDOMAIN,
      //      crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);

    //bool nip_found = nearest_interface_point(coords_ghost,comp,front->interf,
      //      INCLUDE_BOUNDARIES,nullptr,crx_coords,intrp_coeffs,&hsurf_elem,&hsurf);

    bool nip_found = nearest_interface_point(coords_ghost,comp,front->interf,
            NO_SUBDOMAIN,nullptr,crx_coords,intrp_coeffs,&hsurf_elem,&hsurf);

    if (!nip_found)
    {
        printf("ERROR G_CARTESIAN::setSlipBoundaryNIP(): "
                "can't find nearest interface point\n");
        LOC(); clean_up(EXIT_FAILURE);
    }



    double dist_ghost = distance_between_positions(coords_ghost,crx_coords,dim);
    
    //compute the normal and velocity vectors at the interface point
    double vel_intfc[MAXD] = {0.0};
    switch (dim)
	{
        case 2:
            {
                double ns[MAXD] = {0.0};
                double ne[MAXD] = {0.0};
                
                normal(Bond_of_hse(hsurf_elem)->start,hsurf_elem,hsurf,ns,front);
                normal(Bond_of_hse(hsurf_elem)->end,hsurf_elem,hsurf,ne,front);

                STATE* ss;
                STATE* se;

                if (gas_comp(negative_component(hsurf)))
                {
                    ss = (STATE*)left_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)left_state(Bond_of_hse(hsurf_elem)->end);
                }
                else if (gas_comp(positive_component(hsurf)))
                {
                    ss = (STATE*)right_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)right_state(Bond_of_hse(hsurf_elem)->end);
                }
                else
                {
                    printf("ERROR setSlipBoundaryNIP(): "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                    vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
                }
            }
            break;

        case 3:
            {
                TRI* nearTri = Tri_of_hse(hsurf_elem);
                const double* tnor = Tri_normal(nearTri);
                
                STATE* st[3];

                if (gas_comp(negative_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)left_state(Point_of_tri(nearTri)[j]);
                }
                else if (gas_comp(positive_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)right_state(Point_of_tri(nearTri)[j]);
                }
                else
                {
                    printf("ERROR setSlipBoundaryNIP(): "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = tnor[i];

                    vel_intfc[i] = 0.0;
                    for (int j = 0; j < 3; ++j)
                        vel_intfc[i] += intrp_coeffs[j]*st[j]->vel[i];
                }
            }
            break;
	}

    //NOTE: Tri_normal() does not return a unit vector
    double mag_nor = Magd(nor,dim);
    for (int i = 0; i < dim; ++i)
        nor[i] /= mag_nor;
        
    if (comp == negative_component(hsurf))
	{
	    for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
	}
    
    //NOTE: must use unit-length vectors with FT_GridSizeInDir()
        //double dist_reflect = FT_GridSizeInDir(nor,front);
    double dist_reflect = dist_ghost;

    //The desired reflected point
    for (int j = 0; j < dim; ++j)
    {
        coords_reflect[j] = crx_coords[j] + dist_reflect*nor[j];
    }

    
    if (debugging("slip_boundary"))
    {
        printf("\nsetSlipBoundaryNIP() DEBUGGING\n");
        printf("idir = %d nb = %d\n",idir,nb);
        fprint_int_vector(stdout,"icoords",icoords,dim,", ");
        fprint_int_vector(stdout,"ghost_ic",ghost_ic,dim,"\n");
        fprint_general_vector(stdout,"coords",coords,dim,"\n");
        fprint_general_vector(stdout,"coords_ghost",coords_ghost,dim,"\n");
        fprint_general_vector(stdout,"coords_nip",crx_coords,dim,"\n");
        fprint_general_vector(stdout,"normal",nor,dim,"\n");
        fprint_general_vector(stdout,"coords_reflect",coords_reflect,dim,"\n");
        printf("dist_ghost = %g , dist_reflect = %g\n",dist_ghost,dist_reflect);
        printf("dist_ghost/dist_reflect = %g  dist_reflect - dist_ghost = %g\n",
                dist_ghost/dist_reflect, dist_reflect - dist_ghost);
    }

    
    // Interpolate the velocity at the reflected point
    int index = d_index(icoords,top_gmax,dim);
    double vel_reflect[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,vel[j],
                getStateVel[j],&vel_reflect[j],&vel[j][index]);
    }
 
    double vel_rel[MAXD] = {0.0};
    double vn = 0.0;

    for (int j = 0; j < dim; ++j)
    {
        vel_rel[j] = vel_reflect[j] - vel_intfc[j];
        vn += vel_rel[j]*nor[j];
    }

    double vel_rel_tan[MAXD] = {0.0};
    double vel_rel_nor[MAXD] = {0.0};
    double vel_ghost_nor[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
	    vel_rel_tan[j] = vel_rel[j] - vn*nor[j];
	    vel_rel_nor[j] = vn*nor[j];
	    vel_ghost_nor[j] = -1.0*(dist_ghost/dist_reflect)*vn*nor[j];
    }

    if (debugging("slip_boundary"))
    {
        fprint_general_vector(stdout,"vel_reflect",vel_reflect,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_tan",vel_rel_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_nor",vel_rel_nor,dim,"\n");
    }

    double vel_ghost_tan[MAXD] = {0.0};
    double vel_ghost_rel[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] = vel_rel_tan[j];
        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        v_sym[j] = vel_ghost_rel[j] + vel_intfc[j];
    }


    if (std::isnan(v_sym[0]) || std::isinf(v_sym[0]) ||
        std::isnan(v_sym[1]) || std::isinf(v_sym[1]) ||
        std::isnan(v_sym[2]) || std::isinf(v_sym[2]))
    {
        printf("\nsetsymBoundaryNIP() ERROR: nan/inf v_sym\n");
        printf("\nidir = %d nb = %d\n",idir,nb);
        fprint_int_vector(stdout,"icoords",icoords,dim,"\n");
        fprint_general_vector(stdout,"v_sym",v_sym,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        LOC(); //clean_up(EXIT_FAILURE);
    }

    if (debugging("sym_boundary"))
    {
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_rel",vel_ghost_rel,dim,"\n");
        fprint_general_vector(stdout,"v_sym",v_sym,dim,"\n");
    }
    
    /*
    //store data to avoid recomputing values in the fluid solver -- see iFluid handling of this
    int fid = face_index(idir,nb);
    for (int i = 0; i < dim; ++i)
    {
        ghost_data[fid][index].vel[i] = v_sym[i];
        ghost_data[fid][index].force[i] = 0.0;
    }
    */
}

double computeWallShearStress(
        double u_tan,
        double walldist,
        double mu,
        double rho,
        double u_wall_initial_guess)
{
    if (u_tan < MACH_EPS || mu < MACH_EPS) return 0.0;
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

