#include "cFluid.h"
#include "cFturb.h"

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};

    
//TODO: Add other SGS terms.
void G_CARTESIAN::computeSGSTerms()
{
    computeEddyViscosity();
}

void G_CARTESIAN::computeEddyViscosity()
{
    EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
    mu_max = std::max(eqn_params->mu1,eqn_params->mu2);

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

    for (int j = imin[1]; j <= imax[1]; ++j)
    for (int i = imin[0]; i <= imax[0]; ++i)
    {
        int icoords[MAXD] = {i, j, 0};
        int index = d_index(icoords,top_gmax,dim);
        
        COMPONENT comp = top_comp[index];
        if (!gas_comp(comp))
        {
            mu[index] = 0.0;
            continue;
        }

        double mu_molecular;
        switch (comp)
        {
            case GAS_COMP1:
                mu_molecular = eqn_params->mu1;
                break;
            case GAS_COMP2:
                mu_molecular = eqn_params->mu2;
                break;
            default:
                printf("\nERROR computeEddyViscosity2d(): unrecognized component!\n");
                printf("\t\tcomp = %d\n", comp);
                LOC(); clean_up(EXIT_FAILURE);
        }

        //TODO: Use model specified by eqn_params->eddy_viscosity_model
        mu[index] = mu_molecular + computeEddyViscosityVremanModel_BdryAware(icoords);
            //mu[index] = mu_molecular + computeEddyViscosityVremanModel(icoords);

        //TODO: Store computed eddy viscosity and use to compute the turbulent kinetic energy
        //      (Vreman Model) k_turb = 2.0*mu_turb*|S|
        //      The computed k_turb is then used in other SGS terms.
        
        if (mu[index] > mu_max)
        {
            mu_max = mu[index];
        }
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
        
        COMPONENT comp = top_comp[index];
        if (!gas_comp(comp))
        {
            mu[index] = 0.0;
            continue;
        }

        double mu_molecular;
        switch (comp)
        {
            case GAS_COMP1:
                mu_molecular = eqn_params->mu1;
                break;
            case GAS_COMP2:
                mu_molecular = eqn_params->mu2;
                break;
            default:
                printf("\nERROR computeEddyViscosity3d(): unrecognized component!\n");
                printf("\t\tcomp = %d\n", comp);
                LOC(); clean_up(EXIT_FAILURE);
        }

        //TODO: Use model specified by eqn_params->eddy_viscosity_model
        mu[index] = mu_molecular + computeEddyViscosityVremanModel_BdryAware(icoords);
            //mu[index] = mu_molecular + computeEddyViscosityVremanModel(icoords);
    
        //TODO: Store computed eddy viscosity and use to compute the turbulent kinetic energy
        //      (Vreman Model) k_turb = 2.0*mu_turb*|S|
        //      The computed k_turb is then used in other SGS terms.
        
        if (mu[index] > mu_max)
        {
            mu_max = mu[index];
        }
    }
}

double G_CARTESIAN::computeEddyViscosityVremanModel(int* icoords)
{
    double **vel = field.vel;

    double C_v = eqn_params->C_v;
        //double C_v = 0.025;
    
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
    double **vel = field.vel;

    double beta[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};

    double C_v = eqn_params->C_v;
        //double C_v = 0.025;
    

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
            else if (wave_type(hs) == ELASTIC_BOUNDARY)
            {
                //NOTE: zeroing/reducing the relative fluid velocity normal to
                //      canopy via the slip boundary condition causes the pressure
                //      jump imposed by the ghost fluid method porosity model
                //      to vanish -- destroying the porosity of the canopy
                //
                //      Above referring to iFluid solver -- what can/can't we do here?
                int index_nb = next_index_in_dir(icoords,dir[m][nb],top_gmax,dim);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                     wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                double v_slip[MAXD] = {0.0};
                setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,vel,v_slip);
                vel_nb[nb] = v_slip[l];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                if (boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),"cF_flowThroughBoundaryState") == 0)
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

//TODO: Make function return the computed slip velocity,
//      v_slip, as a std::vector<double>
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
    //TODO: dir is unused here?
    GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

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
    double range = 2;

    //TODO: Why does this fail for INCLUDE_BOUNDARIES and NO_SUBDOMAIN values?
    //      Conversely, why does it work with NO_BOUNDARIES in the backward facing
    //      step scenario -- to what degree is it working?
    
    /*
    FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,NO_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */
    /*      
    FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,INCLUDE_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */
    FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,NO_SUBDOMAIN,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);

    //TODO: We should get the ring of tris around the nearest interface point,
    //      and possible consider other nearby interface points that are within
    //      range. For a complex interface such as the human vtk model, there appears
    //      to be some error in the fluid region between the head and the hands.
    //      Excessively large velocities -- maybe not enough drag from the other
    //      nearby interface points that aren't being taken into account.

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

    /////////////////////////////////////////////////////////////////////////
    if (eqn_params->use_eddy_viscosity == NO)
    {
        /*
        for (int j = 0; j < dim; ++j)
            v_slip[j] = vel_reflect[j] + vel_ghost_nor[j];
        */

        //What about this?
        double vel_ghost_rel[MAXD] = {0.0};
        for (int j = 0; j < dim; ++j)
        {
            vel_ghost_rel[j] = vel_rel_tan[j] + vel_ghost_nor[j];
            v_slip[j] = vel_ghost_rel[j] + vel_intfc[j];
        }

        return;
    }
    /////////////////////////////////////////////////////////////////////////

    if (debugging("slip_boundary"))
    {
        fprint_general_vector(stdout,"vel_reflect",vel_reflect,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_tan",vel_rel_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_nor",vel_rel_nor,dim,"\n");
        printf("Magd(vel_rel_tan,dim) = %g\n",mag_vtan);
    }

    double mu_l;
    double rho_l;

    switch (comp)
    {
        case GAS_COMP1:
            mu_l = m_mu[0];
            rho_l = m_dens[0];
            break;
        case GAS_COMP2:
            mu_l = m_mu[1];
            rho_l = m_dens[1];
            break;
        default:
            printf("Unknown fluid COMPONENT: %d\n",comp);
            LOC(); clean_up(EXIT_FAILURE);
            break;
    }
    
    double tau_wall[MAXD] = {0.0};
    double mag_tau_wall = computeWallShearStress(mag_vtan,
                    dist_reflect,mu_l,rho_l,45.0);
    //NOTE: In all numerical experiments, Newton's method converged
    //      when the initial guess for the dimensionless wall velocity
    //      was in the range of 40-50.

    if (mag_vtan > MACH_EPS)
    {
        for (int j = 0; j < dim; ++j)
            tau_wall[j] = mag_tau_wall*vel_rel_tan[j]/mag_vtan;
    }

    // Interpolate the effective viscosity at the reflected point
    double mu_reflect;
    /*
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu,
                getStateMu,&mu_reflect,nullptr);
    if (mu_reflect < MACH_EPS) mu_reflect = field.mu[index]; //TODO: Need this?
    */
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu,
                getStateMu,&mu_reflect,&field.mu[index]);
    
    double vel_ghost_tan[MAXD] = {0.0};
    double vel_ghost_rel[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] = vel_rel_tan[j]
            - (dist_reflect - dist_ghost)/mu_reflect*tau_wall[j];

        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        v_slip[j] = vel_ghost_rel[j] + vel_intfc[j];
    }


    if (debugging("slip_boundary"))
    {
        printf("mu_reflect = %g , mu_[%d] = %g\n",mu_reflect,index,field.mu[index]);
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

