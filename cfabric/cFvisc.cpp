#include "cFluid.h"
#include "cFturb.h"

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};


void G_CARTESIAN::addViscousFlux(
        SWEEP* m_vst,
        FSWEEP *m_flux,
        double delta_t)
{
    switch (dim)
    {
    case 2:
        {
            for (int j = imin[1]; j <= imax[1]; j++)
            for (int i = imin[0]; i <= imax[0]; i++)
            {   
                int icoords[MAXD] = {i,j,0};
                int index = d_index(icoords,top_gmax,dim);
                if (!gas_comp(top_comp[index])) continue;
                
                VStencil2d vsten;
                fillViscousFluxStencil2d(icoords,m_vst,&vsten);
                
                VFLUX v_flux;
                computeViscousFlux2d(icoords,m_vst,&v_flux,delta_t,&vsten);

                for (int k = 0; k < dim; ++k)
                    m_flux->momn_flux[k][index] += v_flux.momn_flux[k];
                m_flux->engy_flux[index] += v_flux.engy_flux;
            }
            break;
        }
    case 3:
        {
            for (int k = imin[2]; k <= imax[2]; k++)
            for (int j = imin[1]; j <= imax[1]; j++)
            for (int i = imin[0]; i <= imax[0]; i++)
            {   
                int icoords[MAXD] = {i,j,k};
                int index = d_index(icoords,top_gmax,dim);
                if (!gas_comp(top_comp[index])) continue;
                
                VStencil3d vsten;
                fillViscousFluxStencil3d(icoords,m_vst,&vsten);

                VFLUX v_flux;
                computeViscousFlux3d(icoords,m_vst,&v_flux,delta_t,&vsten);

                for (int k = 0; k < dim; ++k)
                    m_flux->momn_flux[k][index] += v_flux.momn_flux[k];
                m_flux->engy_flux[index] += v_flux.engy_flux;
            }
            break;
        }
    default:
        {
            printf("addViscousFlux() ERROR: Invalid Dimension\n");
            LOC(); clean_up(EXIT_FAILURE);
        }
    }
}

void G_CARTESIAN::fillViscousFluxStencil2d(
        int* icoords,
        SWEEP* m_vst,
        VStencil2d* vsten)
{
    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];

    for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
    {
        VSWEEP* vs = &vsten->st[j][i];
        vs->icoords[0] = icoords[0] + i - 1;
        vs->icoords[1] = icoords[1] + j - 1;
        
        int idx_nb = d_index(vs->icoords,top_gmax,dim);
        vs->comp = top_comp[idx_nb];
        
        if (vs->comp != comp)
        {
            setViscousGhostState(icoords,comp,vs,m_vst);
        }
        else
        {
            for (int l = 0; l < dim; ++l)
            {
                vs->vel[l] = m_vst->momn[l][idx_nb]/m_vst->dens[idx_nb];
            }
        }
    }
}

void G_CARTESIAN::fillViscousFluxStencil3d(
        int* icoords,
        SWEEP* m_vst,
        VStencil3d* vsten)
{
    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];

    for (int k = 0; k < 3; ++k)
    for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
    {
        VSWEEP* vs = &vsten->st[k][j][i];
        vs->icoords[0] = icoords[0] + i - 1;
        vs->icoords[1] = icoords[1] + j - 1;
        vs->icoords[2] = icoords[2] + k - 1;
        
        int idx_nb = d_index(vs->icoords,top_gmax,dim);
        vs->comp = top_comp[idx_nb];
        
        if (vs->comp != comp)
        {
            setViscousGhostState(icoords,comp,vs,m_vst);
        }
        else
        {
            for (int l = 0; l < dim; ++l)
            {
                vs->vel[l] = m_vst->momn[l][idx_nb]/m_vst->dens[idx_nb];
            }
        }
    }
}

void G_CARTESIAN::setViscousGhostState(
        int* icoords,
        COMPONENT comp,
        VSWEEP* vs,
        SWEEP* m_vst)
{
    double nip_coords[MAXD];
    double intrp_coeffs[MAXD];
    HYPER_SURF* hs;
    HYPER_SURF_ELEMENT* hse;
    
    int ghost_index = d_index(vs->icoords,top_gmax,dim);
    auto ghost_coords = cell_center[ghost_index].getCoords();
    
    bool nip_found = nearest_interface_point(&ghost_coords[0],
                comp,front->interf,INCLUDE_BOUNDARIES,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);
    /*bool nip_found = nearest_interface_point(&ghost_coords[0],
                comp,front->interf,NO_SUBDOMAIN,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);*/
    
    if (!nip_found)
    {
        printf("setViscousGhostState() ERROR: "
                "can't find nearest interface point\n");
        LOC(); clean_up(EXIT_FAILURE);
    }

    switch (wave_type(hs))
    {
        case DIRICHLET_BOUNDARY:
        {
            setDirichletViscousGhostState(vs,comp,intrp_coeffs,hse,hs);
            break;
        }
        case NEUMANN_BOUNDARY:
        case MOVABLE_BODY_BOUNDARY:
        {
            setNeumannViscousGhostState(icoords,m_vst,vs,&ghost_coords[0],
                    nip_coords,comp,intrp_coeffs,hse,hs);
            break;
        }
        /*case ELASTIC_BOUNDARY:
        {
            //setElasticViscousGhostState(vs,comp,intrp_coeffs,hse,hs,state);
            break;
        }*/
        default:
        {
            printf("\n\nsetViscousGhostState() ERROR: unknown boundary type\n\n");
            LOC(); clean_up(EXIT_FAILURE);
        }
    }

}

void G_CARTESIAN::setDirichletViscousGhostState(
        VSWEEP* vs,
        COMPONENT comp,
        double* intrp_coeffs,
        HYPER_SURF_ELEMENT* hse,
        HYPER_SURF* hs)
{
    STATE* state;
    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));

    if (boundary_state(hs) != nullptr)
        ft_assign((POINTER)state,boundary_state(hs),front->sizest);
    else
        state_along_hypersurface_element(comp,intrp_coeffs,hse,hs,(POINTER)state);

    for (int i = 0; i < dim; ++i)
        vs->vel[i] = state->vel[i];

    FT_FreeThese(1,state);
}

void G_CARTESIAN::setNeumannViscousGhostState(
        int* icoords,
        SWEEP* m_vst,
        VSWEEP* vs,
        double* ghost_coords,
        double* crx_coords,
        COMPONENT comp,
        double* intrp_coeffs,
        HYPER_SURF_ELEMENT* hse,
        HYPER_SURF* hs)
{
    double nor[MAXD] = {0.0};
    double vel_intfc[MAXD] = {0.0};

    switch (dim)
    {
    case 2:
        {
            STATE* ss;
            STATE* se;

            if (gas_comp(negative_component(hs)))
            {
                ss = (STATE*)left_state(Bond_of_hse(hse)->start);
                se = (STATE*)left_state(Bond_of_hse(hse)->end);
            }
            else if (gas_comp(positive_component(hs)))
            {
                ss = (STATE*)right_state(Bond_of_hse(hse)->start);
                se = (STATE*)right_state(Bond_of_hse(hse)->end);
            }
            else
            {
                printf("setNeumannViscousGhostState() ERROR: "
                        "no gas component on hypersurface\n");
                LOC(); clean_up(EXIT_FAILURE);
            }

            double ns[MAXD] = {0.0};
            double ne[MAXD] = {0.0};
            
            normal(Bond_of_hse(hse)->start,hse,hs,ns,front);
            normal(Bond_of_hse(hse)->end,hse,hs,ne,front);

            for (int i = 0; i < dim; ++i)
            {
                nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
            }
        }
        break;

    case 3:
        {
            STATE* st[3];
            TRI* nearTri = Tri_of_hse(hse);

            if (gas_comp(negative_component(hs)))
            {
                for (int j = 0; j < 3; ++j)
                    st[j] = (STATE*)left_state(Point_of_tri(nearTri)[j]);
            }
            else if (gas_comp(positive_component(hs)))
            {
                for (int j = 0; j < 3; ++j)
                    st[j] = (STATE*)right_state(Point_of_tri(nearTri)[j]);
            }
            else
            {
                printf("setNeumannViscousGhostState() ERROR: "
                        "no gas component on hypersurface\n");
                LOC(); clean_up(EXIT_FAILURE);
            }

            //NOTE: Tri_normal() does not return a unit vector
            const double* tnor = Tri_normal(nearTri);

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

    //NOTE: must use unit-length vectors with FT_GridSizeInDir()
    double mag_nor = Magd(nor,dim);
    for (int i = 0; i < dim; ++i)
        nor[i] /= mag_nor;
        
    if (comp == negative_component(hs))
	{
	    for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
	}
    
    double dist_reflect = FT_GridSizeInDir(nor,front);
    double dist_ghost = distance_between_positions(ghost_coords,crx_coords,dim);
    
    //The desired reflected point
    double coords_reflect[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
    {
        coords_reflect[j] = crx_coords[j] + dist_reflect*nor[j];
    }

    int index = d_index(icoords,top_gmax,dim);
    //TODO: ADD DEBUG BLOCK HERE -- LIKE setSlipBoundaryNIP()

    //Interpolate Density and Momentum at the reflected point
    //and compute the velocity.
    double dens_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,nullptr);
    /*FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,&m_vst->dens[index]);*/

    double mom_reflect[MAXD];
    double vel_reflect[MAXD];
    for (int j = 0; j < dim; ++j)
    {
        //TODO: Make sure the interface value has been set in case it is
        //      needed as the default value (since we are passing nullptr
        //      as the defailt value). Same for above.
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&mom_reflect[j],nullptr);
        /*FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&momn_reflect[j],&m_vst->momn[j][index]);*/
        vel_reflect[j] = mom_reflect[j]/dens_reflect;
    }


    /*
    /////////////////////////////////////////////////////////////////////////////////////
    for (int j = 0; j < dim; ++j)
    {
        vs->vel[j] = vel_reflect[j];
    }
    /////////////////////////////////////////////////////////////////////////////////////
    */

    
    
    /////////////////////////////////////////////////////////////////////////////////////
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
    EQN_PARAMS *eqn_params = (EQN_PARAMS*)front->extra1;
    if (eqn_params->use_eddy_viscosity == NO)
    {
        for (int j = 0; j < dim; ++j)
        {
            vs->vel[j] = vel_reflect[j] + vel_ghost_nor[j];
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
    
    //NOTE: In all numerical experiments, Newton's method converged
    //      when the initial guess for the dimensionless wall velocity
    //      was in the range of 40-50.
    double tau_wall[MAXD] = {0.0};
    double mag_tau_wall = computeWallShearStress(mag_vtan,
                    dist_reflect,mu_l,rho_l,45.0);
    if (mag_vtan > MACH_EPS)
    {
        for (int j = 0; j < dim; ++j)
            tau_wall[j] = mag_tau_wall*vel_rel_tan[j]/mag_vtan;
    }

    // Interpolate the effective viscosity at the reflected point
    double mu_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu,
                getStateMu,&mu_reflect,nullptr);
        //FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu,
        //      getStateMu,&mu_reflect,&field.mu[index]);
    if (mu_reflect < MACH_EPS) mu_reflect = field.mu[index]; //TODO: Need this?
    
    double vel_ghost_tan[MAXD] = {0.0};
    double vel_ghost_rel[MAXD] = {0.0};
    double v_slip[MAXD] = {0.0};
    
    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] = vel_rel_tan[j]
            - (dist_reflect - dist_ghost)/mu_reflect*tau_wall[j];

        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        
        v_slip[j] = vel_ghost_rel[j] + vel_intfc[j];
        
        vs->vel[j] = v_slip[j];
        //vs->vel[j] = vel_ghost_rel[j] + vel_intfc[j];
    }
    /////////////////////////////////////////////////////////////////////////////////////
    
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
}

void G_CARTESIAN::computeViscousFlux2d(
        int* icoords,
        SWEEP* m_vst,
        VFLUX* v_flux,
        double delta_t,
        VStencil2d* vsten)
{
    auto sten = vsten->st;

    double u = sten[1][1].vel[0];
    double v = sten[1][1].vel[1];

    double u_x = 0.5*(sten[1][2].vel[0] - sten[1][0].vel[0])/top_h[0];
    double u_y = 0.5*(sten[2][1].vel[0] - sten[0][1].vel[0])/top_h[1];
    double v_x = 0.5*(sten[1][2].vel[1] - sten[1][0].vel[1])/top_h[0];
    double v_y = 0.5*(sten[2][1].vel[1] - sten[0][1].vel[1])/top_h[1];

    double u_xx = (sten[1][2].vel[0] - 2.0*sten[1][1].vel[0]
            + sten[1][0].vel[0])/sqr(top_h[0]);
    double u_yy = (sten[2][1].vel[0] - 2.0*sten[1][1].vel[0] 
            + sten[0][1].vel[0])/sqr(top_h[1]);
    double v_xx = (sten[1][2].vel[1] - 2.0*sten[1][1].vel[1]
            + sten[1][0].vel[1])/sqr(top_h[0]);
    double v_yy = (sten[2][1].vel[1] - 2.0*sten[1][1].vel[1]
            + sten[0][1].vel[1])/sqr(top_h[1]);
    
    double u_xy = 0.25*(sten[2][2].vel[0] - sten[2][0].vel[0]
            - sten[0][2].vel[0] + sten[0][0].vel[0])/top_h[0]/top_h[1];
    double v_xy = 0.25*(sten[2][2].vel[1] - sten[2][0].vel[1]
            - sten[0][2].vel[1] + sten[0][0].vel[1])/top_h[0]/top_h[1];
    
    //TODO: Get mu from m_vst instead?
    double* mu = field.mu;
    int index = d_index(icoords,top_gmax,dim);
    
    double tauxx = 2.0/3.0*mu[index]*(2.0*u_x - v_y);
    double tauyy = 2.0/3.0*mu[index]*(2.0*v_y - u_x);
    double tauxy = mu[index]*(u_y + v_x);

    double tauxx_x = 2.0/3.0*mu[index]*(2.0*u_xx - v_xy);
    double tauyy_y = 2.0/3.0*mu[index]*(2.0*v_yy - u_xy);
    double tauxy_y = mu[index]*(u_yy + v_xy);
    double tauxy_x = mu[index]*(u_xy + v_xx);
    
    v_flux->momn_flux[0] = delta_t*(tauxx_x + tauxy_y);
    v_flux->momn_flux[1] = delta_t*(tauxy_x + tauyy_y);
    
    v_flux->engy_flux = delta_t*(u_x*tauxx + u*tauxx_x + v_x*tauxy
            + v*tauxy_y + u_y*tauxy + u*tauxy_y + v_y*tauyy + v*tauyy_y);
}

void G_CARTESIAN::computeViscousFlux3d(
        int* icoords,
        SWEEP* m_vst,
        VFLUX* v_flux,
        double delta_t,
        VStencil3d* vsten)
{
    auto sten = vsten->st;

    double u = sten[1][1][1].vel[0];
    double v = sten[1][1][1].vel[1];
    double w = sten[1][1][1].vel[2];

    double u_x = 0.5*(sten[1][1][2].vel[0] - sten[1][1][0].vel[0])/top_h[0];
    double u_y = 0.5*(sten[1][2][1].vel[0] - sten[1][0][1].vel[0])/top_h[1];
    double u_z = 0.5*(sten[2][1][1].vel[0] - sten[0][1][1].vel[0])/top_h[2];

    double v_x = 0.5*(sten[1][1][2].vel[1] - sten[1][1][0].vel[1])/top_h[0];
    double v_y = 0.5*(sten[1][2][1].vel[1] - sten[1][0][1].vel[1])/top_h[1];
    double v_z = 0.5*(sten[2][1][1].vel[1] - sten[0][1][1].vel[1])/top_h[2];
    
    double w_x = 0.5*(sten[1][1][2].vel[2] - sten[1][1][0].vel[2])/top_h[0];
    double w_y = 0.5*(sten[1][2][1].vel[2] - sten[1][0][1].vel[2])/top_h[1];
    double w_z = 0.5*(sten[2][1][1].vel[2] - sten[0][1][1].vel[2])/top_h[2];

    double u_xx = (sten[1][1][2].vel[0] - 2.0*sten[1][1][1].vel[0]
            + sten[1][1][0].vel[0])/sqr(top_h[0]);
    double u_yy = (sten[1][2][1].vel[0] - 2.0*sten[1][1][1].vel[0] 
            + sten[1][0][1].vel[0])/sqr(top_h[1]);
    double u_zz = (sten[2][1][1].vel[0] - 2.0*sten[1][1][1].vel[0]
            + sten[0][1][1].vel[0])/sqr(top_h[2]);

    double v_xx = (sten[1][1][2].vel[1] - 2.0*sten[1][1][1].vel[1]
            + sten[1][1][0].vel[1])/sqr(top_h[0]);
    double v_yy = (sten[1][2][1].vel[1] - 2.0*sten[1][1][1].vel[1]
            + sten[1][0][1].vel[1])/sqr(top_h[1]);
    double v_zz = (sten[2][1][1].vel[1] - 2.0*sten[1][1][1].vel[1]
            + sten[0][1][1].vel[1])/sqr(top_h[2]);
    
    double w_xx = (sten[1][1][2].vel[2] - 2.0*sten[1][1][1].vel[2]
            + sten[1][1][0].vel[2])/sqr(top_h[0]);
    double w_yy = (sten[1][2][1].vel[2] - 2.0*sten[1][1][1].vel[2]
            + sten[1][0][1].vel[2])/sqr(top_h[1]);
    double w_zz = (sten[2][1][1].vel[2] - 2.0*sten[1][1][1].vel[2]
            + sten[0][1][1].vel[2])/sqr(top_h[2]);

    double u_xy = 0.25*(sten[1][2][2].vel[0] - sten[1][2][0].vel[0]
            - sten[1][0][2].vel[0] + sten[1][0][0].vel[0])/top_h[0]/top_h[1];
    double v_xy = 0.25*(sten[1][2][2].vel[1] - sten[1][2][0].vel[1]
            - sten[1][0][2].vel[1] + sten[1][0][0].vel[1])/top_h[0]/top_h[1];
    double w_xy = 0.25*(sten[1][2][2].vel[2] - sten[1][2][0].vel[2]
            - sten[1][0][2].vel[2] + sten[1][0][0].vel[2])/top_h[0]/top_h[1];

    double u_xz = 0.25*(sten[2][1][2].vel[0] - sten[2][1][0].vel[0]
            - sten[0][1][2].vel[0] + sten[0][1][0].vel[0])/top_h[0]/top_h[2];
    double v_xz = 0.25*(sten[2][1][2].vel[1] - sten[2][1][0].vel[1]
            - sten[0][1][2].vel[1] + sten[0][1][0].vel[1])/top_h[0]/top_h[2];
    double w_xz = 0.25*(sten[2][1][2].vel[2] - sten[2][1][0].vel[2]
            - sten[0][1][2].vel[2] + sten[0][1][0].vel[2])/top_h[0]/top_h[2];

    double u_yz = 0.25*(sten[2][2][1].vel[0] - sten[2][0][1].vel[0]
            - sten[0][2][1].vel[0] + sten[0][0][1].vel[0])/top_h[1]/top_h[2];
    double v_yz = 0.25*(sten[2][2][1].vel[1] - sten[2][0][1].vel[1]
            - sten[0][2][1].vel[1] + sten[0][0][1].vel[1])/top_h[1]/top_h[2];
    double w_yz = 0.25*(sten[2][2][1].vel[2] - sten[2][0][1].vel[2]
            - sten[0][2][1].vel[2] + sten[0][0][1].vel[2])/top_h[1]/top_h[2];

    //TODO: Get mu from m_vst instead?
    double* mu = field.mu;
    int index = d_index(icoords,top_gmax,dim);
    
    double tauxx = 2.0/3.0*mu[index]*(2.0*u_x - v_y - w_z);
    double tauyy = 2.0/3.0*mu[index]*(2.0*v_y - u_x - w_z);
    double tauzz = 2.0/3.0*mu[index]*(2.0*w_z - u_x - v_y);
    
    double tauxy = mu[index]*(u_y + v_x);
    double tauxz = mu[index]*(u_z + w_x);
    double tauyz = mu[index]*(v_z + w_y);

    double tauxx_x = 2.0/3.0*mu[index]*(2.0*u_xx - v_xy - w_xz);
    double tauyy_y = 2.0/3.0*mu[index]*(2.0*v_yy - u_xy - w_yz);
    double tauzz_z = 2.0/3.0*mu[index]*(2.0*w_zz - u_xz - v_yz);

    double tauxy_x = mu[index]*(u_xy + v_xx);
    double tauxz_x = mu[index]*(u_xz + w_xx);
    double tauxy_y = mu[index]*(u_yy + v_xy);
    double tauyz_y = mu[index]*(v_yz + w_yy);
    double tauxz_z = mu[index]*(u_zz + w_xz);
    double tauyz_z = mu[index]*(v_zz + w_yz);
    
    v_flux->momn_flux[0] = delta_t*(tauxx_x + tauxy_y + tauxz_z);
    v_flux->momn_flux[1] = delta_t*(tauxy_x + tauyy_y + tauyz_z);
    v_flux->momn_flux[2] = delta_t*(tauxz_x + tauyz_y + tauzz_z);
    
    v_flux->engy_flux = delta_t*(u_x*tauxx + u*tauxx_x + v_x*tauxy
            + v*tauxy_y + w_x*tauxz + w*tauxz_x + u_y*tauxy + u*tauxy_y
            + v_y*tauyy + v*tauyy_y + w_y*tauyz + w*tauyz_y + u_z*tauxz
            + u*tauxz_z + v_z*tauyz + v*tauyz_z + w_z*tauzz + w*tauzz_z);
}


