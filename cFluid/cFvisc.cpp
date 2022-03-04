#include "cFluid.h"
#include "cFturb.h"

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};


//TODO: Should rename to reflect that temperature and dynamic viscosity
//      are being computed/updated.
void G_CARTESIAN::computeDynamicViscosity()
{
    mu_max = 0.0;
    switch (dim)
    {
        case 2:
            computeDynamicViscosity2d();
            break;
        case 3:
            computeDynamicViscosity3d();
            break;
        default:
            printf("\nERROR computeDynamicViscosity(): invalid dim\n");
            printf("\t\t dim = %d\n\n",dim);
            LOC(); clean_up(EXIT_FAILURE);
    }
}

void G_CARTESIAN::computeDynamicViscosity2d()
{
    double* dens = field.dens;
    double* pres = field.pres;
    double* temp = field.temp;
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
            temp[index] = 0.0;
            continue;
        }

        STATE state;
        state.eos = &eqn_params->eos[comp];
        state.dens = dens[index];
        state.pres = pres[index];
        
        state.temp = EosTemperature(&state);
        mu[index] = EosViscosity(&state);
        temp[index] = state.temp;

        if (mu[index] > mu_max)
        {
            mu_max = mu[index];
        }
    }
}

void G_CARTESIAN::computeDynamicViscosity3d()
{
    double* dens = field.dens;
    double* pres = field.pres;
    double* temp = field.temp;
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
            temp[index] = 0.0;
            continue;
        }

        STATE state;
        state.eos = &eqn_params->eos[comp];
        state.dens = dens[index];
        state.pres = pres[index];

        state.temp = EosTemperature(&state);
        mu[index] = EosViscosity(&state);
        temp[index] = state.temp;

        if (mu[index] > mu_max)
        {
            mu_max = mu[index];
        }
    }
}

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
                
                COMPONENT comp = top_comp[index];
                if (!gas_comp(comp)) continue;
                
                VStencil2d vsten;
                vsten.comp = comp;
                fillViscousFluxStencil2d(icoords,m_vst,&vsten);
                
                VFLUX v_flux;
                computeViscousFlux2d(icoords,m_vst,&v_flux,delta_t,&vsten);

                /*
                //TODO: Need to deal with stencil at domain boundary for this to work.
                //      Can we just use the i+-1 state for the i+-2 state beyond the domain?
                //      Or can use a 3pt stencil near the domain boundary ???

                VStencil2d_5pt vsten;
                fillViscousFluxStencil2d_5pt(icoords,m_vst,&vsten);
                
                VFLUX v_flux;
                computeViscousFlux2d_5pt(icoords,m_vst,&v_flux,delta_t,&vsten);
                */

                for (int k = 0; k < dim; ++k)
                {
                    m_flux->momn_flux[k][index] += v_flux.momn_flux[k];
                }
                
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
                
                COMPONENT comp = top_comp[index];
                if (!gas_comp(comp)) continue;
                
                VStencil3d vsten;
                vsten.comp = comp;
                fillViscousFluxStencil3d(icoords,m_vst,&vsten);

                VFLUX v_flux;
                computeViscousFlux3d(icoords,m_vst,&v_flux,delta_t,&vsten);

                for (int k = 0; k < dim; ++k)
                {
                    m_flux->momn_flux[k][index] += v_flux.momn_flux[k];
                }

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
    COMPONENT comp = vsten->comp;

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
            
            vs->temp = m_vst->temp[idx_nb];
        }
    }
}

//For 4th order approximation to viscous flux 
//add points i-2 and i+2 to stencil in each direction.
void G_CARTESIAN::fillViscousFluxStencil2d_5pt(
        int* icoords,
        SWEEP* m_vst,
        VStencil2d_5pt* vsten)
{
    COMPONENT comp = vsten->comp;

    //TODO: buffer not wide enough near boundary
    //
    //      need to copy the i+1 state to the i+2 state
    //      if i+1 point is on rect domain boundary
    //
    //      Or can we just revert to a 3pt stencil near boundary?
    //      Can add a flag indicating that the 3 pt stencil is in use.
    
    for (int j = 0; j < 5; ++j)
    for (int i = 0; i < 5; ++i)
    {
        VSWEEP* vs = &vsten->st[j][i];
        vs->icoords[0] = icoords[0] + i - 2;
        vs->icoords[1] = icoords[1] + j - 2;
        
        int idx_nb = d_index(vs->icoords,top_gmax,dim);
        vs->comp = top_comp[idx_nb];
        
        //TODO: This will only see rigid body wall boundaries -- no porous walls etc.
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

            vs->temp = m_vst->temp[idx_nb];
        }
    }
}

void G_CARTESIAN::fillViscousFluxStencil3d(
        int* icoords,
        SWEEP* m_vst,
        VStencil3d* vsten)
{
    /*
    GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    INTERFACE* grid_intfc = front->grid_intfc;
    HYPER_SURF* hs;
    STATE* state;
    double crx_coords[MAXD];
    */

    COMPONENT comp = vsten->comp;

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

        //TODO: Checking all the necessary crossings could be very expensive.
        //      Better if we can find a way to have the index_coating algorithm
        //      components persist for the fluid solver to use.
        //
        //TODO: Can try using this to check:
        //
        //      Table* T = front->interf->table;
        //      if (T->compon3d[iz][iy][ix] == ONFRONT)
        
        /*
        int status = FT_StateStructAtGridCrossing(front,grid_intfc,
                icoords,dir[idir][nb],comp,(POINTER*)&state,&hs,crx_coords);
        */

        //TODO: This will only see rigid body wall boundaries -- no porous walls etc.
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

            vs->temp = m_vst->temp[idx_nb];
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
    COMPONENT ghost_comp = vs->comp;

    //TODO: Why can't we use nearest_intfc_point_in_range()???
    int range = 2;
    //int range = 3;
    

    /*
    bool nip_found =
    FT_FindNearestIntfcPointInRange(front,ghost_comp,&ghost_coords[0],
            INCLUDE_BOUNDARIES,nip_coords,intrp_coeffs,&hse,&hs,range);
    */


    /*
    bool nip_found = nearest_interface_point_within_range(&ghost_coords[0],
                comp,front->interf,NO_BOUNDARIES,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs,3);
    */

    /*
    bool nip_found = nearest_interface_point_within_range(&ghost_coords[0],
                comp,front->interf,NO_SUBDOMAIN,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs,3);
    */
    
    /*
    bool nip_found = nearest_interface_point(&ghost_coords[0],
                comp,front->interf,NO_BOUNDARIES,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);
    */
    
    /*
    bool nip_found = nearest_interface_point(&ghost_coords[0],
                comp,front->interf,INCLUDE_BOUNDARIES,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);
    */
    
    bool nip_found = nearest_interface_point(&ghost_coords[0],
                comp,front->interf,NO_SUBDOMAIN,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);


    if (!nip_found)
    {
        printf("ERROR G_CARTESIAN::setViscousGhostState(): "
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
        case MOVABLE_BODY_BOUNDARY:
        case NEUMANN_BOUNDARY:
        {
            setNeumannViscousGhostState(icoords,m_vst,vs,&ghost_coords[0],
                    nip_coords,comp,intrp_coeffs,hse,hs);
            break;
        }
        case SYMMETRY_BOUNDARY:
        {
            setSymmetryViscousGhostState(icoords,m_vst,vs,&ghost_coords[0],
                    nip_coords,comp,intrp_coeffs,hse,hs);
            break;
        }
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

    //TODO: Need to deal with the case when these conditions are both true.
    //      e.g. Nodes/Curves of domain boundary interface where where a
    //       constant inlet boundary meets a flowthru/outlet boundary.
    if (boundary_state_function(hs) != nullptr)
    {
        state_along_hypersurface_element(comp,intrp_coeffs,hse,hs,(POINTER)state);
    }
    else if (boundary_state(hs) != nullptr)
    {
        ft_assign((POINTER)state,boundary_state(hs),front->sizest);
    }
    else
    {
        printf("\n\nERROR setDirichletViscousGhostState(): unrecognized DIRICHLET_BOUNDARY\n");
        LOC(); clean_up(EXIT_FAILURE);
    }
    
    for (int i = 0; i < dim; ++i)
        vs->vel[i] = state->vel[i];

    vs->temp = state->temp;

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
    double dens_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,&m_vst->dens[index]);
    /*
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,nullptr);
    */

    double mom_reflect[MAXD];
    double vel_reflect[MAXD];
    for (int j = 0; j < dim; ++j)
    {
        //TODO: Make sure the interface value has been set in case it is
        //      needed as the default value (since we are passing nullptr
        //      as the defailt value). Same for above.
        /*
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&mom_reflect[j],nullptr);
        */
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&mom_reflect[j],&m_vst->momn[j][index]);
        
        vel_reflect[j] = mom_reflect[j]/dens_reflect;
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

    EOS_PARAMS eos = eqn_params->eos[comp];
    double R_specific = eos.R_specific;
    double gamma = eos.gamma;
    double Pr = eos.Pr;
    
    double Cp = gamma/(gamma - 1.0)*R_specific;
    double sqrmag_vel_tan = Dotd(vel_rel_tan, vel_rel_tan, dim);

    double temp_reflect;
    double temp_wall = eqn_params->fixed_wall_temp;
    if (!eqn_params->use_fixed_wall_temp)
    {
        //Interpolate the temperature at the reflected point
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->temp,
                getStateTemp,&temp_reflect,&m_vst->temp[index]);

        //Compute Wall Temperature 
        temp_wall = temp_reflect + 0.5*pow(Pr,1.0/3.0)*sqrmag_vel_tan/Cp;
    }

    //Interpolate the pressure at the reflected point
    double pres_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->pres,
            getStatePres,&pres_reflect,&m_vst->pres[index]);

    //Compute density near wall using the wall temperature and the pressure at the reflected point
    double dens_wall = pres_reflect/temp_wall/R_specific;

    //Interpolate the viscosity at the reflected point
    double mu_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->mu,
            getStateMu,&mu_reflect,&m_vst->mu[index]);
    /*FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->mu,
            getStateMu,&mu_reflect,nullptr);*/
    
    if (std::isnan(mu_reflect) || std::isinf(mu_reflect))
    {
       printf("\nsetSlipBoundaryNIP() ERROR: nan/inf mu_reflect\n");
       printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,m_vst->mu[index]);
       LOC(); clean_up(EXIT_FAILURE);
    }
    
    /*
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu,
                    getStateMu,&mu_reflect,nullptr);
        if (mu_reflect < MACH_EPS) mu_reflect = field.mu[index]; //TODO: Need this?
    */
    
    double mu_turb_reflect = 0.0;
    if (eqn_params->use_eddy_viscosity)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->mu_turb,
                getStateMuTurb,&mu_turb_reflect,&m_vst->mu_turb[index]);
        mu_reflect += mu_turb_reflect;
    }
    
    double tau_wall[MAXD] = {0.0};
    
    //double mag_tau_wall = computeWallShearStress(mag_vtan,dist_reflect,mu_l,rho_l,45.0);
    double mag_tau_wall = computeWallShearStress(mag_vtan,dist_reflect,mu_reflect,dens_wall,100.0);
        //double mag_tau_wall = computeWallShearStress(mag_vtan,dist_reflect,mu_reflect,dens_reflect,45.0);
    
    if (mag_vtan > MACH_EPS)
    {
        for (int j = 0; j < dim; ++j)
            tau_wall[j] = mag_tau_wall*vel_rel_tan[j]/mag_vtan;
    }

    double vel_ghost_tan[MAXD] = {0.0};
    double vel_ghost_rel[MAXD] = {0.0};
    double v_slip[MAXD] = {0.0};
    
    //TODO: Should it be:
    //                       (dist_reflect + dist_ghost)/mu_reflect;
    //      instead?                 
    double coeff_tau = (mu_reflect == 0) ? 0.0 : (dist_reflect - dist_ghost)/mu_reflect;

    double slip = 1.0;
    if (eqn_params->no_slip_wall)
    {
        slip = 0.0;
    }

    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] = slip*vel_rel_tan[j] - coeff_tau*tau_wall[j];

        /*
        vel_ghost_tan[j] = vel_rel_tan[j]
            - (dist_reflect + dist_ghost)/mu_reflect*tau_wall[j];
        */

        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        
        v_slip[j] = vel_ghost_rel[j] + vel_intfc[j];
        
        vs->vel[j] = v_slip[j]; //vs->vel[j] = vel_ghost_rel[j] + vel_intfc[j];
    }
    
    
    double temp_ghost = eqn_params->fixed_wall_temp;
    if (!eqn_params->use_fixed_wall_temp)
    {
        double sqrmag_vel_ghost_tan = Dotd(vel_ghost_tan, vel_ghost_tan, dim);
        temp_ghost = temp_reflect + 0.5*pow(Pr,1.0/3.0)*(sqrmag_vel_tan - sqrmag_vel_ghost_tan)/Cp;
    }

    vs->temp = temp_ghost;

    //TODO: COMPUTE GHOST DENSITY WITH EQUATION OF STATE?
        //vs->dens = pres_reflect/temp_ghost/R_specific;

    
    if (std::isnan(v_slip[0]) || std::isinf(v_slip[0]) ||
        std::isnan(v_slip[1]) || std::isinf(v_slip[1]) ||
        std::isnan(v_slip[2]) || std::isinf(v_slip[2]))
    {
        printf("\nsetSlipBoundaryNIP() ERROR: nan/inf v_slip\n");
        fprint_int_vector(stdout,"icoords",icoords,dim,"\n");
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,m_vst->mu[index]);
        printf("mu_turb_reflect = %g , mu_turb[%d] = %g\n",mu_turb_reflect,index,m_vst->mu_turb[index]);
        LOC(); clean_up(EXIT_FAILURE);
    }

    if (debugging("slip_boundary"))
    {
        printf("mu_reflect = %g , mu[%d] = %g\n",mu_reflect,index,m_vst->mu[index]);
        printf("mag_tau_wall = %g\n",mag_tau_wall);
        fprint_general_vector(stdout,"tau_wall",tau_wall,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_rel",vel_ghost_rel,dim,"\n");
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
        printf("temp_ghost = %g\n",temp_ghost);
    }
}

void G_CARTESIAN::setSymmetryViscousGhostState(
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
    
    double dist_ghost = distance_between_positions(ghost_coords,crx_coords,dim);
    double dist_reflect = dist_ghost;
        //double dist_reflect = FT_GridSizeInDir(nor,front);
    
    //The desired reflected point
    double coords_reflect[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
    {
        coords_reflect[j] = crx_coords[j] + dist_reflect*nor[j];
    }

    int index = d_index(icoords,top_gmax,dim);

    //Interpolate the temperature at the reflected point
    double temp_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->temp,
            getStateTemp,&temp_reflect,&m_vst->temp[index]);

    vs->temp = temp_reflect;

    //Interpolate Density and Momentum at the reflected point
    double dens_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,&m_vst->dens[index]);

    double mom_reflect[MAXD];
    double vel_reflect[MAXD];
    for (int j = 0; j < dim; ++j)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&mom_reflect[j],&m_vst->momn[j][index]);
        
        vel_reflect[j] = mom_reflect[j]/dens_reflect;
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
    double vel_ghost_tan[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
	    vel_rel_tan[j] = vel_rel[j] - vn*nor[j];
	    vel_ghost_tan[j] = vel_rel_tan[j];
	    vel_rel_nor[j] = vn*nor[j];
	    vel_ghost_nor[j] = -1.0*(dist_ghost/dist_reflect)*vn*nor[j];
    }

    if (debugging("slip_boundary"))
    {
        fprint_general_vector(stdout,"vel_reflect",vel_reflect,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_tan",vel_rel_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_nor",vel_rel_nor,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
    }

    
    double vel_ghost_rel[MAXD] = {0.0};
    double v_sym[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        v_sym[j] = vel_ghost_rel[j] + vel_intfc[j];
        vs->vel[j] = v_sym[j];
    }
    
    if (std::isnan(v_sym[0]) || std::isinf(v_sym[0]) ||
        std::isnan(v_sym[1]) || std::isinf(v_sym[1]) ||
        std::isnan(v_sym[2]) || std::isinf(v_sym[2]))
    {
        printf("\nsetSymmetryViscousGhostState() ERROR: nan/inf v_sym\n");
        fprint_int_vector(stdout,"icoords",icoords,dim,"\n");
        fprint_general_vector(stdout,"v_sym",v_sym,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        LOC(); clean_up(EXIT_FAILURE);
    }

    if (debugging("sym_boundary"))
    {
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_rel",vel_ghost_rel,dim,"\n");
        fprint_general_vector(stdout,"v_sym",v_sym,dim,"\n");
        printf("temp_reflect = %g\n",temp_reflect);
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
    
    //TODO: Can also try using nonlinearized viscosity for the stencil
    //      i.e. every point of the stencil records its viscosity and
    //      use in the discretization (half indices most likely)
    
    int index = d_index(icoords,top_gmax,dim);
    double mu = m_vst->mu[index] + m_vst->mu_turb[index];
    //double mu = m_vst->mu[index] + field.mu_turb[index];
    
    double tauxx = 2.0/3.0*mu*(2.0*u_x - v_y);
    double tauyy = 2.0/3.0*mu*(2.0*v_y - u_x);
    double tauxy = mu*(u_y + v_x);

    double tauxx_x = 2.0/3.0*mu*(2.0*u_xx - v_xy);
    double tauyy_y = 2.0/3.0*mu*(2.0*v_yy - u_xy);
    double tauxy_y = mu*(u_yy + v_xy);
    double tauxy_x = mu*(u_xy + v_xx);

    v_flux->momn_flux[0] = delta_t*(tauxx_x + tauxy_y);
    v_flux->momn_flux[1] = delta_t*(tauxy_x + tauyy_y);
    
    v_flux->engy_flux = delta_t*(u_x*tauxx + u*tauxx_x + v_x*tauxy
            + v*tauxy_y + u_y*tauxy + u*tauxy_y + v_y*tauyy + v*tauyy_y);

    //if (debugging("no_heatflux")) return;
    
    //double T_x = 0.5*(sten[1][2].temp - sten[1][0].temp)/top_h[0];
    //double T_y = 0.5*(sten[2][1].temp - sten[0][1].temp)/top_h[1];

    double T_xx = (sten[1][2].temp - 2.0*sten[1][1].temp
            + sten[1][0].temp)/sqr(top_h[0]);
    double T_yy = (sten[2][1].temp - 2.0*sten[1][1].temp 
            + sten[0][1].temp)/sqr(top_h[1]);
    
    EOS_PARAMS eos = eqn_params->eos[vsten->comp];
    double gamma = eos.gamma;
    double R_specific = eos.R_specific;
    double Pr = eos.Pr;
    double Pr_turb = eos.Pr_turb;

    double Cp = gamma/(gamma - 1.0)*R_specific;

    double lambda = Cp*m_vst->mu[index]/Pr; //thermal conductivity
    double lambda_turb = Cp*m_vst->mu_turb[index]/Pr_turb; //thermal conductivity

    v_flux->engy_flux += delta_t*(lambda + lambda_turb)*(T_xx + T_yy);
    /////////////////////////////////////////////////////////////////////////////////////
}

//4th order centered difference
void G_CARTESIAN::computeViscousFlux2d_5pt(
        int* icoords,
        SWEEP* m_vst,
        VFLUX* v_flux,
        double delta_t,
        VStencil2d_5pt* vsten)
{
    auto sten = vsten->st;

    double u = sten[2][2].vel[0];
    double v = sten[2][2].vel[1];

    /*
    double u = sten[1][1].vel[0];
    double v = sten[1][1].vel[1];
    */

    double u_x = (-sten[2][4].vel[0] + 8.0*sten[2][3].vel[0] 
            - 8.0*sten[2][1].vel[0] + sten[2][0].vel[0])/12.0/top_h[0];
    
    double u_y = (-sten[4][2].vel[0] + 8.0*sten[3][2].vel[0] 
            - 8.0*sten[1][2].vel[0] + sten[0][2].vel[0])/12.0/top_h[1];
    
    double v_x = (-sten[2][4].vel[1] + 8.0*sten[2][3].vel[1] 
            - 8.0*sten[2][1].vel[1] + sten[2][0].vel[1])/12.0/top_h[0];

    double v_y = (-sten[4][2].vel[1] + 8.0*sten[3][2].vel[1] 
            - 8.0*sten[1][2].vel[1] + sten[0][2].vel[1])/12.0/top_h[1];
    
    /*
    double u_x = 0.5*(sten[1][2].vel[0] - sten[1][0].vel[0])/top_h[0];
    double u_y = 0.5*(sten[2][1].vel[0] - sten[0][1].vel[0])/top_h[1];
    double v_x = 0.5*(sten[1][2].vel[1] - sten[1][0].vel[1])/top_h[0];
    double v_y = 0.5*(sten[2][1].vel[1] - sten[0][1].vel[1])/top_h[1];
    */

    double u_xx = (-sten[2][4].vel[0] + 16.0*sten[2][3].vel[0] - 30.0*sten[2][2].vel[0]
            + 16.0*sten[2][1].vel[0] - sten[2][0].vel[0])/12.0/sqr(top_h[0]);

    double u_yy = (-sten[4][2].vel[0] + 16.0*sten[3][2].vel[0] - 30.0*sten[2][2].vel[0]
            + 16.0*sten[1][2].vel[0] - sten[0][2].vel[0])/12.0/sqr(top_h[1]);

    double v_xx = (-sten[2][4].vel[1] + 16.0*sten[2][3].vel[1] - 30.0*sten[2][2].vel[1]
            + 16.0*sten[2][1].vel[1] - sten[2][0].vel[1])/12.0/sqr(top_h[0]);

    double v_yy = (-sten[4][2].vel[1] + 16.0*sten[3][2].vel[1] - 30.0*sten[2][2].vel[1]
            + 16.0*sten[1][2].vel[1] - sten[0][2].vel[1])/12.0/sqr(top_h[1]);

    /*
    double u_xx = (sten[1][2].vel[0] - 2.0*sten[1][1].vel[0]
            + sten[1][0].vel[0])/sqr(top_h[0]);
    double u_yy = (sten[2][1].vel[0] - 2.0*sten[1][1].vel[0] 
            + sten[0][1].vel[0])/sqr(top_h[1]);
    double v_xx = (sten[1][2].vel[1] - 2.0*sten[1][1].vel[1]
            + sten[1][0].vel[1])/sqr(top_h[0]);
    double v_yy = (sten[2][1].vel[1] - 2.0*sten[1][1].vel[1]
            + sten[0][1].vel[1])/sqr(top_h[1]);
    */
    
    
    double u_xy = 
          (sten[0][3].vel[0] + sten[1][4].vel[0] + sten[3][0].vel[0] + sten[4][1].vel[0])/18.0/top_h[0]/top_h[1]
        - (sten[0][1].vel[0] + sten[1][0].vel[0] + sten[4][3].vel[0] + sten[3][4].vel[0])/18.0/top_h[0]/top_h[1]
        - (sten[0][4].vel[0] + sten[4][0].vel[0] - sten[0][0].vel[0] - sten[4][4].vel[0])/144.0/top_h[0]/top_h[1]
        + (sten[1][1].vel[0] + sten[3][3].vel[0] - sten[1][3].vel[0] - sten[3][1].vel[0])*4.0/9.0/top_h[0]/top_h[1];
    
    double v_xy = 
          (sten[0][3].vel[1] + sten[1][4].vel[1] + sten[3][0].vel[1] + sten[4][1].vel[1])/18.0/top_h[0]/top_h[1]
        - (sten[0][1].vel[1] + sten[1][0].vel[1] + sten[4][3].vel[1] + sten[3][4].vel[1])/18.0/top_h[0]/top_h[1]
        - (sten[0][4].vel[1] + sten[4][0].vel[1] - sten[0][0].vel[1] - sten[4][4].vel[1])/144.0/top_h[0]/top_h[1]
        + (sten[1][1].vel[1] + sten[3][3].vel[1] - sten[1][3].vel[1] - sten[3][1].vel[1])*4.0/9.0/top_h[0]/top_h[1];
    
    
    /*
    double u_xy = 0.25*(sten[2][2].vel[0] - sten[2][0].vel[0]
            - sten[0][2].vel[0] + sten[0][0].vel[0])/top_h[0]/top_h[1];
    double v_xy = 0.25*(sten[2][2].vel[1] - sten[2][0].vel[1]
            - sten[0][2].vel[1] + sten[0][0].vel[1])/top_h[0]/top_h[1];
    */

    int index = d_index(icoords,top_gmax,dim);
    double mu = m_vst->mu[index] + m_vst->mu_turb[index];
    //double mu = m_vst->mu[index] + field.mu_turb[index];
    
    double tauxx = 2.0/3.0*mu*(2.0*u_x - v_y);
    double tauyy = 2.0/3.0*mu*(2.0*v_y - u_x);
    double tauxy = mu*(u_y + v_x);

    double tauxx_x = 2.0/3.0*mu*(2.0*u_xx - v_xy);
    double tauyy_y = 2.0/3.0*mu*(2.0*v_yy - u_xy);
    double tauxy_y = mu*(u_yy + v_xy);
    double tauxy_x = mu*(u_xy + v_xx);
    /*
    double* mu = m_vst->mu;
    double* mu_turb = m_vst->mu_turb;
    int index = d_index(icoords,top_gmax,dim);
    
    double tauxx = 2.0/3.0*mu[index]*(2.0*u_x - v_y);
    double tauyy = 2.0/3.0*mu[index]*(2.0*v_y - u_x);
    double tauxy = mu[index]*(u_y + v_x);

    double tauxx_x = 2.0/3.0*mu[index]*(2.0*u_xx - v_xy);
    double tauyy_y = 2.0/3.0*mu[index]*(2.0*v_yy - u_xy);
    double tauxy_y = mu[index]*(u_yy + v_xy);
    double tauxy_x = mu[index]*(u_xy + v_xx);
    */

    v_flux->momn_flux[0] = delta_t*(tauxx_x + tauxy_y);
    v_flux->momn_flux[1] = delta_t*(tauxy_x + tauyy_y);
    
    v_flux->engy_flux = delta_t*(u_x*tauxx + u*tauxx_x + v_x*tauxy
            + v*tauxy_y + u_y*tauxy + u*tauxy_y + v_y*tauyy + v*tauyy_y);

    
    //if (debugging("no_heatflux")) return;
    /////////////////////////////////////////////////////////////////////////////////////
    
    /*
    double T_x = 0.5*(sten[1][2].temp - sten[1][0].temp)/top_h[0];
    double T_y = 0.5*(sten[2][1].temp - sten[0][1].temp)/top_h[1];
    */

    double T_xx = (sten[1][2].temp - 2.0*sten[1][1].temp
            + sten[1][0].temp)/sqr(top_h[0]);
    double T_yy = (sten[2][1].temp - 2.0*sten[1][1].temp 
            + sten[0][1].temp)/sqr(top_h[1]);
    
    EOS_PARAMS eos = eqn_params->eos[vsten->comp];
    double gamma = eos.gamma;
    double R_specific = eos.R_specific;
    double Pr = eos.Pr;
    double Pr_turb = eos.Pr_turb;

    double Cp = gamma/(gamma - 1.0)*R_specific;
    
    double lambda = Cp*m_vst->mu[index]/Pr; //thermal conductivity
    double lambda_turb = Cp*m_vst->mu_turb[index]/Pr_turb; //thermal conductivity

    v_flux->engy_flux += delta_t*(lambda + lambda_turb)*(T_xx + T_yy);
    //v_flux->engy_flux += delta_t*lambda*(T_xx + T_yy);

    /////////////////////////////////////////////////////////////////////////////////////
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

    
    int index = d_index(icoords,top_gmax,dim);
    double mu = m_vst->mu[index] + m_vst->mu_turb[index];
    //double mu = m_vst->mu[index] + field.mu_turb[index];
    
    double tauxx = 2.0/3.0*mu*(2.0*u_x - v_y - w_z);
    double tauyy = 2.0/3.0*mu*(2.0*v_y - u_x - w_z);
    double tauzz = 2.0/3.0*mu*(2.0*w_z - u_x - v_y);
    
    double tauxy = mu*(u_y + v_x);
    double tauxz = mu*(u_z + w_x);
    double tauyz = mu*(v_z + w_y);

    double tauxx_x = 2.0/3.0*mu*(2.0*u_xx - v_xy - w_xz);
    double tauyy_y = 2.0/3.0*mu*(2.0*v_yy - u_xy - w_yz);
    double tauzz_z = 2.0/3.0*mu*(2.0*w_zz - u_xz - v_yz);

    double tauxy_x = mu*(u_xy + v_xx);
    double tauxz_x = mu*(u_xz + w_xx);
    double tauxy_y = mu*(u_yy + v_xy);
    double tauyz_y = mu*(v_yz + w_yy);
    double tauxz_z = mu*(u_zz + w_xz);
    double tauyz_z = mu*(v_zz + w_yz);

    v_flux->momn_flux[0] = delta_t*(tauxx_x + tauxy_y + tauxz_z);
    v_flux->momn_flux[1] = delta_t*(tauxy_x + tauyy_y + tauyz_z);
    v_flux->momn_flux[2] = delta_t*(tauxz_x + tauyz_y + tauzz_z);
    
    v_flux->engy_flux = delta_t*(u_x*tauxx + u*tauxx_x + v_x*tauxy
            + v*tauxy_y + w_x*tauxz + w*tauxz_x + u_y*tauxy + u*tauxy_y
            + v_y*tauyy + v*tauyy_y + w_y*tauyz + w*tauyz_y + u_z*tauxz
            + u*tauxz_z + v_z*tauyz + v*tauyz_z + w_z*tauzz + w*tauzz_z);


    //if (debugging("no_heatflux")) return;
    /////////////////////////////////////////////////////////////////////////////////////
    
    /*
    double T_x = 0.5*(sten[1][1][2].temp - sten[1][1][0].temp)/top_h[0];
    double T_y = 0.5*(sten[1][2][1].temp - sten[1][0][1].temp)/top_h[1];
    double T_z = 0.5*(sten[2][1][1].temp - sten[0][1][1].temp)/top_h[2];
    */

    double T_xx = (sten[1][1][2].temp - 2.0*sten[1][1][1].temp
            + sten[1][1][0].temp)/sqr(top_h[0]);
    double T_yy = (sten[1][2][1].temp - 2.0*sten[1][1][1].temp
            + sten[1][0][1].temp)/sqr(top_h[1]);
    double T_zz = (sten[2][1][1].temp - 2.0*sten[1][1][1].temp
            + sten[0][1][1].temp)/sqr(top_h[2]);
    
    EOS_PARAMS eos = eqn_params->eos[vsten->comp];
    double gamma = eos.gamma;
    double R_specific = eos.R_specific;
    double Pr = eos.Pr;
    double Pr_turb = eos.Pr_turb;

    double Cp = gamma/(gamma - 1.0)*R_specific;
    
    double lambda = Cp*m_vst->mu[index]/Pr; //thermal conductivity
    double lambda_turb = Cp*m_vst->mu_turb[index]/Pr_turb; //thermal conductivity

    v_flux->engy_flux += delta_t*(lambda + lambda_turb)*(T_xx + T_yy + T_zz);
    /////////////////////////////////////////////////////////////////////////////////////
}


