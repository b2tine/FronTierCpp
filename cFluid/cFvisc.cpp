#include "cFluid.h"
//#include "cFvisc.h"

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

    for (int jj = 0; jj < 3; ++jj)
    for (int ii = 0; ii < 3; ++ii)
    {
        VSWEEP* vs = &vsten->st[jj][ii];
        vs->icoords[0] = icoords[0] + ii - 1;
        vs->icoords[1] = icoords[1] + jj - 1;
        int idx_nb = d_index(vs->icoords,top_gmax,dim);
        vs->comp = top_comp[idx_nb];
        if (vs->comp != comp)
        {
            setViscousGhostState(comp,vs,m_vst);
        }
        else
        {
            for (int i = 0; i < dim; ++i)
            {
                vs->vel[i] = m_vst->momn[i][idx_nb]/m_vst->dens[idx_nb];
            }
        }
    }

    /*
    //TESTING
    bool needghost = false;
    for (int jj = 0; jj < 3; ++jj)
    {
        for (int ii = 0; ii < 3; ++ii)
        {
            VSWEEP vs = vsten->st[ii][jj];
            if (vs.comp != comp)
            {
                needghost = true;
            }
        }
    }

    if (needghost)
    {
        printf("\nicoords = (%d,%d,%d) index = %d\n\n",
                icoords[0],icoords[1],icoords[2],index);
        
        for (int jj = 2; jj >= 0; --jj)
        {
            for (int ii = 0; ii < 3; ++ii)
            {
                VSWEEP vs = vsten->st[ii][jj];
                printf("  %d  %s",vs.comp,
                        ((ii+1) % 3 == 0) ? "" : "|");
            }
            printf("\n-----------------\n");
        }
    }

    return;
    ////////////////////////////////////////////////////////
    */

    //OLD DEBUG INFO
    /*
    int icx0[MAXD], icx1[MAXD];
    int icy0[MAXD], icy1[MAXD];
    int icx0y0[MAXD], icx1y0[MAXD];
    int icx0y1[MAXD], icx1y1[MAXD];
    for (int i = 0; i < dim; ++i)
    {
        icx0[i] = icoords[i];
        icx1[i] = icoords[i];
        icy0[i] = icoords[i];
        icy1[i] = icoords[i];
        icx0y0[i] = icoords[i];
        icx1y0[i] = icoords[i];
        icx0y1[i] = icoords[i];
        icx1y1[i] = icoords[i];
    }

    icx0[0] = icoords[0] - 1;
    icx1[0] = icoords[0] + 1;
    icy0[1] = icoords[1] - 1;
    icy1[1] = icoords[1] + 1;
    
    icx0y0[0] = icoords[0] - 1; icx0y0[1] = icoords[1] - 1;
    icx1y0[0] = icoords[0] + 1; icx1y0[1] = icoords[1] - 1;
    icx0y1[0] = icoords[0] - 1; icx0y1[1] = icoords[1] + 1;
    icx1y1[0] = icoords[0] + 1; icx1y1[1] = icoords[1] + 1;

    int index_x0 = d_index(icx0,top_gmax,dim);
    int index_x1 = d_index(icx1,top_gmax,dim);
    int index_y0 = d_index(icy0,top_gmax,dim);
    int index_y1 = d_index(icy1,top_gmax,dim);
    
    int index_x0y0 = d_index(icx0y0,top_gmax,dim);
    int index_x1y0 = d_index(icx1y0,top_gmax,dim);
    int index_x0y1 = d_index(icx0y1,top_gmax,dim);
    int index_x1y1 = d_index(icx1y1,top_gmax,dim);

    COMPONENT comp_x0 = top_comp[index_x0];
    COMPONENT comp_x1 = top_comp[index_x1];
    COMPONENT comp_y0 = top_comp[index_y0];
    COMPONENT comp_y1 = top_comp[index_y1];
    
    COMPONENT comp_x0y0 = top_comp[index_x0y0];
    COMPONENT comp_x1y0 = top_comp[index_x1y0];
    COMPONENT comp_x0y1 = top_comp[index_x0y1];
    COMPONENT comp_x1y1 = top_comp[index_x1y1];

    if (comp_x0 != comp || comp_x1 != comp ||
        comp_y0 != comp || comp_y1 != comp ||
        comp_x0y0 != comp || comp_x1y0 != comp ||
        comp_x0y1 != comp || comp_x1y1 != comp)
    {
        printf("\nicoords = (%d,%d,%d) index = %d\n\n",
                icoords[0],icoords[1],icoords[2],index);

        printf("  %d  |  %d  |  %d  \n",comp_x0y1,comp_y1,comp_x1y1);
        printf("--------------------\n");
        printf("  %d  |  %d  |  %d  \n",comp_x0,comp,comp_x1);
        printf("--------------------\n");
        printf("  %d  |  %d  |  %d  \n",comp_x0y0,comp_y0,comp_x1y0);
    }
    */
}

void G_CARTESIAN::fillViscousFluxStencil3d(
        int* icoords,
        SWEEP* m_vst,
        VStencil3d* vsten)
{
    printf("ERROR: fillViscousFluxStencil3d() not implemented yet!\n");
    LOC(); clean_up(EXIT_FAILURE);
}

void G_CARTESIAN::setViscousGhostState(
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
                comp,front->interf,NO_SUBDOMAIN,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);
    
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
            setNeumannViscousGhostState(m_vst,vs,&ghost_coords[0],
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
            printf("setViscousGhostState() ERROR: "
                    "unknown boundary type\n");
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
    
    double coords_reflect[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
        coords_reflect[j] = crx_coords[j] + dist_reflect*nor[j];

    double dens_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,nullptr);
    /*FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,&m_vst->dens[index]);*/

    double mom_reflect[MAXD];
    double vel_reflect[MAXD];
    for (int j = 0; j < dim; ++j)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&mom_reflect[j],nullptr);
        /*FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&momn_reflect[j],&m_vst->momn[j][index]);*/
        vel_reflect[j] = mom_reflect[j]/dens_reflect;
    }

    //TODO: modify vel_reflect -- slip velocity etc. 

    for (int j = 0; j < dim; ++j)
    {
        vs->vel[j] = vel_reflect[j];
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
        VStencil3d* sten)
{
    printf("ERROR: computeViscousFlux3d() not implemented yet!\n");
    LOC(); clean_up(EXIT_FAILURE);

    /*
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
    */
}


