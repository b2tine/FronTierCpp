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
                computeViscousFlux(icoords,m_vst,&v_flux,delta_t);
            }
            break;
        }
    /*case 3:
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
                //VFLUX v_flux;
                //computeViscousFlux(icoords,m_vst,&v_flux,delta_t);
            }
            break;
        }*/
    default:
        {
            printf("addViscousFlux() ERROR: Invalid Dimension\n");
            LOC(); clean_up(EXIT_FAILURE);
        }
    }

    LOC(); clean_up(0);
}

void G_CARTESIAN::fillViscousFluxStencil2d(
        int* icoords,
        SWEEP* m_vst,
        VStencil2d* vsten)
{
    /*
    POINTER intfc_state;
    double nip_coords[MAXD];
    double intrp_coeffs[MAXD];
    HYPER_SURF_ELEMENT* hse;
    HYPER_SURF* hs;
    
    int crx_status;

    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };
    */

    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];

    for (int jj = 0; jj < 3; ++jj)
    for (int ii = 0; ii < 3; ++ii)
    {
        VSWEEP* vs = &vsten->st[ii][jj];
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
            //vs->mu = m_vst->mu;;
            //vs->temp = m_vst->temp;
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
                //TODO: generate a ghost state
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

/*
void G_CARTESIAN::fillViscousFluxStencil3d(
        int* icoords,
        SWEEP* m_vst,
        VStencil3d* vsten)
{
}
*/

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
    //TODO: need to add viscosity to STATE?
    //vs->mu = state->mu;
    //vs->temp = state->temp;

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
            double ns[MAXD] = {0.0};
            double ne[MAXD] = {0.0};
            
            normal(Bond_of_hse(hse)->start,hse,hs,ns,front);
            normal(Bond_of_hse(hse)->end,hse,hs,ne,front);

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

            for (int i = 0; i < dim; ++i)
            {
                nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
            }
        }
        break;

    case 3:
        {
            TRI* nearTri = Tri_of_hse(hse);
            
            STATE* st[3];

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

//TODO: See gvisc.c : g_ns_soln()
void G_CARTESIAN::computeViscousFlux(
        int* icoords,
        SWEEP* m_vst,
        VFLUX* v_flux,
        double delta_t)
{
    LOC(); clean_up(0);

    /*
    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];
    if (!gas_comp(comp)) return;

    POINTER state;
    HYPER_SURF* hs;
        //HYPER_SURF_ELEMENT* hse;
    double crx_coords[MAXD];
    int crx_status;

    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    
    int ic[MAXD];
    VStencil3d vsten;
    VSWEEP vstate = vsten.st[1][1][1];//center of stencil
    
    for (int l = 0; l < dim; ++l)
    {
        vstate.icoords[l] = icoords[l];
        vstate.vel[l] = m_vst->momn[l][index]/m_vst->dens[index];
    }

    //From grid point i,j,k need to check for interface crossings
    //between immediate neighbors. Will then need to move to each
    //immediate neighbor and check for interface crossings between
    //its neighbors involved in computing the mixed derivatives.
    
    //NOTE: Not possible to precompute every state ahead of time,
    //      must check crossing and generate ghost states when computing
    //      a specific derivative -- a single point that is involved
    //      in the computation of two separate derivatives may have different
    //      states in the two computations...

    //TODO: The non-mixed second derivatives can be computed with
    //      a loop like this one, and the mixed ones can be dealt with
    //      separately.
    
    double u_nb[2], v_nb[2], w_nb[2], mu_nb[2];

    for (int l = 0; l < dim; ++l)
    {
        for (int nb = 0; nb < 2; ++nb)
        {
            for (int j = 0; j < dim; ++j)
                ic[j] = icoords[j];
            ic[l] = (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;
            
            vstate = vsten.st[ic[0]][ic[1]][ic[2]];
    
        }

    }

//TODO: For mixed partials write functions that accept an icoords array and
//      a direction. Then can compute individual ghost states per derivative.
    */

}


