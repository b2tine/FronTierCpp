#include "cFluid.h"
//#include "cFvisc.h"


//TODO: Finish Implementation
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
                //VFLUX v_flux;
                //computeViscousFlux(icoords,m_vst,&v_flux,delta_t);
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
    POINTER state;
    HYPER_SURF_ELEMENT* hse;
    HYPER_SURF* hs;
    double crx_coords[MAXD];
    */
    
    /*
    int crx_status;

    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };
    */

    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];

    ////////////////////////////////////////////////////////
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
            //TODO: generate a ghost state based on reflection
            //      through the nearest_interface_point.
            //
            //      FT_FindNearestIntfcPointInRange(front,comp,
            //          coords_ghost,NO_BOUNDARIES,crx_coords,
            //          intrp_coeffs,&hse,&hs,range);
            
            //setViscousGhostState(vs);
        }
    }

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

//TODO: This function should become the fillStencil() like function
void G_CARTESIAN::computeViscousFlux(
        int* icoords,
        SWEEP* m_vst,
        VFLUX* v_flux,
        double delta_t)
{
    LOC(); clean_up(0);

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
    
            /*
            FT_StateStructAtGridCrossing(front,front->grid_intfc,
                    icoords,dir[l][nb],comp,&intfc_state,&hs,crx_coords);
            */
            
                /*
                FT_StateStructAtGridCrossing2(front,icoords,dir[l][nb],
                        comp,&intfc_state,&hs,&hse,crx_coords);
                */
        }

    }

            /*
            //NOTE icoords replaced with ic in this call
            FT_StateStructAtGridCrossing(front,front->grid_intfc,
                    ic,dir[l][nb],comp,&intfc_state,&hs,crx_coords);
            */


//TODO: For mixed partials write functions that accept an icoords array and
//      a direction. Then can compute individual ghost states per derivative.


}


