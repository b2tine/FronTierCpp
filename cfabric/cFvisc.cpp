//#include "cFvisc.h"
#include "cFluid.h"



//TODO: Probably don't need these stencil functions
/*
void initStencil(Stencil_3d* stencil)
{
    //FT_QuadArrayMemoryAlloc();
}
*/

/*
void fillStencil(Stencil_3d* stencil)
{
}
*/

//TODO: Finish Implementation
void G_CARTESIAN::addViscousFlux(
        SWEEP* m_vst,
        FSWEEP *m_flux,
        double delta_t)
{
    for (int k = imin[2]; k <= imax[2]; k++)
    for (int j = imin[1]; j <= imax[1]; j++)
    for (int i = imin[0]; i <= imax[0]; i++)
    {   
        int icoords[MAXD] = {i,j,k};
        computeViscousFlux(icoords,m_vst,m_flux,delta_t);
        //TODO: Could have computeViscousFlux() return a
        //      FSWEEP like structure which we then add the 
        //      contents of to m_flux here. More faithful
        //      implementation to the function name, and allows
        //      for single responsibility.
    }
    
}

//TODO: Finish Implementation
void G_CARTESIAN::computeViscousFlux(
        int* icoords,
        SWEEP* m_vst,
        FSWEEP *m_flux,
        double delta_t)
{
    int index = d_index3d(i,j,k,top_gmax);
    COMPONENT comp = top_comp[index];
    if (!gas_comp(comp))
    {
        return;
    }

    POINTER state;
    HYPER_SURF* hs;
        //HYPER_SURF_ELEMENT* hse;
    double crx_coords[MAXD];
    int crx_status;

    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    
    //TODO: create FSWEEP var and alloc data

    //From grid point i,j,k need to check for interface crossings
    //between immediate neighbors. Will then need to move to each
    //immediate neighbor and check for interface crossings between
    //its neighbors involved in computing the mixed derivatives.
    
    for (int l = 0; l < dim; ++l)
    {

        for (int nb = 0; nb < dim; ++nb)
        {
    
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


    //TODO: Add computed flux to m_flux, or return an FSWEEP like structure.

}
