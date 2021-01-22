#include "cFvisc.h"



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
    //for (int k = imin[2]; k <= imax[2]; k++)
    for (int j = imin[1]; j <= imax[1]; j++)
    for (int i = imin[0]; i <= imax[0]; i++)
    {
        int index = d_index3d(i,j,k,top_gmax);
        if (!gas_comp(top_comp[index]))
        {
            //zero state vals if necessary
            continue;
        }

        int icoords[MAXD] = {i,j,k};
        computeViscousFlux(icoords,m_vst,m_flux,delta_t);
    }
    
}

//TODO: Finish Implementation
void G_CARTESIAN::computeViscousFlux(
        int* icoords,
        SWEEP* m_vst,
        FSWEEP *m_flux,
        double delta_t)
{
    GRID_DIRECTION ldir[3] = {WEST,SOUTH,LOWER};
    GRID_DIRECTION rdir[3] = {EAST,NORTH,UPPER};

    //From grid point i,j,k need to check for interface crossings
    //between immediate neighbors. Will then need to move to each
    //immediate neighbor and check for interface crossings between
    //its neighbors involved in computing the mixed derivatives.
    
    for (int l = 0; l < dim; ++l)
    {

        for (int nb = 0; nb < dim; ++nb)
        {
    
            //FT_StateStructAtGridCrossing();
        }

    }
}
