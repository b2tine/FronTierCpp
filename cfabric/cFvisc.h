#ifndef CFVISC_H
#define CFVISC_H

//#include "cFluid.h"

struct VFLUX
{
    int icoords[MAXD];
    double momn_flux[MAXD];
    double engy_flux;
    double heat_flux;
};

struct VSWEEP
{
    int icoords[MAXD];
    COMPONENT comp;
    double vel[MAXD];
    double mu;
    double temp;
};

struct VStencil2d
{
    VSWEEP st[MAXD][MAXD]; //[1][1] is center of stencil (i,j)
};

struct VStencil3d
{
    VSWEEP st[MAXD][MAXD][MAXD]; //[1][1][1] is center of stencil (i,j,k)
};


#endif
