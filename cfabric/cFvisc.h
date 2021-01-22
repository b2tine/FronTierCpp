#ifndef CFVISC_H
#define CFVISC_H

#include "cFluid.h"


//TODO: Probably won't need these stencil functions
/*
struct Stencil_3d
{
    int**** icoords;
    //STATE**** sts3d;
    //Locstate*** sts3d;
};


void initStencil(Stencil_3d* stencil);
void fillStencil(Stencil_3d* stencil);
*/


struct VSWEEP
{
    int icoords[MAXD];
    double vel[MAXD];
    double mu;
    double temp;
};

struct VStencil3d
{
    VSWEEP st[MAXD][MAXD][MAXD]; //[1][1][1] is center of stencil (i,j,k)
};

#endif
