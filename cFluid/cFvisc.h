#ifndef CFVISC_H
#define CFVISC_H


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
    double dens;
};

struct VStencil2d
{
    COMPONENT comp;
    VSWEEP st[3][3]; //[1][1] is center of stencil (i,j)
};

struct VStencil2d_5pt
{
    COMPONENT comp;
    VSWEEP st[5][5]; //[1][1] is center of stencil (i,j)
};

struct VStencil3d
{
    COMPONENT comp;
    VSWEEP st[3][3][3]; //[1][1][1] is center of stencil (i,j,k)
};


#endif
