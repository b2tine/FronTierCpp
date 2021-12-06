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
    double temp;
    double mu;
    double mu_lam;
    double mu_turb;
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
