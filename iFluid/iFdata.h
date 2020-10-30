#ifndef IFDATA_H
#define IFDATA_H

//For generating PINN training data



//2d
struct VENTRY2d
{
    int icoords[2];
    double vel[2];
    double vort;
};

struct VDATA2d
{
    int tstep;
    double dt;
    double time;
    std::vector<VENTRY2d> data;
};

//3d
struct VENTRY3d
{
    int icoords[3];
    double vel[3];
    double vort[3];
};

struct VDATA3d
{
    int tstep;
    double dt;
    double time;
    std::vector<VENTRY3d> data;
};



#endif
