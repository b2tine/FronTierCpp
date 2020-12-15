#ifndef IFDATA_H
#define IFDATA_H

//For generating PINN training data

//TODO: Add constructors, at least for the
//      VENTRYxd and IENTRYxd objects. 

//fluid
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

//interface
struct IENTRY2d
{
    double coords[2];
    double vel[2];
    double vort;
};

struct IDATA2d
{
    int tstep;
    double dt;
    double time;
    std::vector<IENTRY2d> data;
};

//fluid
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

//interface
struct IENTRY3d
{
    double coords[3];
    double vel[3];
    //double vort[3];
    //double pres
};

struct IDATA3d
{
    int tstep;
    double dt;
    double time;
    std::vector<IENTRY3d> data;
};


#endif
