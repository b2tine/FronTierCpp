#include <iostream>
#include <vector>

#define d_index2d(ix,iy,gmax) \
    ((iy)*((gmax)[0] + 1) + (ix))

struct VENTRY
{
    int icoords[2];
    double vel[2];
};

struct VDATA
{
    int tstep;
    std::vector<VENTRY> data;
};

int main()
{
    int top_gmax[2] = {3,3};
    int maxidx = (top_gmax[0] + 1)*(top_gmax[1] + 1);

    double vel0[2][maxidx] = {
        {0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0},
        {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3}
    };


    std::vector<VDATA> velmat;

    VDATA veldata0;
    veldata0.tstep = 0;
    veldata0.data.reserve(50);

    for (int i = 0; i <= top_gmax[0]; ++i)
    for (int j = 0; j <= top_gmax[1]; ++j)
    {
        int index = d_index2d(i,j,top_gmax);
        printf("%g %g",vel0[0][index],vel0[1][index]);
        printf("\t (%d,%d) index = %d\n",i,j,index);
        VENTRY ventry = {i,j,vel0[0][index],vel0[1][index]};
        veldata0.data.push_back(ventry);
    }
    printf("\n\n");

    /*
    for (auto it : veldata0.data)
    {
        auto ic = it.icoords;
        auto iv = it.vel;
        printf("%g %g",iv[0],iv[1]);
        printf("\t (%d,%d) index = %d\n",
            ic[0],ic[1],d_index2d(ic[0],ic[1],top_gmax));
    }
    printf("\n\n");
    */

    velmat.push_back(veldata0);


    double vel1[2][maxidx] = {
        {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3},
        {0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0}
    };

    VDATA veldata1;
    veldata1.tstep = 1;
    veldata1.data.reserve(50);

    for (int i = 0; i <= top_gmax[0]; ++i)
    for (int j = 0; j <= top_gmax[1]; ++j)
    {
        int index = d_index2d(i,j,top_gmax);
        printf("%g %g",vel1[0][index],vel1[1][index]);
        printf("\t (%d,%d) index = %d\n",i,j,index);
        VENTRY ventry = {i,j,vel1[0][index],vel1[1][index]};
        veldata1.data.push_back(ventry);
    }
    printf("\n\n");

    velmat.push_back(veldata1);



    int T = 2;
    for (int t = 0; t < T; ++t)
    {
        int tstep = velmat[t].tstep;
        printf("\n\ttstep = %d\n",tstep);
        
        for (auto it : velmat[t].data)
        {
            auto iv = it.vel;
            auto ic = it.icoords;
            printf("%g %g",iv[0],iv[1]);
            printf("\t (%d,%d) index = %d\n",
                ic[0],ic[1],d_index2d(ic[0],ic[1],top_gmax));
        }
    }
    
    return 0;
}
