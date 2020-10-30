#include "iFluid.h"
#include "iFdata.h"

#include <iostream>
#include <vector>




//For generating PINN training data

/*
void Incompress_Solver_Smooth_Basis::writeTimeFile()
{
    char tv_name[100];
    sprintf(tv_name,"%s/time-%d.txt",out_name,(int)timevec.size());
    FILE* tv_file = fopen(tv_name,"w");
    for (int i = 0; i < timevec.size(); ++i)
    {
        fprintf(tv_file,"%20.14f\n",timevec[i]);
    }
}
*/

void Incompress_Solver_Smooth_Basis::writeMeshFile2d()
{
    char mesh_name[250];
    sprintf(mesh_name,"%s/mesh-%d-%d",OutName(front),
            imax,jmax);
    if (pp_numnodes() > 1)
	{
        sprintf(mesh_name,"%s-nd%s",mesh_name,
                right_flush(pp_mynode(),4));
	}
    sprintf(mesh_name,"%s.txt",mesh_name);
    FILE* mesh_file = fopen(mesh_name,"w");

    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    {
        int index = d_index2d(i,j,top_gmax);
        auto coords = cell_center[index].getCoords();
        fprintf(mesh_file,"%20.14f %20.14f\n",coords[0],coords[1]);
            //fprintf(mesh_file,"%20.14f %20.14f",coords[0],coords[1]);
            //fprintf(mesh_file,"\t (%d,%d) index = %d\n",i,j,
            //      d_index2d(i,j,top_gmax));
    }
    fclose(mesh_file);
}

void Incompress_Solver_Smooth_Basis::writeMeshFile3d()
{
    char mesh_name[250];
    sprintf(mesh_name,"%s/mesh-%d-%d-%d",OutName(front),
            imax,jmax,kmax);
    if (pp_numnodes() > 1)
	{
        sprintf(mesh_name,"%s-nd%s",mesh_name,
                right_flush(pp_mynode(),4));
    }
    sprintf(mesh_name,"%s.txt",mesh_name);
    FILE* mesh_file = fopen(mesh_name,"w");

    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    for (int k = kmin; k <= kmax; ++k)
    {
        int index = d_index3d(i,j,k,top_gmax);
        auto coords = cell_center[index].getCoords();
        fprintf(mesh_file,"%20.14f %20.14f %20.14f\n",
                coords[0],coords[1],coords[2]);
    }
    fclose(mesh_file);
}

VDATA2d Incompress_Solver_Smooth_Basis::getVelData2d()
{
    double **vel = field->vel;
    double *vort = field->vort;

    VDATA2d veldata;
    veldata.tstep = front->step;
    veldata.dt = front->dt;
    veldata.time = front->time;
    veldata.data.reserve(imax*jmax);
    
    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    {
        int index  = d_index2d(i,j,top_gmax);
        VENTRY2d ventry = {i,j,vel[0][index],vel[1][index],vort[index]};
        veldata.data.push_back(ventry);
    }
    return veldata;
}

VDATA3d Incompress_Solver_Smooth_Basis::getVelData3d()
{
    double **vel = field->vel;
    double **vort = field->vorticity;

    VDATA3d veldata;
    veldata.tstep = front->step;
    veldata.dt = front->dt;
    veldata.time = front->time;
    veldata.data.reserve(imax*jmax*kmax);
    
    for (int i = imin; i <= imax; ++i)
    for (int j = jmin; j <= jmax; ++j)
    for (int k = kmin; k <= kmax; ++k)
    {
        int index  = d_index3d(i,j,k,top_gmax);
        VENTRY3d ventry = {i,j,k,
            vel[0][index],vel[1][index],vel[2][index],
            vort[0][index],vort[1][index],vort[2][index]};
        veldata.data.push_back(ventry);
    }
    return veldata;
}

//2d
std::vector<int> Incompress_Solver_Smooth_Basis::getMaxIJ()
{
    return std::vector<int>{imax,jmax};
}

//3d
std::vector<int> Incompress_Solver_Smooth_Basis::getMaxIJK()
{
    return std::vector<int>{imax,jmax,kmax};
}

std::vector<int> Incompress_Solver_Smooth_Basis::getTopGMax()
{
    switch(dim)
    {
        case 2:
            return std::vector<int>(top_gmax,top_gmax+1);
        case 3:
            return std::vector<int>(top_gmax,top_gmax+2);
    }
}

/*
std::vector<int> Incompress_Solver_Smooth_Basis::getTopGMax()
{
    return std::vector<int>(top_gmax,top_gmax+1);
}
*/

