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

//TODO: Write and use constructors for xENTRYXd objects

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

IDATA2d Incompress_Solver_Smooth_Basis::getIntfcData2d()
{
    IDATA2d idata;
    idata.tstep = front->step;
    idata.dt = front->dt;
    idata.time = front->time;

    INTERFACE* intfc = front->interf;
    CURVE** c;
    BOND* b;

    intfc_curve_loop(intfc,c)
    {
        if (is_bdry(*c)) continue;

        curve_bond_loop(*c,b)
        {
            POINT* p = b->end;
            double* coords = Coords(p);
            double* vel = p->vel;
            STATE* sl = (STATE*)left_state(p);
            double vort = sl->vort;
            
            IENTRY2d ientry = {coords[0],coords[1],vel[0],vel[1],vort};
           
            /* 
            IENTRY2d ientry;
            for (int l = 0; l < dim; ++l)
            {
                ientry.coords[l] = coords[l];
                ientry.vel[l] = vel[l];
            }
            ientry.vort = vort;
            */

            idata.data.push_back(ientry);
        }
    }

    return idata;
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

IDATA3d Incompress_Solver_Smooth_Basis::getIntfcData3d()
{
    IDATA3d idata;
    idata.tstep = front->step;
    idata.dt = front->dt;
    idata.time = front->time;

    //NOTE: Starting with just strings/curves, add surfaces later
    INTERFACE* intfc = front->interf;
    CURVE** c;
    BOND* b;

    intfc_curve_loop(intfc,c)
    {
        //TODO: check if string hsbdry
        if (is_bdry(*c)) continue;

        curve_bond_loop(*c,b)
        {
            POINT* p = b->end;
            double* coords = Coords(p);
            double* vel = p->vel;
            STATE* sl = (STATE*)left_state(p);
                //double vort = sl->vort;
            
            IENTRY3d ientry = {
                coords[0],coords[1],coords[2],
                vel[0],vel[1],vel[2]};
           
            /* 
            IENTRY3d ientry;
            for (int l = 0; l < dim; ++l)
            {
                ientry.coords[l] = coords[l];
                ientry.vel[l] = vel[l];
                    //ientry.vort[l] = vort[l];
            }
            */

            idata.data.push_back(ientry);
        }
    }

    return idata;
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

