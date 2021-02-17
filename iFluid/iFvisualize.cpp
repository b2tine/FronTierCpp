#include "iFluid.h"

//TODO: probably better to take the grid as input so there is no
//      confusion as to what grid we are plotting...

//TODO: binary format

//Corresponds to the locations of fluid solver cell centers
void Incompress_Solver_Smooth_Basis::writeMeshFileVTK()
{
    char dirname[250];
    sprintf(dirname,"%s/vtk",OutName(front));

    if (pp_mynode() == 0)
    {
        if (!create_directory(dirname,YES))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }
    }

    char mesh_name[250];
    sprintf(mesh_name,"%s/top_grid.vtk",dirname);
    FILE* file = fopen(mesh_name,"w");

    fprintf(file,"# vtk DataFile Version 3.0\n"
            "Topological Grid\n"
            "ASCII\n"
            "DATASET RECTILINEAR_GRID\n");

    int NX = top_gmax[0] + 1;
    int NY = top_gmax[1] + 1;
    int NZ = top_gmax[2] + 1;

    fprintf(file,"DIMENSIONS %d %d %d\n",NX,NY,NZ);

    fprintf(file,"X_COORDINATES %d float\n",NX);
    for (int i = 0; i <= top_gmax[0]; ++i)
    {
        fprintf(file,"%f ",top_L[0] + top_h[0]*i);
    }
    fprintf(file,"\n");

    fprintf(file,"Y_COORDINATES %d float\n",NY);
    for (int j = 0; j <= top_gmax[1]; ++j)
    {
        fprintf(file,"%f ",top_L[1] + top_h[1]*j);
    }
    fprintf(file,"\n");

    fprintf(file,"Z_COORDINATES %d float\n",NZ);
    for (int k = 0; k <= top_gmax[2]; ++k)
    {
        fprintf(file,"%f ",top_L[2] + top_h[2]*k);
    }

    fclose(file);
}
