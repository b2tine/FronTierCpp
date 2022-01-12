#include "cFluid.h"

//Corresponds to the locations of fluid solver cell centers
void G_CARTESIAN::writeMeshFileVTK()
{
    if (pp_mynode() != 0) return;

    RECT_GRID grid = front->pp_grid->Global_grid;
    int* gmax = grid.gmax;
    double* L = grid.L;
    double* h = grid.h;

    char dirname[250];
    sprintf(dirname,"%s/vtk",OutName(front));

    if (!create_directory(dirname,YES))
    {
        screen("Cannot create directory %s\n",dirname);
        LOC(); clean_up(ERROR);
    }

    char mesh_name[250];
    sprintf(mesh_name,"%s/top_grid.vtk",dirname);
    FILE* file = fopen(mesh_name,"w");

    fprintf(file,"# vtk DataFile Version 3.0\n"
            "Topological Grid\n"
            "ASCII\n"
            "DATASET RECTILINEAR_GRID\n");

    int NX = (gmax[0] > 0) ? gmax[0] : 1;
    int NY = (gmax[1] > 0) ? gmax[1] : 1;
    int NZ = (gmax[2] > 0) ? gmax[2] : 1;

    fprintf(file,"DIMENSIONS %d %d %d\n",NX,NY,NZ);

    fprintf(file,"X_COORDINATES %d float\n",NX);
    for (int i = 0; i < NX; ++i)
    {
        fprintf(file,"%f ", L[0] + (i + 0.5)*h[0]);
    }
    fprintf(file,"\n");

    fprintf(file,"Y_COORDINATES %d float\n",NY);
    for (int j = 0; j < NY; ++j)
    {
        fprintf(file,"%f ", L[1] + (j + 0.5)*h[1]);
    }
    fprintf(file,"\n");

    fprintf(file,"Z_COORDINATES %d float\n",NZ);
    for (int k = 0; k < NZ; ++k)
    {
        fprintf(file,"%f ", L[2] + (k + 0.5)*h[2]);
    }

    fclose(file);
}

void G_CARTESIAN::writeCompGridMeshFileVTK()
{
    if (pp_mynode() != 0) return;

    RECT_GRID grid = front->pp_grid->Global_grid;
    int* gmax = grid.gmax;
    double* L = grid.L;
    double* h = grid.h;

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
    sprintf(mesh_name,"%s/comp_grid.vtk",dirname);
    FILE* file = fopen(mesh_name,"w");

    fprintf(file,"# vtk DataFile Version 3.0\n"
            "Computational Grid\n"
            "ASCII\n"
            "DATASET RECTILINEAR_GRID\n");

    int NX = gmax[0] + 1;
    int NY = gmax[1] + 1;
    int NZ = gmax[2] + 1;

    fprintf(file,"DIMENSIONS %d %d %d\n",NX,NY,NZ);

    fprintf(file,"X_COORDINATES %d float\n",NX);
    for (int i = 0; i < NX; ++i)
    {
        fprintf(file,"%f ", L[0] + i*h[0]);
    }
    fprintf(file,"\n");

    fprintf(file,"Y_COORDINATES %d float\n",NY);
    for (int j = 0; j < NY; ++j)
    {
        fprintf(file,"%f ", L[1] + j*h[1]);
    }
    fprintf(file,"\n");

    fprintf(file,"Z_COORDINATES %d float\n",NZ);
    for (int k = 0; k < NZ; ++k)
    {
        fprintf(file,"%f ", L[2] + k*h[2]);
    }

    fclose(file);
}

//TODO: Write parallel version.
//      Follow the pattern used in vtk_plot_scalar_field() in src/front/fprint.c
//      so the components can be observed over the course of the run.
void G_CARTESIAN::writeMeshComponentsVTK()
{
    if (pp_numnodes() > 1)
    {
        printf("\n\tWARNING writeMeshComponentsVTK() not implemented in parallel yet\n\n");
        return;
    }

    serialWriteMeshComponentsVTK();
}

void G_CARTESIAN::serialWriteMeshFileVTK()
{
    if (pp_numnodes() > 1)
    {
        printf("\n\tWARNING: writeMeshFileVTK() not implemented in parallel yet\n\n");
        return;
    }

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

void G_CARTESIAN::serialWriteCompGridMeshFileVTK()
{
    if (pp_numnodes() > 1)
    {
        printf("\n\tWARNING writeMeshFileVTK() not implemented in parallel yet\n\n");
        return;
    }

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
    sprintf(mesh_name,"%s/comp_grid.vtk",dirname);
    FILE* file = fopen(mesh_name,"w");

    fprintf(file,"# vtk DataFile Version 3.0\n"
            "Computational Grid\n"
            "ASCII\n"
            "DATASET RECTILINEAR_GRID\n");

    RECT_GRID* comp_grid = computational_grid(front->interf);
    int* ctop_gmax = comp_grid->gmax;
    double* ctop_L = comp_grid->L;
    double* ctop_U = comp_grid->U;
    double* ctop_h = comp_grid->h;

    int NX = ctop_gmax[0] + 1;
    int NY = ctop_gmax[1] + 1;
    int NZ = ctop_gmax[2] + 1;

    fprintf(file,"DIMENSIONS %d %d %d\n",NX,NY,NZ);

    fprintf(file,"X_COORDINATES %d float\n",NX);
    for (int i = 0; i <= ctop_gmax[0]; ++i)
    {
        fprintf(file,"%f ",ctop_L[0] + ctop_h[0]*i);
    }
    fprintf(file,"\n");

    fprintf(file,"Y_COORDINATES %d float\n",NY);
    for (int j = 0; j <= ctop_gmax[1]; ++j)
    {
        fprintf(file,"%f ",ctop_L[1] + ctop_h[1]*j);
    }
    fprintf(file,"\n");

    fprintf(file,"Z_COORDINATES %d float\n",NZ);
    for (int k = 0; k <= ctop_gmax[2]; ++k)
    {
        fprintf(file,"%f ",ctop_L[2] + ctop_h[2]*k);
    }

    fclose(file);
}

//NOTE: Using POINT_DATA corresponding to nodes of the topological grid
//      correctly places the components at the center of the computational
//      grid cells since the two grids are dual to each other.
//      By overlaying the POINT_DATA on the computational grid during visualization
//      and increasing the point size, the computational grid cells will appear
//      to be colored by their component as intended.

//TODO: Write parallel version
void G_CARTESIAN::serialWriteMeshComponentsVTK()
{
    if (pp_numnodes() > 1)
    {
        printf("\n\tWARNING writeMeshFileVTK() not implemented in parallel yet\n\n");
        return;
    }

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
    sprintf(mesh_name,"%s/grid_comps-%d.vtk",dirname,front->step);
    FILE* file = fopen(mesh_name,"w");

    fprintf(file,"# vtk DataFile Version 3.0\n"
            "Grid Components\n"
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
    fprintf(file,"\n");

    int num_points = 1;
    for (int i = 0; i < dim; ++i)
    {
        num_points *= top_gmax[i] + 1;
    }

    fprintf(file,"POINT_DATA %d\n",num_points);
    fprintf(file,"SCALARS comp int 1\n");
    fprintf(file,"LOOKUP_TABLE component\n");

    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file,"%d\n",cell_center[i].comp);
            //fprintf(file,"%d\n",top_comp[i]);
    }
    
    fclose(file);
}

