#include "cFluid.h"

//Corresponds to the locations of fluid solver cell centers
void G_CARTESIAN::writeMeshFileVTK()
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

void G_CARTESIAN::writeCompGridMeshFileVTK()
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

void G_CARTESIAN::writeMeshComponentsVTK()
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

    
    /*
    //TODO: Using CELL_DATA to color the grid cells requires mapping
    //      the indices of the topological grid nodes to the corresponding
    //      cell centers of the computational grid.
    //
    //      Using POINT_DATA, as above, correctly places the components
    //      and if the point size is made large they appears as if they
    //      were using the CELL_DATA field of the computational grid
    //        
    //
    ////////////////////////////////////////////////////
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
    ////////////////////////////////////////////////////
    */


    /*
    int num_cells = 1;
    for (int i = 0; i < dim; ++i)
    {
        num_cells *= (top_gmax[i] == 0) ? 1 : top_gmax[i];
        //num_cells *= (ctop_gmax[i] == 0) ? 1 : ctop_gmax[i];
    }

    fprintf(file,"CELL_DATA %d\n",num_cells);
    fprintf(file,"SCALARS cell_scalars int 1\n");
    fprintf(file,"LOOKUP_TABLE component\n");
    */

    /*
    for (int i = 0; i < num_cells; ++i)
    {
        fprintf(file,"%d\n",cell_center[i].comp);
            //fprintf(file,"%d\n",top_comp[i]);
    }

    for (int i = 0; i <= top_gmax[0]; ++i)
    for (int j = 0; j <= top_gmax[1]; ++j)
    for (int k = 0; k <= top_gmax[2]; ++k)
    {
        int icoords[MAXD] = {i,j,k};
        int index = d_index(icoords,top_gmax,dim);
        //fprintf(file,"%d\n",cell_center[index].comp);
        fprintf(file,"%d\n",top_comp[index]);
    }

    fclose(file);
    */
}
