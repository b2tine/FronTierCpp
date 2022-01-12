#include "iFluid.h"


/*
void Incompress_Solver_Smooth_Basis::vtkPlotSurfacePressure()
{
    if (pp_numnodes() > 1)
    {
        printf("\n\tWARNING vtkPlotSurfacePressure() \
                not implemented in parallel\n\n");
        return;
    }

    INTERFACE *intfc = front->interf;
    SURFACE **s;
    TRI *tri;
    POINT *p;

    FILE *vfile;
    char dirname[200], fname[200];

    sprintf(dirname,"%s/%s%s",outname,"vtk.ts",right_flush(front->step,7));
    if (!create_directory(dirname,NO))
    {
        printf("Cannot create directory %s\n",dirname);
        clean_up(ERROR);
    }
    sprintf(fname,"%s/%s",dirname,"SURFACE_PRESSURE.vtk");

    vfile = fopen(fname,"w");
    fprintf(vfile,"# vtk DataFile Version 2.0\n");
    fprintf(vfile,"Surface stress\n");
    fprintf(vfile,"ASCII\n");
    fprintf(vfile,"DATASET UNSTRUCTURED_GRID\n");

    int num_tri = 0;
    int num_pts = 0;

    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        unsort_surf_point(*s);
        surf_tri_loop(*s,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                if (sorted(p)) continue;
                sorted(p) = YES;
                num_pts++;
            }
            //num_tri++;
        }
        num_tri += (*s)->num_tri;
    }

    //NOTE: sorted(p) and Index_of_point(p) refer to variables
    //      belonging to the same union data structure.
    //      Therefore setting one overwrites the other.

    fprintf(vfile,"POINTS %d float\n", num_pts);

    int n = 0;
    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        surf_tri_loop(*s,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                Index_of_point(Point_of_tri(tri)[i]) = -1;
            }
        }

        surf_tri_loop(*s,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                if (Index_of_point(p) == -1)
                {
                    fprintf(vfile,"%f %f %f\n",Coords(p)[0],Coords(p)[1],Coords(p)[2]);
                    Index_of_point(p) = n++;
                }
            }
        }
    }

    fprintf(vfile,"CELLS %d %d\n",num_tri,4*num_tri);

    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        surf_tri_loop(*s,tri)
        {
            fprintf(vfile,"3 %d %d %d\n",
                    Index_of_point(Point_of_tri(tri)[0]),
                    Index_of_point(Point_of_tri(tri)[1]),
                    Index_of_point(Point_of_tri(tri)[2]));
        }
    }

    fprintf(vfile, "CELL_TYPES %i\n",num_tri);

    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        surf_tri_loop(*s,tri)
        {
            fprintf(vfile,"5\n");
        }
    }

    fprintf(vfile, "CELL_DATA %i\n", num_tri);
    fprintf(vfile, "SCALARS SURFACE_PRESSURE float 1\n");
    fprintf(vfile, "LOOKUP_TABLE default\n");

    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        surf_tri_loop(*s,tri)
        {
            //TODO: Write tri pressure differental on canopy.
            //      See force_on_hse3d() in iFluid for idea of how
            //      to compute the pressure differential acting on
            //      the face of a triangle.
            //
            //fprintf(vfile,"%f\n",tri->color);
        }
    }
}
*/


//FT_GetStatesAtPoint(point,hse,hs,&sl,&sr);
//if (pos_side)
//    *pres += getStatePres(sr);
//else
//    *pres += getStatePres(sl);



//TODO: probably better to take the grid as input so there is no
//      confusion as to what grid we are plotting...

//TODO: binary format

//Corresponds to the locations of fluid solver cell centers
void Incompress_Solver_Smooth_Basis::writeMeshFileVTK()
{
    if (pp_numnodes() > 1)
    {
        //TODO: See cFluid/cFvisualization.cpp
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

void Incompress_Solver_Smooth_Basis::writeCompGridMeshFileVTK()
{
    if (pp_numnodes() > 1)
    {
        //TODO: See cFluid/cFvisualization.cpp
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

/*
void Incompress_Solver_Smooth_Basis::writeGridComponentsVTK()
{
    if (pp_numnodes() > 1)
    {
        printf("\n\tWARNING writeGridComponentsVTK() not implemented in parallel yet\n\n");
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
    sprintf(mesh_name,"%s/mesh_components.vtk",dirname);
    FILE* file = fopen(mesh_name,"w");

    //TODO: Write vtk file
    //
    for (icoords[0] = 1; icoords[0] < top_gmax[0]; ++icoords[0])
    for (icoords[1] = 1; icoords[1] < top_gmax[1]; ++icoords[1])
    for (icoords[2] = 1; icoords[2] < top_gmax[2]; ++icoords[2])
    {
        index = d_index(icoords,top_gmax,dim);
        for (i = 0; i < dim; ++i)
            coords[i] = L[i] + icoords[i]*h[i];

        int comp = top_comp[index]; //print this
    }


    fclose(file);
}
*/


