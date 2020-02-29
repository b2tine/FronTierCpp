#ifndef VTK_H
#define VTK_H

#include <FronTier.h>

#include <fstream>
#include <vector>
#include <string>


void createDirectory(std::string new_dir);


void vtk_write_pointset(
        const std::vector<POINT*>& points,
        const std::string& fname,
        int setID);

/*
void vtk_write_pointset(
        const POINT** points,
        int npts,
        const std::string& fname,
        int setID);
*/




#endif
