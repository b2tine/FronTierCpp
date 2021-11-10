#include "vtk.h"


static void mkdirTree(std::string sub, std::string dir);

void vtk_write_pointset(
        const std::vector<POINT*>& points,
        const std::string& fname, const int setID)
{
    std::ofstream outfile(fname+".vtk");
    int npts = points.size();

    //TODO: should we use Version 2.0 ?
	outfile << "# vtk DataFile Version 3.0\n";
    outfile << "annotated point set\n";
    outfile << "ASCII\n";
    outfile << "DATASET POLYDATA\n";

	outfile << "POINTS " << npts << " float\n";
    for (int i = 0; i < npts; ++i)
    {
        double* coords = Coords(points[i]);
        outfile << coords[0] << " "
                << coords[1] << " "
                << coords[2] << "\n";
    }
    
	outfile << "VERTICES " << npts << " " << 2*npts << "\n";
    for (int i = 0; i < npts; ++i)
        outfile << 1 << " " << i << "\n";
    
	outfile << "CELL_DATA " << npts << "\n";
    outfile << "SCALARS collision int 1\n";
    outfile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < npts; ++i)
	    outfile << setID << "\n";
    
    outfile.close();
}	/* end vtk_write_pointset */

/*
void vtk_write_pointset(
        const POINT** points,
        const std::string& fname,
        int setID)
{
    std::ofstream outfile(fname);
    int npts = points.size();

	outfile << "# vtk DataFile Version 3.0\n";
    outfile << "annotated point set\n";
    outfile << "ASCII\n";
    outfile << "DATASET POLYDATA\n";

	outfile << "POINTS " << npts << " float\n";
    for (int i = 0; i < npts; ++i)
    {
        double* coords = Coords(points[i]);
        outfile << coords[0] << " "
                << coords[1] << " "
                << coords[2] << "\n";
    }
    
	outfile << "VERTICES " << npts << " " << 2*npts << "\n";
    for (int i = 0; i < npts; ++i)
        outfile << 1 << " " << i << "\n";
    
	outfile << "CELL_DATA " << npts << "\n";
    outfile << "SCALARS collision int 1\n";
    outfile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < npts; ++i)
	    outfile << setID << "\n";
    
    outfile.close();
}*/	/* end vtk_write_pointset */

void createDirectory(std::string new_dir)
{
    struct stat st;
    int status = stat(new_dir.c_str(), &st);
    if( status != 0 && !S_ISDIR(st.st_mode) )
        mkdirTree(new_dir, "");
}

void mkdirTree(std::string sub, std::string dir)
{
    if(sub.length() == 0)
        return;

    int i = 0;
    for( i; i < sub.length(); i++)
    {
        dir += sub[i];
        if (sub[i] == '/')
            break;
    }

    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if(i+1 < sub.length())
        mkdirTree(sub.substr(i+1), dir);
}

