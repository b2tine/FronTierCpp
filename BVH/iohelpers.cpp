#include "BVH.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>


static void mkdirTree(std::string sub, std::string dir)
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


void createDirectory(std::string new_dir)
{
    struct stat st;
    int status = stat(new_dir.c_str(), &st);
    if( status != 0 && !S_ISDIR(st.st_mode) )
        mkdirTree(new_dir, "");
}

void BVH::writeHilbertCurveFile(std::string outdir) const
{
    auto leafvec = getSortedLeafData();
    int num_hse = leafvec.size();
    //assert(!bvMap.empty());
    //int num_hse = bvMap.size();

    outdir += "/";
    std::string geomdir("OOGL/");
    createDirectory(outdir + geomdir);

    //Write geomview OOGL file for the curve
    std::ofstream outfile(outdir + geomdir + "hilbert_curve.vect");
    outfile << "VECT\n\n";
    //num_polylines total_num_vertices num_colors
    outfile << 1 << " " << num_hse << " " << 1 << "\n\n";
    //number vertices in each polyline (would be a list if more than one polyline)
    outfile << num_hse << "\n\n";
    //number of colors in each polyline (would be a list if more than one polyline)
    outfile << 1 << "\n";
    //coordinates of curve vertices (same as the points/boxcenters)
    for( int i = 0; i < num_hse; i++ )
    {
        CGAL_Point p = leafvec[i].first;
        outfile << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    outfile << "\n";
    //red, green, blue, alpha(opacity)
    outfile << 1.0 << " " << 0.0 << " " <<  0.0 << " " <<  1.0;//curve color (red)
    outfile.close();


    //Write geomview OOGL file for the points of the hilbert
    //curve so we can display in a different color
    outfile.open(outdir + geomdir + "hilbert_curve_points.vect");
    outfile << "VECT\n\n";
    //num_polylines total_num_vertices num_colors
    outfile << num_hse << " " << num_hse << " " << 1 << "\n\n";
    //number vertices in each polyline (degenerate polyline is a point)
    for( int i = 0; i < num_hse; i++ )
        outfile << 1 << "\n";
    outfile << "\n";
    //first point is assigned the color (indexed by 1) given below
    outfile << 1 << "\n";
    //remaining points get same color as first point (indicated by 0)
    for( int i = 1; i < num_hse; i++ )
        outfile << 0 << "\n";
    outfile << "\n";
    //coordinates of the curve points
    for( int i = 0; i < num_hse; i++ )
    {
        CGAL_Point p = leafvec[i].first;
        outfile << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    outfile << "\n";
    //red, green, blue, alpha(opacity)
    outfile << 0.0 << " " << 1.0 << " " <<  0.0 << " " <<  1.0;//point color (green)
    outfile.close();


    //write sorted bboxes to geomview OOOGL files
    char geomfile[100];
    for( int i = 0; i < num_hse; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile.open(outdir + geomdir + geomfile);
        outfile << "BBOX\n";

        BV_Point lcoords = leafvec[i].second->getBV().lower;
        //BV_Point lcoords = bvMap[ctrVec[i]]->getBV().lower;
        for( int k = 0; k < 3; k++ )
            outfile << lcoords[k] << " ";
        outfile << "\n";

        BV_Point ucoords = leafvec[i].second->getBV().upper;
        //BV_Point ucoords = bvMap[ctrVec[i]]->getBV().upper;
        for( int k = 0; k < 3; k++ )
            outfile << ucoords[k] << " ";
        outfile.close();
    }


    //OOGL file to display mesh + curve + points
    outfile.open(outdir + "hcurve.list");
    outfile << "LIST\n\n";

    //input mesh
    outfile << "{\n" << "appearance {\n\t+edge +face +transparent\n\t";
    outfile << "material {\n\t\tedgecolor 0 0 1\n\t\t";
    outfile << "ambient 0 0 0.5\n\t\talpha 0.45\n \t\t}\n\n}";
    outfile << "\n\n{ < input-mesh.off }\n\n}\n\n";

    //hilbert curve
    outfile << "{\n" << "appearance {linewidth 2}\n\n";
    outfile <<  "{ < " << geomdir << "hilbert_curve.vect }\n";
    outfile << "\n}\n\n";

    //points of hilbert curve
    outfile << "{\n" << "appearance {linewidth 5}\n\n";
    outfile <<  "{ < " << geomdir << "hilbert_curve_points.vect }\n";
    outfile << "\n}\n";
    outfile.close();


    //OOGL file display the bounding boxes + curve + points
    outfile.open(outdir + "hcurve_bbox.list");
    outfile << "LIST\n\n";

    //hilbert curve
    outfile << "{\n" << "appearance {linewidth 2}\n\n";
    outfile <<  "{ < " << geomdir << "hilbert_curve.vect }\n";
    outfile << "\n}\n\n";

    //points of hilbert curve
    outfile << "{\n" << "appearance {linewidth 5}\n\n";
    outfile <<  "{ < " << geomdir << "hilbert_curve_points.vect }\n";
    outfile << "\n}\n\n";
    
    //bounding boxes
    outfile << "{\n" << "appearance {+edge}\n\n";
    outfile << "LIST\n";
    for( int i = 0; i < num_hse; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile << "{ < " << geomdir + geomfile << " }\n";
    }
    outfile << "\n}\n";
    outfile.close();

}



