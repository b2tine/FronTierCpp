#include "BVH.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>


static void mkdirTree(std::string sub, std::string dir);
static void createDirectory(std::string new_dir);


/*
void BVH::drawHeirarchy() const
{
    //TODO: Need to write a SpaceFillingCurve class.
    //      Can't just postorder traverse because the
    //      order depends on the type of curve.
    //      Want to be able to use different curves on
    //      each level of the heirarchy also.
}
*/


void BVH::DrawUnlock() noexcept
{
    isDrawTime = true;
}

void BVH::setDrawDirectory(std::string dir) noexcept
{
    outdir = dir;
}


void BVH::drawHeirarchyLevel() const
{
    printf("\nnumnodes = %d\n",children.size());
    for (int i = 0; i < children.size(); ++i)
    {
        auto node = children[i].second;
        printf("node->level = %d\n",node->level);
    }
    printf("\n\n");

    int dummytstep = 0;
    if( isDrawTime == true )
        writeHilbertCurveFiles(dummytstep);
}


//TODO: Break out into several functions for reuse/flexibility
void BVH::writeHilbertCurveFiles(int step) const
{
    assert(!drawdir.empty());
    int level = this->sort_iter;

    std::string lvlid = std::string(2,'0').append(std::to_string(level));
    std::string lvldir = drawdir + "/level-" + lvlid + "/";

    std::string geomdir("OOGL/");
    createDirectory(lvldir + geomdir);

    int numBV = children.size();

    //Write geomview OOGL file for the curve
    std::ofstream outfile(lvldir + geomdir + "hilbert_curve.vect");
    
    outfile << "VECT\n\n";
    
    //num_polylines total_num_vertices num_colors
    outfile << 1 << " " << numBV << " " << 1 << "\n\n";

    //number vertices in each polyline (would be a list if more than one polyline)
    outfile << numBV << "\n\n";
    
    //number of colors in each polyline (would be a list if more than one polyline)
    outfile << 1 << "\n";
    
    //coordinates of curve vertices (same as the points/boxcenters)
    for( int i = 0; i < numBV; i++ )
    {
        CGAL_Point p = children[i].first;
        outfile << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    outfile << "\n";
             //red,          green,         blue,     alpha(opacity)
    outfile << 1.0 << " " << 0.0 << " " <<  0.0 << " " <<  1.0;      //curve color (red)
    outfile.close();


    //Write geomview OOGL file for the points of the hilbert
    //curve so we can display in a different color
    outfile.open(lvldir + geomdir + "hilbert_curve_points.vect");
    
    outfile << "VECT\n\n";
    
    //num_polylines total_num_vertices num_colors
    outfile << numBV << " " << numBV << " " << 1 << "\n\n";
    
    //number vertices in each polyline (degenerate polyline is a point)
    for( int i = 0; i < numBV; i++ )
        outfile << 1 << "\n";
    outfile << "\n";
    
    //first point is assigned the color (indexed by 1) given below
    outfile << 1 << "\n";
    
    //remaining points get same color as first point (indicated by 0)
    for( int i = 1; i < numBV; i++ )
        outfile << 0 << "\n";
    outfile << "\n";
    
    //coordinates of the curve points
    for( int i = 0; i < numBV; i++ )
    {
        CGAL_Point p = children[i].first;
        outfile << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    outfile << "\n"; 
              //red,        green,         blue,        alpha(opacity)
    outfile << 0.0 << " " << 1.0 << " " <<  0.0 << " " <<  1.0;        //point color (green)
    outfile.close();


    //write sorted bboxes to geomview OOGL files
    char geomfile[100];
    for( int i = 0; i < numBV; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile.open(lvldir + geomdir + geomfile);
        outfile << "BBOX\n";

        BV_Point lcoords = children[i].second->getBV().lower;
        for( int k = 0; k < 3; k++ )
            outfile << lcoords[k] << " ";
        outfile << "\n";

        BV_Point ucoords = children[i].second->getBV().upper;
        for( int k = 0; k < 3; k++ )
            outfile << ucoords[k] << " ";
        outfile.close();
    }


    //TODO: Generalize this next output routine to the case
    //      where there is no input mesh provided.

    //OOGL file to display mesh + curve + points
    outfile.open(lvldir + "hcurve.list");
    outfile << "LIST\n\n";

    //input mesh
    outfile << "{\n" << "appearance {\n\t+edge +face +transparent\n\t";
    outfile << "material {\n\t\tedgecolor 0 0 1\n\t\t";
    outfile << "ambient 0 0 0.5\n\t\talpha 0.45\n \t\t}\n\n}";
    outfile << "\n\n{ < ../input-mesh.off }\n\n}\n\n";

    //hilbert curve
    outfile << "{\n" << "appearance {linewidth 2}\n\n";
    outfile <<  "{ < " << geomdir << "hilbert_curve.vect }\n";
    outfile << "\n}\n\n";

    //points of hilbert curve
    outfile << "{\n" << "appearance {linewidth 5}\n\n";
    outfile <<  "{ < " << geomdir << "hilbert_curve_points.vect }\n";
    outfile << "\n}\n";
    outfile.close();

    //END TODO


    //OOGL file display the bounding boxes + curve + points
    outfile.open(lvldir + "hcurve_bbox.list");
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
    for( int i = 0; i < numBV; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile << "{ < " << geomdir + geomfile << " }\n";
    }
    outfile << "\n}\n";
    outfile.close();

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


void createDirectory(std::string new_dir)
{
    struct stat st;
    int status = stat(new_dir.c_str(), &st);
    if( status != 0 && !S_ISDIR(st.st_mode) )
        mkdirTree(new_dir, "");
}

