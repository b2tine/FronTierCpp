#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/hilbert_sort.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <limits>
#include <map>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point3;
typedef CGAL::Surface_mesh<Point3> Mesh;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Face_index Face_index;


double minval(double a, double b, double c)
{
    double val = std::min(a,b);
    return std::min(val,c);
}

double maxval(double a, double b, double c)
{
    double val = std::max(a,b);
    return std::max(val,c);
}

using tri_pts = std::vector<Point3>;

struct AABB
{
    std::vector<double> min;
    std::vector<double> max;
    std::vector<double> ctr;

    AABB() = default;
    AABB(const AABB&) = default;
    AABB& operator=(const AABB&) = default;
    ~AABB() = default;

    explicit AABB(const tri_pts& pts)
    {
        min.push_back(minval(pts[0].x(),
                pts[1].x(), pts[2].x()));
        min.push_back(minval(pts[0].y(),
                pts[1].y(), pts[2].y()));
        min.push_back(minval(pts[0].z(),
                pts[1].z(), pts[2].z()));
        max.push_back(maxval(pts[0].x(),
                pts[1].x(), pts[2].x()));
        max.push_back(maxval(pts[0].y(),
                pts[1].y(), pts[2].y()));
        max.push_back(maxval(pts[0].z(),
                pts[1].z(), pts[2].z()));
        ctr.push_back(0.5*(min[0] + max[0]));
        ctr.push_back(0.5*(min[1] + max[1]));
        ctr.push_back(0.5*(min[2] + max[2]));
    }

    Point3 center()
    {
        return Point3(ctr[0],ctr[1],ctr[2]);
    }

};


extern void create_directory(std::string);


int main(int argc, char** argv)
{
    int num_names = 0; // number of filenames
    char* filename[2];
    bool help = false;
    for (int i = 1; i < argc; i++)
    {
	    if( strcmp("-h", argv[i]) == 0 ||
                strcmp("--h", argv[i]) == 0 ||
                strcmp("-help", argv[i]) == 0 ||
                strcmp("--help", argv[i]) == 0 )
        {
            help = true;
        }
        else if( num_names < 2 )
        {
            filename[num_names++] = argv[i];
        }
        else
        {
            ++num_names;
            break;
        }
    }

    if( num_names != 2 || help )
    {
        if( !help )
            std::cerr << "Error: in parameter list\n";
        std::cerr << "Usage: " << argv[0] << " [infile] [outdir]\n";
        std::cerr << "where [infile] is an OFF file.\n";
        exit( !help );
    }

    const char* infile = filename[0];
    std::ifstream input(infile);

    Mesh mesh;
    if ( !input || !(input >> mesh) || mesh.is_empty()
            || !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Not a valid input file." << "\n";
        exit(1);
    }

    std::string outdir(argv[2]);
    outdir += "/";
    std::string geomdir("OOGL/");
    create_directory(outdir + geomdir);

    std::ofstream outfile(outdir + geomdir + "input-mesh.off");
    outfile << mesh;
    outfile.close();


    int face_id = 0;
    const unsigned int num_faces = mesh.number_of_faces();
    std::vector<AABB> boxvec;
    std::vector<Point3> ctrvec;
    std::map<Point3,AABB> boxmap;

    //compute bbox for each triangle of the mesh
    Mesh::Face_range frange = mesh.faces();
    Mesh::Face_range::iterator fit = frange.begin();
    while( fit != frange.end() )
    { 
        CGAL::Vertex_around_face_circulator<Mesh> vb(mesh.halfedge(*fit),mesh); 
        CGAL::Vertex_around_face_circulator<Mesh> ve(vb);
        tri_pts tps;

        do {
            Point3 p = mesh.point(*vb);
            tps.push_back(p);
            vb++;
        } while( vb != ve );

        AABB bv(tps);
        ctrvec.push_back(bv.center());
        boxmap[ctrvec[face_id++]] = bv;
        boxvec.push_back(bv);

        fit++;
    }

    //For Debugging purposes
    outfile.open(outdir + "unsorted_bbox_centers.txt");
    for( int i = 0; i < num_faces; i++ )
    {
        outfile << i << "\t";
        for(int k = 0; k < 3; k++ )
        {
            outfile << boxvec[i].ctr[k] << "\t";
        }
        outfile << "\n";
    }
    outfile.close();



    //Sort bounding box centers on a hilbert curve
    CGAL::hilbert_sort(ctrvec.begin(),ctrvec.end());


    std::ofstream outfile1(outdir + "sorted_bbox_centers.txt");

    //Write geomview OOGL file for the curve
    outfile.open(outdir + geomdir + "hilbert_curve.vect");
    outfile << "VECT\n\n";
    //num_polylines total_num_vertices num_colors
    outfile << 1 << "\t" << num_faces << "\t" << 1 << "\n\n";
    //number vertices in each polyline (would be a list if more than one polyline)
    outfile << num_faces << "\n\n";
    //number of colors in each polyline (would be a list if more than one polyline)
    outfile << 1 << "\n";
    //coordinates of curve vertices (same as the points/boxcenters)
    for( int i = 0; i < num_faces; i++ )
    {
        std::vector<double> coords = boxmap[ctrvec[i]].ctr;
        for( int k = 0; k < 3; k++ )
        {
            outfile << coords[k] << "\t";
            outfile1 << coords[k] << "\t";
        }
        outfile << "\n";
        outfile1 << "\n";
    }
    outfile << "\n";
    //red, green, blue, alpha(opacity)
    outfile << 1.0 << "\t" << 0.0 << "\t" <<  0.0 << "\t" <<  1.0;//curve color (red)
    outfile.close();
    outfile1.close();


    //Write geomview OOGL file for the points of the hilbert
    //curve so we can display in a different color
    outfile.open(outdir + geomdir + "hilbert_curve_points.vect");
    outfile << "VECT\n\n";
    //num_polylines total_num_vertices num_colors
    outfile << num_faces << "\t" << num_faces << "\t" << 1 << "\n\n";
    //number vertices in each polyline (degenerate polyline is a point)
    for( int i = 0; i < num_faces; i++ )
        outfile << 1 << "\n";
    outfile << "\n";
    //first point is assigned the color (indexed by 1) given below
    outfile << 1 << "\n";
    //remaining points get same color as first point (indicated by 0)
    for( int i = 1; i < num_faces; i++ )
        outfile << 0 << "\n";
    outfile << "\n";
    //coordinates of the curve points
    for( int i = 0; i < num_faces; i++ )
    {
        std::vector<double> coords = boxmap[ctrvec[i]].ctr;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << "\t";
        outfile << "\n";
    }
    outfile << "\n";
    //red, green, blue, alpha(opacity)
    outfile << 0.0 << "\t" << 1.0 << "\t" <<  0.0 << "\t" <<  1.0;//point color (green)
    outfile.close();


    //write sorted bboxes to geomview OOOGL files
    char geomfile[100];
    for( int i = 0; i < num_faces; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile.open(outdir + geomdir + geomfile);
        outfile << "BBOX\n";
        std::vector<double> coords = boxmap[ctrvec[i]].min;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << "\t";
        outfile << "\n";
        coords = boxmap[ctrvec[i]].max;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << "\t";
        outfile.close();
    }



    //OOGL file to display mesh + curve + points
    outfile.open(outdir + "hcurve.list");
    outfile << "LIST\n\n";

    //input mesh
    outfile << "{\n" << "appearance {\n\t+edge +face +transparent\n\t";
    outfile << "material {\n\t\tedgecolor 0 0 1\n\t\t";
    outfile << "ambient 0 0 0.5\n\t\talpha 0.45\n \t\t}\n\n}";
    outfile << "\n\n{ < " << geomdir << "input-mesh.off }\n\n}\n\n";

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

    /*
    //input mesh
    outfile << "{\n" << "appearance {\n\t+edge +face +transparent\n\t";
    outfile << "material {\n\t\tedgecolor 0 0 1\n\t\t";
    outfile << "ambient 0 0 0.5\n\t\talpha 0.45\n \t\t}\n\n}";
    outfile << "\n\n{ < " << geomdir << "input-mesh.off }\n\n}\n\n";
    */

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
    for( int i = 0; i < num_faces; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile << "{ < " << geomdir + geomfile << " }\n";
    }
    outfile << "\n}\n";
    outfile.close();





    return 0;
}





