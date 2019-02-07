#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/hilbert_sort.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <limits>
#include <map>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 CGAL_Point3;
typedef CGAL::Surface_mesh<CGAL_Point3> Mesh;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Face_index Face_index;

//struct MyPoint3
class MyPoint3
{
    private:
    double x{0}, y{0}, z{0};

    public:
    MyPoint3() = default;
    MyPoint3(const MyPoint3&) = default;
    MyPoint3& operator=(const MyPoint3&) = default;
    MyPoint3(MyPoint3&&) = default;
    MyPoint3& operator=(MyPoint3&&) = default;
    ~MyPoint3() = default;

    MyPoint3(double X, double Y, double Z)
        : x{X}, y{Y}, z{Z}
    {}

    const double& operator[](const std::size_t i) const
    {
        assert( i >= 0 && i <= 2 );
        switch(i)
        {
            case 0: return x; break;
            case 1: return y; break;
            case 2: return z; break;
        }
    }

    double& operator[](const std::size_t i)
    {
        //call const version and cast off the const,
        //but need to cast *this as const first
        return const_cast<double&>(
                static_cast<const MyPoint3&>(*this)[i]);
    }

    //needed for use as key in std::map<key,val>
    bool operator < (const MyPoint3& p) const
    {
        if( x == p[0] )
        {
            if( y == p[1] )
            {
                return z < p[2];
            }
            return y < p[1];
        }
        return x < p[0];
    }

    //not necessary for std::map, but doesn't hurt
    bool operator == (const MyPoint3& p) const
    {
        return x == p[0] && y == p[1] && z == p[2];
    }
};

//these needed for use in hilbert_sort()
struct MyLessX
{
    bool operator()(const MyPoint3& p, const MyPoint3& q) const
    {
        return p[0] < q[0];
        //return p.x < q.x;
    }
};

struct MyLessY
{
    bool operator()(const MyPoint3& p, const MyPoint3& q) const
    {
        return p[1] < q[1];
        //return p.y < q.y;
    }
};

struct MyLessZ
{
    bool operator()(const MyPoint3& p, const MyPoint3& q) const
    {
        return p[2] < q[2];
        //return p.z < q.z;
    }
};

struct MyHilbertSortingTraits
{
    using Point_3 = MyPoint3;
    using Less_x_3 = MyLessX;
    using Less_y_3 = MyLessY;
    using Less_z_3 = MyLessZ;

    Less_x_3 less_x_3_object() const
    {
        return Less_x_3();
    }

    Less_y_3 less_y_3_object() const
    {
        return Less_y_3();
    }

    Less_z_3 less_z_3_object() const
    {
        return Less_z_3();
    }
};


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

using tri_pts = std::vector<CGAL_Point3>;

struct AABB
{
    MyPoint3 min;
    MyPoint3 max;
    MyPoint3 ctr;

    AABB() = default;
    AABB(const AABB&) = default;
    AABB& operator=(const AABB&) = default;
    ~AABB() = default;

    explicit AABB(const tri_pts& pts)
    {
        min[0] = minval(pts[0].x(),pts[1].x(), pts[2].x());
        min[1] = minval(pts[0].y(),pts[1].y(), pts[2].y());
        min[2] = minval(pts[0].z(),pts[1].z(), pts[2].z());

        max[0] = maxval(pts[0].x(),pts[1].x(), pts[2].x());
        max[1] = maxval(pts[0].y(),pts[1].y(), pts[2].y());
        max[2] = maxval(pts[0].z(),pts[1].z(), pts[2].z());

        ctr[0] = 0.5*(min[0] + max[0]);
        ctr[1] = 0.5*(min[1] + max[1]);
        ctr[2]  = 0.5*(min[2] + max[2]);
        /*
        min.x = minval(pts[0].x(),pts[1].x(), pts[2].x());
        min.y = minval(pts[0].y(),pts[1].y(), pts[2].y());
        min.z = minval(pts[0].z(),pts[1].z(), pts[2].z());

        max.x = maxval(pts[0].x(),pts[1].x(), pts[2].x());
        max.y = maxval(pts[0].y(),pts[1].y(), pts[2].y());
        max.z = maxval(pts[0].z(),pts[1].z(), pts[2].z());

        ctr.x = 0.5*(min.x + max.x);
        ctr.y = 0.5*(min.y + max.y);
        ctr.z  = 0.5*(min.z + max.z);*/
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
    std::vector<MyPoint3> ctrvec;
    std::map<MyPoint3,AABB> boxmap;

    //compute bbox for each triangle of the mesh
    Mesh::Face_range frange = mesh.faces();
    Mesh::Face_range::iterator fit = frange.begin();
    while( fit != frange.end() )
    { 
        CGAL::Vertex_around_face_circulator<Mesh> vb(mesh.halfedge(*fit),mesh); 
        CGAL::Vertex_around_face_circulator<Mesh> ve(vb);
        tri_pts tps;

        do {
            CGAL_Point3 p = mesh.point(*vb);
            tps.push_back(p);
            vb++;
        } while( vb != ve );

        AABB bv(tps);
        ctrvec.push_back(bv.ctr);
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



    MyHilbertSortingTraits hst;
    //Sort bounding box centers on a hilbert curve
    CGAL::hilbert_sort(ctrvec.begin(),ctrvec.end(),hst);



    //Write geomview OOGL file for the curve
    outfile.open(outdir + geomdir + "hilbert_curve.vect");
    outfile << "VECT\n\n";
    //num_polylines total_num_vertices num_colors
    outfile << 1 << " " << num_faces << " " << 1 << "\n\n";
    //number vertices in each polyline (would be a list if more than one polyline)
    outfile << num_faces << "\n\n";
    //number of colors in each polyline (would be a list if more than one polyline)
    outfile << 1 << "\n";
    //coordinates of curve vertices (same as the points/boxcenters)
    for( int i = 0; i < num_faces; i++ )
    {
        MyPoint3 coords = boxmap[ctrvec[i]].ctr;
        //std::vector<double> coords = boxmap[ctrvec[i]].ctr;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile << "\n";
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
    outfile << num_faces << " " << num_faces << " " << 1 << "\n\n";
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
        MyPoint3 coords = boxmap[ctrvec[i]].ctr;
        //std::vector<double> coords = boxmap[ctrvec[i]].ctr;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile << "\n";
    }
    outfile << "\n";
    //red, green, blue, alpha(opacity)
    outfile << 0.0 << " " << 1.0 << " " <<  0.0 << " " <<  1.0;//point color (green)
    outfile.close();


    //write sorted bboxes to geomview OOOGL files
    char geomfile[100];
    for( int i = 0; i < num_faces; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile.open(outdir + geomdir + geomfile);
        outfile << "BBOX\n";
        MyPoint3 coords = boxmap[ctrvec[i]].min;
        //std::vector<double> coords = boxmap[ctrvec[i]].min;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile << "\n";
        coords = boxmap[ctrvec[i]].max;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
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


    /*
    geomdir = "HSEGS/";
    create_directory(outdir + geomdir);
    //write hilbert curve segments OOOGL files for animation
    for( int i = 0; i < num_faces-1; i++ )
    {
        sprintf(geomfile, "hseg%05d.vect", i);
        outfile.open(outdir + geomdir + geomfile);
        outfile << "VECT\n\n";
        //num_polylines total_num_vertices num_colors
        outfile << 1 << " " << 2 << " " << 1 << "\n\n";
        //number vertices in each polyline (segment)
        outfile << 2 << "\n\n";
        //number of colors in each polyline (would be a list if more than one polyline)
        outfile << 1 << "\n";
        //coordinates of segment end points
        for( int j = 0; j < 2; j++ )
        {
            std::vector<double> coords = boxmap[ctrvec[i+j]].ctr;
            for( int k = 0; k < 3; k++ )
                outfile << coords[k] << " ";
            outfile << "\n";
        }
        outfile << "\n";
        //red, green, blue, alpha(opacity)
        outfile << 1.0 << " " << 0.0 << " " <<  0.0 << " " <<  1.0;//curve color (red)
        outfile.close();
    }
    */


    return 0;
}





