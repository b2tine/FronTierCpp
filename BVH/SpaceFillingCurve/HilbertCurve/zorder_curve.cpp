#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <limits>
#include <map>


typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Face_index Face_index;


float minval(float a, float b, float c)
{
    float val = std::min(a,b);
    return std::min(val,c);
}

float maxval(float a, float b, float c)
{
    float val = std::max(a,b);
    return std::max(val,c);
}

using tri_pts = std::vector<Point>;

struct AABB
{
    std::vector<float> min;
    std::vector<float> max;
    std::vector<float> ctr;

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

    void translate(const std::vector<float>& vec)
    {
        for( int i = 0; i < 3; i++ )
        {
            min[i] += vec[i];
            max[i] += vec[i];
            ctr[i] += vec[i];
        }
    }
};


extern uint32_t morton3d(float, float, float);
extern void radix_sort(uint32_t*, const unsigned);
extern void check_order(uint32_t*, const unsigned);
extern int check_duplicates(uint32_t*, const unsigned);
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
    //create_directory(outdir);

    std::string geomdir("OOGL/");
    create_directory(outdir + geomdir);

    std::ofstream outfile(outdir + geomdir + "input-mesh.off");
    outfile << mesh;
    outfile.close();


    int face_id = 0;
    const unsigned int num_faces = mesh.number_of_faces();
    //std::map<Face_index,int> fmap;
    std::vector<float> bbox_ctr_min(3,std::numeric_limits<float>::max());
    std::vector<AABB> boxvec;

    //compute bbox for each triangle of the mesh
    Mesh::Face_range frange = mesh.faces();
    Mesh::Face_range::iterator fit = frange.begin();
    while( fit != frange.end() )
    { 
        CGAL::Vertex_around_face_circulator<Mesh> vb(mesh.halfedge(*fit),mesh); 
        CGAL::Vertex_around_face_circulator<Mesh> ve(vb);
        tri_pts tps;

        do {
            Point p = mesh.point(*vb);
            tps.push_back(p);
            vb++;
        } while( vb != ve );

        AABB bv(tps);
        boxvec.push_back(bv);
        for( int i = 0; i < 3; i++ )
        {
            if( boxvec[face_id].ctr[i] < bbox_ctr_min[i] )
                bbox_ctr_min[i] = boxvec[face_id].ctr[i];
        }
   
        //fmap[*fit] = face_id++;        
        face_id++;
        fit++;
    }

    outfile.open(outdir + "boxcenters.txt");
    for( int i = 0; i < num_faces; i++ )
    {
        outfile << i << "\t";
        for(int k = 0; k < 3; k++ )
        {
            outfile << boxvec[i].ctr[k] << "\t";
        }
        outfile << "\n";
    }

    outfile << "\n\n";
    outfile << "min coord: \t";
    for( int i = 0; i < 3; i++ )
        outfile << bbox_ctr_min[i] << "\t";
    outfile << "\n\n";

    //calculate displacement of translation
    std::vector<float> shift(3,0.0);
    outfile << "translation vector: \t";
    for(int i = 0; i < 3; i++ )
    {
        if( bbox_ctr_min[i] < 0.0 )
            shift[i] += fabs(bbox_ctr_min[i]);//may need a rounding tolerance pad
        outfile << shift[i] << "\t";
    }
    outfile.close();

    //translate coords of bbox to 1st octant to ensure uniqueness of morton codes
    std::vector<AABB> tbvec = boxvec;
    outfile.open(outdir + "translated_centers.txt");
    for( int i = 0; i < num_faces; i++ )
    {
        outfile << i << "\t";
        tbvec[i].translate(shift);
        for( int k = 0; k < 3; k++ )
            outfile << tbvec[i].ctr[k] << "\t";
        outfile << "\n";
    }
    outfile.close();

    /*
    float** ctrs = new float*[num_faces];
    for( int i = 0; i < num_faces; i++ )
    {
        ctrs[i] = new float[3];
        for( int k = 0; k < 3; k++ )
            ctrs[i][k] = boxvec[i].ctr[k] + shift[k];
    }

    unsigned int* mcode = new unsigned int[num_faces];
    for( int i = 0; i < num_faces; i++ )
        mcode[i] = morton3d(ctrs[i][0], ctrs[i][1], ctrs[i][2]);
    */

    //compute morton code for each translated bbox center
    std::map<unsigned int,AABB> codemap;
    std::map<unsigned int,AABB> t_codemap;
    unsigned int* mcode = new unsigned int[num_faces];
    for( int i = 0; i < num_faces; i++ )
    {
        for( int k = 0; k < 3; k++ )
            assert(tbvec[i].ctr[k] >= 0.0);

        mcode[i] = morton3d(tbvec[i].ctr[0],
                tbvec[i].ctr[1], tbvec[i].ctr[2]);
        codemap[mcode[i]] = boxvec[i];
        t_codemap[mcode[i]] = tbvec[i];
    }

    radix_sort(mcode, num_faces);   //sort morton codes
    check_order(mcode, num_faces);  //check success of the radix sort

    outfile.open(outdir + "sorted_mcodes.txt");
    for( int i = 0; i < num_faces; i++ )
        outfile << i << "\t" << mcode[i] << "\n";

    //check for duplicate morton codes
    int duplicate_index = check_duplicates(mcode, num_faces);
    if( duplicate_index > 0 )
    {
        std::cout << "duplicate morton code at index " << duplicate_index << "\n";
        outfile << "\n" << "duplicate morton code at index " << duplicate_index << "\n";
    }
    outfile.close();


    //write the z-order curve graphics file from sorted morton codes
    outfile.open(outdir + geomdir + "z-order_curve.vect");
    outfile << "VECT\n";
    //num_polylines total_num_vertices num_colors
    outfile << 1 << " " << num_faces << " " << 1 << "\n";
    //number vertices in each polyline (would be a list if more than one polyline)
    outfile << num_faces << "\n";
    //number of colors in each polyline (would be a list if more than one polyline)
    outfile << 1 << "\n";
    for( int i = 0; i < num_faces; i++ )
    {
        std::vector<float> coords = codemap[mcode[i]].ctr;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile << "\n";
    }
    //red, green, blue, alpha(opacity)
    outfile << 1.0 << " " << 0.0 << " " <<  0.0 << " " <<  1.0;
    outfile.close();

    //write sorted bboxes graphics files
    char geomfile[100];
    for( int i = 0; i < num_faces; i++ )
    {
        sprintf(geomfile, "aabb%05d.bbox", i);
        outfile.open(outdir + geomdir + geomfile);
        outfile << "BBOX\n";
        std::vector<float> coords = codemap[mcode[i]].min;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile << "\n";
        coords = codemap[mcode[i]].max;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile.close();
    }

    //write file to load all graphics file
    outfile.open(outdir + "zcurve_bbox.list");
    outfile << "LIST\n\n";

    //input mesh
    outfile << "{\n" << "appearance {\n\t+edge +face +transparent\n\t";
    outfile << "material {\n\t\tedgecolor 0 0 1\n\t\t";
    outfile << "ambient 0 0 0.5\n\t\talpha 0.45\n \t\t}\n\n}";
    outfile << "\n\n{ < " << geomdir << "input-mesh.off }\n\n}\n\n";

    //z-order curve
    outfile << "{\n" << "appearance {linewidth 2}\n\n";
    outfile <<  "{ < " << geomdir << "z-order_curve.vect }\n";
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

    //FOR DEBUGGING
    //write the z-order curve graphics file for translated boxes
    outfile.open(outdir + geomdir + "z-order_curve-t.vect");
    outfile << "VECT\n";
    outfile << 1 << " " << num_faces << " " << 1 << "\n";
    outfile << num_faces << "\n";
    outfile << 1 << "\n";
    for( int i = 0; i < num_faces; i++ )
    {
        std::vector<float> coords = t_codemap[mcode[i]].ctr;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile << "\n";
    }
    //red, green, blue, alpha(opacity)
    outfile << 1.0 << " " << 0.0 << " " <<  0.0 << " " <<  1.0;
    outfile.close();

    //write sorted translated boxes graphics files
    for( int i = 0; i < num_faces; i++ )
    {
        sprintf(geomfile, "aabb-t%05d.bbox", i);
        outfile.open(outdir + geomdir + geomfile);
        outfile << "BBOX\n";
        std::vector<float> coords = t_codemap[mcode[i]].min;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile << "\n";
        coords = t_codemap[mcode[i]].max;
        for( int k = 0; k < 3; k++ )
            outfile << coords[k] << " ";
        outfile.close();
    }

    //write file to load all graphics file
    outfile.open(outdir + "zcurve_bbox-t.list");
    outfile << "LIST\n\n";

    /*
    //TODO: need to translate mesh as well
    //input mesh
    outfile << "{\n" << "appearance {\n\t+edge +face +transparent\n\t";
    outfile << "material {\n\t\tedgecolor 0 0 1\n\t\t";
    outfile << "ambient 0 0 0.5\n\t\talpha 0.45\n \t\t}\n\n}";
    outfile << "\n\n{ < " << geomdir << "input-mesh.off }\n\n}\n\n";
    */

    //z-order curve
    outfile << "{\n" << "appearance {linewidth 2}\n\n";
    outfile <<  "{ < " << geomdir << "z-order_curve-t.vect }\n";
    outfile << "\n}\n\n";
    
    //bounding boxes
    outfile << "{\n" << "appearance {+edge}\n\n";
    outfile << "LIST\n";
    for( int i = 0; i < num_faces; i++ )
    {
        sprintf(geomfile, "aabb-t%05d.bbox", i);
        outfile << "{ < " << geomdir + geomfile << " }\n";
    }
    outfile << "\n}\n";
    outfile.close();

    delete[] mcode;
    return 0;
}


