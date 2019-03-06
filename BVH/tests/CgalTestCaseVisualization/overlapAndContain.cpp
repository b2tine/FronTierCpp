#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

int main()
{
    Mesh mesh;

    vertex_descriptor b = mesh.add_vertex(Point(1,1,1));
    vertex_descriptor d = mesh.add_vertex(Point(10,0,0));
    vertex_descriptor e = mesh.add_vertex(Point(0,10,0));
    vertex_descriptor f = mesh.add_vertex(Point(0,0,10));

    mesh.add_face(d,e,f);
    mesh.add_face(e,d,b);

    Mesh::Vertex_range vrange = mesh.vertices();
    Mesh::Vertex_range::iterator vit = vrange.begin();

    while( vit != vrange.end() )
    {
        Point p = mesh.point(*vit);

        std::cout << p.x() << "\t" << p.y() << "\t" << p.z() << "\n";
        vit++;
    }

    std::cout << "\n\n\n";

    Mesh::Face_range frange = mesh.faces();
    Mesh::Face_range::iterator fit = frange.begin();

    while( fit != frange.end() )
    { 
        CGAL::Vertex_around_face_circulator<Mesh> vb(mesh.halfedge(*fit),mesh); 
        CGAL::Vertex_around_face_circulator<Mesh> ve(vb);

        std::cout << "Facet" << "\n\n";

        do {
            Point p = mesh.point(*vb);
            std::cout << p.x() << "\t" << p.y() << "\t" << p.z() << "\n";
            vb++;
        } while( vb != ve );

        std::cout << "\n";
        fit++;
    }


    std::ofstream outfile("tris.off");
    outfile << mesh;
    outfile.close();

    return 0;
}


