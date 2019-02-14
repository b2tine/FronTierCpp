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

    vertex_descriptor a = mesh.add_vertex(Point(0,0,0));
    vertex_descriptor b = mesh.add_vertex(Point(1,0,0));
    vertex_descriptor c = mesh.add_vertex(Point(0,1,0));
    vertex_descriptor d = mesh.add_vertex(Point(1,1,0));
    vertex_descriptor e = mesh.add_vertex(Point(0.75,0,0.75));
    //vertex_descriptor e = mesh.add_vertex(Point(1.5,0.5,1));
    vertex_descriptor f = mesh.add_vertex(Point(2,0,1));
    vertex_descriptor g = mesh.add_vertex(Point(0.75,0.75,1));
    //vertex_descriptor g = mesh.add_vertex(Point(-0.5,0.5,1.5));
    vertex_descriptor h = mesh.add_vertex(Point(0.75,-0.5,0.25));

    mesh.add_face(a,b,c);
    mesh.add_face(b,d,c);
    mesh.add_face(b,e,d);
    mesh.add_face(b,f,e);
    mesh.add_face(a,c,g);
    mesh.add_face(a,g,h);

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


