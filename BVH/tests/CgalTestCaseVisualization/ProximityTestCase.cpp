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
    vertex_descriptor c = mesh.add_vertex(Point(0,1,1));
    vertex_descriptor u = mesh.add_vertex(Point(-1,1,0));
    vertex_descriptor v = mesh.add_vertex(Point(-0.5,2,1));
    vertex_descriptor w = mesh.add_vertex(Point(-0.25,2.5,0.5));
    //vertex_descriptor d = mesh.add_vertex(Point(-1,1,0));
    //vertex_descriptor e = mesh.add_vertex(Point(-0.5,2,1));

    mesh.add_face(a,b,c);
    mesh.add_face(a,c,u);
    //mesh.add_face(a,c,d);
    mesh.add_face(u,c,v);
    mesh.add_face(v,c,w);
    

    vertex_descriptor d = mesh.add_vertex(Point(1,1,0));    //g
    vertex_descriptor e = mesh.add_vertex(Point(0,1,0.75)); //h
    vertex_descriptor f = mesh.add_vertex(Point(2,0,1));    //i
    vertex_descriptor g = mesh.add_vertex(Point(1.5,1.5,1));//j
    vertex_descriptor h = mesh.add_vertex(Point(2.5,1,0));  //k
    vertex_descriptor i = mesh.add_vertex(Point(2.75,2,0.25));//l

    mesh.add_face(e,d,f);//h,g,i
    mesh.add_face(f,g,e);//i,j,h
    mesh.add_face(f,h,g);//i,k,j
    mesh.add_face(h,i,g);//k,l,j


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


