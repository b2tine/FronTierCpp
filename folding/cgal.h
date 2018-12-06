#ifndef CGAL_H_
#define CGAL_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <iostream>
#include <fstream>
#include <FronTier.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,Itag> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Cgal_Point;

bool findAndLocate(std::ifstream&, const char*);

class cgalSurf
{
    INTERFACE* c_intfc;
    SURFACE** c_surf; 
    CDT cdt; 
    int num_finite_face; 
    double k_dx; 
    // used to provide upper bound for minimum angle 
    double c_bound; 
    double height_; 
    std::vector<double> extConPoint_;
    bool hole_ = false; 
    // non-virtual interface idiom (Template Method design patter)
    // Only make different part virtual and keep 
    // other common part directly inherited by derived class
    virtual void addHole(std::vector<bool>&, int&);
    virtual void getSpecialParaFromFile(std::ifstream&) = 0;
protected :
    void height(double h) { height_ = h; } 
    double height() const { return height_; }
    void CGALCoeffRestricSize(double k) { k_dx = k; }
    double CGALCoeffRestricSize() const { return k_dx; }
    void CGALMinAngleUb(double c) { c_bound = c; }
    double CGALMinAngleUb() const { return c_bound; }
    Vertex_handle insertPointToCDT(Cgal_Point a) { return cdt.insert(a); }
    void insertConstraintToCDT(Vertex_handle v1, Vertex_handle v2) { 
	cdt.insert_constraint(v1, v2); }
    CDT readCDT() const { return cdt; }
    CDT& readCDT() { return cdt; }
    std::vector<double> extConPoint() const { return extConPoint_; }
    void extConPoint(double p, int index) { extConPoint_[index] = p; }
    INTERFACE* interface() { return c_intfc; }
    int numFinFace() const { return num_finite_face; }
    void numFinFace(int num) { num_finite_face = num; }
    SURFACE** surface() { return c_surf; }
    bool hole() const { return hole_; }
    void hole(bool b) { hole_ = b; }
public :
    // std::ifstream& can't be used in constructor
    // any constructor involving it been added =delete in the end
    cgalSurf(INTERFACE*, SURFACE**);   
    // add extra constraint
    void getExtraConstPoint(std::ifstream&); 
    void setSurfZeroMesh(); 
    void setCurveZeroLength(CURVE*, double); 
    void setMonoCompBdryZeroLength(); 
    // by generated tri mesh generates surface
    void cgalGenSurf(); 
    // add boundary constraint. different for differernt shape
    virtual void addCgalConst() = 0;
    // preprocess and generate tri mesh 
    void getParaFromFile(std::ifstream&); 
    void cgalTriMesh(std::ifstream&); 
    virtual ~cgalSurf() {}; 
};

class cgalRectangleSurf : public cgalSurf
{
    double lower[2]; 
    double upper[2]; 
    void getSpecialParaFromFile(std::ifstream&);
public :
    cgalRectangleSurf(INTERFACE* intfc, SURFACE** surf) : 
    		cgalSurf(intfc, surf) {}
    void addCgalConst(); 
    ~cgalRectangleSurf() {}
};

class cgalCircleSurf : public cgalSurf
{
    double cen[2]; 
    double radius; 
    int num_reg_const;
    static const double eps;
    virtual void getSpecialParaFromFile(std::ifstream&);
protected : 
    double* getCenter() { return cen; }
    void setCenter(double c1, double c2) { cen[0] = c1; cen[1] = c2; }
    double getRadius() const { return radius; }
    void setRadius(double r) { radius = r; }
    int numRegConst() const { return num_reg_const; }
    void numRegConst(int num) { num_reg_const = num; } 
public : 
    cgalCircleSurf(INTERFACE* intfc, SURFACE** surf) : 
	cgalSurf(intfc, surf) { num_reg_const = 150; }
    void addCgalConst();
    double distance(double*, double*, int); 
    ~cgalCircleSurf() {}
};

class cgalParaSurf : public cgalCircleSurf {
    int num_lines; 
    int num_cons; 
    double innerRad;
    void addHole(std::vector<bool>&, int&); 
    void getSpecialParaFromFile(std::ifstream&);
public : 
    cgalParaSurf(INTERFACE* intfc, SURFACE** surf) :  
		cgalCircleSurf(intfc, surf) {}
    void addCgalConst(); 
    ~cgalParaSurf() {}; 
};

#endif
