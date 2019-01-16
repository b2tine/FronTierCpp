#ifndef ORIGAMI_H
#define ORIGAMI_H

#include "spring_solver.h"
#include "drag.h"
#include "nlopt.hpp"
#include <unordered_map>
#include <armadillo>
#include "singleton.h"

namespace origamiSurface {
    typedef std::unordered_map<std::string, int> MapStrInt; 
    typedef enum {
        POLYGON,
        FACEONEARC,
        FACETWOARC
    } faceType;
}

class Crease; 

class Vertex {
    int index; 
    arma::vec coords;
    bool nv; 
public : 
    arma::vec getCoords() { return coords; }
    int getIdx() { return index; }
    bool judgeNonVirtual() { return nv; }
    Vertex(const std::vector<double>& cds, int idx, bool nonVir = false) : 
        coords(cds), index(idx), nv(nonVir) {};
};

class Nvvertex : public Vertex {
    std::vector<Crease*> crsFromVertex; 
public : 
    Nvvertex(const std::vector<double>& cds, int idx) : 
        Vertex(cds, idx, true) {}
    void insertCrease(Crease* crs) { crsFromVertex.push_back(crs); }
    std::vector<Crease*> getCreaseOnVertex() { return crsFromVertex; }
}; 

class Crease {
    Vertex* v1_; 
    Vertex* v2_; 
    // normalized vector v1v2
    arma::vec dir; 
    int index1; 
    int index2; 
    double rho_T;
    arma::mat rot_matrix; 
    void calRotMatrix(double);
public : 
    Crease(Vertex* v1, Vertex* v2, int idx1, int idx2, double); 
    void updateRotMatrix(double); 
    arma::mat getRotMatrix() { return rot_matrix; }
    arma::vec getDir() { return dir; }
    double getRotAngle() { return rho_T; }
};

class FaceType : public Singleton<FaceType>{
    origamiSurface::MapStrInt mymap;
protected:
    friend class Singleton<FaceType>;
    FaceType();
public:
    origamiSurface::MapStrInt& getMap() { return mymap; } 
};

class OgmPoint;

class Face {
    // vertices constructing the face
    std::vector<Vertex*> vertices_;
    // crease lines vertices on the face will go through
    // when do origami
    std::vector<Crease*> crsAlongPath;
    // folding matrix corresponds to crease lines l1, l2, ..., lk in 
    // crsAlongPath
    arma::mat fd_matrix; 
    origamiSurface::faceType type; 
protected: 
    void setType(origamiSurface::faceType ft) { type = ft; }
    std::vector<Vertex*> getVertices() { return vertices_; }
    bool leftOnStraightLine(arma::vec, arma::vec, arma::vec);
    bool leftOnArc(arma::vec, arma::vec, arma::vec, arma::vec);
public :
    Face(const std::vector<Vertex*>&, const std::vector<Crease*>&); 
    virtual bool poInside(OgmPoint*) = 0;
    std::vector<Crease*> getCrsAlongPath() { return crsAlongPath; }
    void crsFoldCrds(std::vector<double>&); 
    arma::mat& getFoldingMatrix() { return fd_matrix; }
    void updateFoldingMatrix(); 
    virtual ~Face() = 0; 
};

class Polygon : public Face {
public : 
    Polygon(const std::vector<Vertex*>& vv, const std::vector<Crease*>& cv) : 
        Face(vv, cv) { setType(origamiSurface::POLYGON); }
    virtual bool poInside(OgmPoint*); 
}; 

// the face with only 1 arc and others are straight lines
class FaceOneArc : public Face {
    arma::vec center; 
public : 
    FaceOneArc(const std::vector<Vertex*>& vv, const std::vector<Crease*>& cv, 
        const std::vector<double>& cen) : Face(vv, cv) { 
            center.resize(3); std::copy(cen.begin(), cen.end(), 
            center.begin()); setType(origamiSurface::FACEONEARC); }
    virtual bool poInside(OgmPoint*); 
};

class FaceTwoArc : public Face {
    arma::vec center; 
public : 
    FaceTwoArc(std::vector<Vertex*>& vv, const std::vector<Crease*>& cv, 
        const std::vector<double>& cen) : Face(vv, cv) {
        center.resize(3); std::copy(cen.begin(), cen.end(), center.begin()); 
        setType(origamiSurface::FACETWOARC); }
    virtual bool poInside(OgmPoint*); 
};

class OgmPoint : public SpringVertex
{
    std::vector<Face*> ogmFaces; 
    std::vector<double> x0;
public:
    OgmPoint(SpringVertex&);
    std::vector<Face*> getOgmFaces() { return ogmFaces; }
    void addOgmFace(Face* f) { ogmFaces.push_back(f); }
    std::vector<double> getInitialCoords() { return x0; }
};

class OptAlgorithm : public Singleton<OptAlgorithm> {
    origamiSurface::MapStrInt mymap;
protected:
    friend class Singleton<OptAlgorithm>;
    OptAlgorithm();
public : 
    origamiSurface::MapStrInt& getMap() { return mymap; }
};

class OrigamiFold : public Drag
{
    const double stepSize;
    const double wplus;
    const double wminus;
    double weight = 0.8;
    std::string outname; 
public:
    void preprocess(std::vector<SpringVertex*>&);
    void postprocess(std::vector<SpringVertex*>&);
    void setVel(SpringVertex*);
    void updateVel(std::vector<SpringVertex*>& pts, double t) {}
    void setAccel(SpringVertex*);
    std::string id() {return "OrigamiDrag";}
    OrigamiFold(std::string outname, const std::vector<std::vector<double>>&,
                const std::vector<double>&, const std::vector<int>&,
		const std::vector<std::vector<int>>&,
                const std::vector<std::pair<int, int>>&,
                const std::vector<std::vector<int>>&, 
                const std::vector<double>&, int);
    OrigamiFold();
    ~OrigamiFold();
    size_t dataSize();
private:
    bool first = true;
    std::vector<Vertex*> vertices_;
    std::vector<Crease*> creases_; 
    std::vector<Face*> faces_; 
    std::vector<double> rho, rho_delta, rho_tau, rho_T;
    nlopt::opt *m_opt = NULL;
    void assignFacesOgmPoints(std::vector<SpringVertex*>&); 
    std::vector<double> findFoldable(std::vector<double>&, double&);
    static double staticTargetFunction(const std::vector<double>&, std::vector<double>&, void*);
    double targetFunction(const std::vector<double>&);
    void findNextFoldingAngle();
    void ogmComputeNewPosition(SpringVertex* sv, std::vector<double>&);
    bool isValid(const std::vector<double>&);
    Drag * clone(const Info&);
    size_t totalDataSize;
    int optAlgoType_; 
};

#endif
