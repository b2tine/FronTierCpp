#include "origami.h"
#include <iostream>
#include <iterator>
#include <fstream>
#include <typeinfo>
#include <algorithm>
#include <deque>

static std::ostream& operator << (std::ostream& os, const std::vector<double>&);

//class OgmPoint
OgmPoint::OgmPoint(SpringVertex& sv) : SpringVertex(sv)
{
    x0.resize(3,0);
    std::copy(sv.getCoords(), sv.getCoords()+3, x0.begin());
}

// class optAlgoSingleton
OptAlgorithm::OptAlgorithm() {
    mymap.insert({"GN_DIRECT_L_RAND", 2});
    mymap.insert({"LN_COBYLA", 25});
    mymap.insert({"LN_NELDERMEAD", 28});
    mymap.insert({"LD_TNEWTON_PRECOND", 17});
}

//class faceTypeSingleton
FaceType::FaceType() {
    mymap.insert({"POLYGON", 0}); 
    mymap.insert({"FACEONEARC", 1});
    mymap.insert({"FACETWOARC", 2});
}

//class OrigamiFold
void OrigamiFold::preprocess(std::vector<SpringVertex*>& pts)
{
     if (first)
     { 
         first = false;
	 for (size_t i = 0; i < pts.size(); ++i)
         {
              SpringVertex* sv = new OgmPoint(*pts[i]);
              delete pts[i];
              pts[i] = sv;
              pts[i]->setRegistered();
         }
	 assignFacesOgmPoints(pts); 
     }
     findNextFoldingAngle(); 
     for (int i = 0; i < creases_.size(); i++) 
	  creases_[i]->updateRotMatrix(rho_delta[i]); 
     for (auto it : faces_) 
          it->updateFoldingMatrix(); 
}

void OrigamiFold::assignFacesOgmPoints(std::vector<SpringVertex*>& pts) {
    for (size_t i = 0; i < pts.size(); i++) {
         OgmPoint* op = static_cast<OgmPoint*>(pts[i]); 
         
         for (auto it : faces_) 
              if (it->poInside(op)) {
                  op->addOgmFace(it); 
              }
    }     
}

std::vector<double> OrigamiFold::findFoldable(std::vector<double>& input, 
double& err)
{
    if (m_opt == NULL)
    {
        size_t N = rho.size();
        nlopt::algorithm a = static_cast<nlopt::algorithm>(optAlgoType_); 
        //nlopt::algorithm a = nlopt::LN_COBYLA; 
        m_opt = new nlopt::opt(a, N);        
        this->m_opt->set_stopval(1e-6);
        this->m_opt->set_ftol_abs(1e-6);
        this->m_opt->set_xtol_rel(1e-6);
        this->m_opt->set_maxeval(3000);
        nlopt::srand(1234);
        std::vector<double> lowerbounds(N, -M_PI);
        std::vector<double> upperbounds(N, M_PI);

        for (size_t i = 0; i < N; ++i)
        {
            if (rho_T[i] > 0)
                lowerbounds[i] = 0;
            if (rho_T[i] < 0)
                upperbounds[i] = 0;
        }
        m_opt->set_lower_bounds(lowerbounds);
        m_opt->set_upper_bounds(upperbounds);
        this->m_opt->set_min_objective(staticTargetFunction, this);
    }
    std::vector<double> ans = m_opt->optimize(input);
    err = targetFunction(ans);
    return ans;
}

double OrigamiFold::staticTargetFunction(
                const std::vector<double>& x,
                std::vector<double> &grad,
                void *my_func_data)
{
    OrigamiFold* obj = static_cast<OrigamiFold*>(my_func_data);
    return obj->targetFunction(x);
}

double OrigamiFold::targetFunction(
                const std::vector<double>& x)
{
    double ans = 0;
    
    for (size_t i = 0; i < creases_.size(); i++) {
         creases_[i]->updateRotMatrix(x[i]); 
    }
    Nvvertex* nvId; 
    for (auto it : vertices_) {
         if (!it->judgeNonVirtual()) continue; 
         Nvvertex* nv = static_cast<Nvvertex*>(it); 

         std::vector<Crease*> crsOnVertex = nv->getCreaseOnVertex();    
         arma::mat I(4, 4, arma::fill::eye);
         arma::mat M(4, 4, arma::fill::eye);

         for (auto it1 : crsOnVertex) {

              M *= it1->getRotMatrix();
         } 
         M -= I;
         ans += arma::norm(M, "fro"); 
    }
    return ans;
}

void OrigamiFold::findNextFoldingAngle()
{
    const int maxIter = 5000;
    static int iter = 0;
    bool success = false;

    while (!success && iter++ < maxIter)
    {
        int N = rho.size();
        std::vector<double> rho_rand(N); 
        std::vector<double> dir(N);

        dir.reserve(N);
        // use lambda to specify how to generate value
        std::generate_n(rho_rand.begin(), N, 
            []() { return ((double)rand() / RAND_MAX * 2 - 1) * M_PI; });
        // capture this by reference for lambda use
        std::transform(rho_rand.begin(), rho_rand.end(), rho_T.begin(), dir.begin(), 
            [this](double r1, double r2) 
                { return (1 - weight) * r1 + weight * r2; });
        std::transform(rho_delta.begin(), rho_delta.end(), dir.begin(), rho_tau.begin(), 
            [this](double r, double d) { 
                double rho = r + stepSize * d;
                return std::min(M_PI, std::max(rho, -M_PI)); 
            });

        double err = 0;
        
        arma::vec vrt(rho_T);
        arma::vec vrd(rho_delta);

        rho = findFoldable(rho_tau, err);

        arma::vec vr(rho);
        arma::vec v1(vrt-vr);
        arma::vec v2(vrt-vrd);
        double d1 = arma::norm(v1);
        double d2 = arma::norm(v2);

        if (isValid(rho) && d1 < d2)
        {
            rho_delta = rho;
            weight = weight + wplus;
            success = true;
            std::cout << "rho = " << rho << std::endl;
            if (d2 < 1.0e-3) {
                std::cout << "convergence reached: err = 1e-3" 
                << "folding angle = [" << rho_delta << "]" << std::endl; 
                std::cout << "Folding process is terminated" << std::endl; 
                m_t = 0; 
            }
        }
        else
        {
            weight = weight - wminus;
            success = false;
        }
        weight = std::max(std::min(weight, 1.0), 0.0);
    }
    if (iter >= maxIter)
    {
        std::cout << "Maxium iteration reached\n"
        << "folding angle = [" << rho_delta << "]"
        << std::endl;
        std::cout << "err = " << targetFunction(rho_delta) << std::endl;
        std::cout << "Folding process is terminated"
        << std::endl;
        m_t = 0;
    }
}

bool OrigamiFold::isValid(const std::vector<double>& rho)
{
    for (auto a : rho)
    {
        if (fabs(a) > M_PI)
            return false;
    }
    return true;
}

void OrigamiFold::postprocess(std::vector<SpringVertex*>& pts) {}

void OrigamiFold:: setVel(SpringVertex* sv)
{
    std::vector<double> new_crds(3, 0);
    ogmComputeNewPosition(sv, new_crds);
    double* vel = sv->getVel();
    for (int i = 0; i < 3; ++i)
        vel[i] = (new_crds[i] - sv->getCoords()[i])/getTimeStepSize();

}

void OrigamiFold::setAccel(SpringVertex*){}
size_t OrigamiFold::dataSize()
{
    return totalDataSize;
}

Drag* OrigamiFold::clone(const Info & info)
{
    if (info.data().size() == 0)
    {
        std::cout << "Warning: OrigamiFold is not created, no data provided"
                  << std::endl;
        return NULL;
    }
    if (info.data().size() > 0 && info.data().size() < 2)
    {
        std::cout << "Warning: OrigamiFold is not created, "
                  << "insufficient data is given"
                  << std::endl;
    }

    std::vector<double> v = info.data();
    size_t n_pt = (size_t)v[0];
    size_t n_crs = (size_t)v[1];
    size_t n_fs = (size_t)v[2]; 
    size_t n_fsp = (size_t)v[3]; 
    size_t n_ncm = (size_t)v[4];
    std::vector<std::vector<double>> points(n_pt);
    std::vector<std::pair<int, int>> creases(n_crs);
    std::vector<std::vector<int>> faces(n_fs); 
    std::vector<int> typeIdx(n_fs); 
    std::vector<double> cen; 
    std::vector<std::vector<int>> mappings(n_fs); 
    std::vector<double> angles(n_crs);
    int optAlgoType; 

    totalDataSize = 5 + n_pt * 3 + n_crs * 2 + n_crs * 1 + 3*n_fs + n_fsp + 
        3 + n_ncm + 1;
    if(!validateData(info))
        return NULL;

    double *it = &(v.front()) + 5;
    for (int i = 0; i < n_pt; ++i)
    {
        points[i] = {*it, *(it+1), *(it+2)};
        it += 3;
    }
    for (int i = 0; i < n_crs; ++i)
    {
        creases[i] = std::make_pair(int(*it), int(*(it+1)));
        it += 2;
    }
    cen.resize(3); 
    std::copy(it, it+3, cen.begin()); 
    it += 3; 
    for (int i = 0; i < n_fs; i++) {
         typeIdx[i] = *(it++); 

	 std::deque<int> temp; 

         for (int j = 0; j < int(*it); j++) {
	      temp.push_back(int(*(it+j+1)));
         }
	 it = it + int(*it)+1; 
	 faces[i] = std::vector<int>(temp.begin(), temp.end());
    }
    for (int i = 0; i < n_fs; i++) {
         std::deque<int> temp; 

         for (int j = 0; j < int(*it); j++)
              temp.push_back(*(it+j+1)); 
         mappings[i] = std::vector<int>(temp.begin(), temp.end()); 
         it = it + int(*it) + 1; 
    }
    for (int i = 0; i < n_crs; ++i)
    {
        angles[i] = *it;
        it++;
    }
    optAlgoType = (int)*it; 
    return new OrigamiFold(points, cen, typeIdx, faces, creases, mappings, 
        angles, optAlgoType);
}

OrigamiFold::OrigamiFold(): m_opt(NULL), stepSize(0), wplus(0), wminus(0) {}

OrigamiFold::OrigamiFold(const std::vector<std::vector<double>>& points,
                         const std::vector<double>& cen, 
                         const std::vector<int>& typeIdx, 
			 const std::vector<std::vector<int>>& fs, 
                         const std::vector<std::pair<int, int>>& creases,
                         const std::vector<std::vector<int>>& mappings, 
                         const std::vector<double>& angles, int optAlgoType) :
                        stepSize(0.015), wplus(0.2), wminus(0.01)
{
    m_t = 3.0;
    m_opt = NULL;
    if (creases.size() != angles.size())
    {
        std::cerr << "Number of creases does not equal to number of angles"
                  << std::endl;
    }

    size_t N = creases.size();

    rho.resize(N, 0);
    rho_delta.resize(N, 0);
    rho_tau.resize(N);
    rho_T.reserve(N);
    vertices_.resize(points.size()); 
    for (int i = 0; i < points.size(); i++) 
	 vertices_[i] = NULL;  
    for (size_t i = 0; i < creases.size(); ++i)
    {
        int vindex = creases[i].first;
        int pindex = creases[i].second;

        if (vertices_[vindex] == NULL)
            vertices_[vindex] = new Nvvertex(points[vindex], vindex);
	if (vertices_[pindex] == NULL) 
	    vertices_[pindex] = new Vertex(points[pindex], pindex); 

        creases_.push_back(new Crease(vertices_[vindex], vertices_[pindex], 
            vindex, pindex, angles[i]));
    }
    for (size_t i = 0; i < points.size(); i++) {
         if (vertices_[i] == NULL) 
             vertices_[i] = new Vertex(points[i], i); 
    }
    for (auto it : creases_)
         rho_T.push_back(it->getRotAngle()); 
    faces_.reserve(fs.size());
    for (size_t i = 0; i < fs.size(); i++) {
	 std::vector<Vertex*> tempFaceVertex(fs[i].size()); 
         std::vector<Crease*> tempFaceCrease(mappings[i].size()); 
         
         for (int j = 0; j < mappings[i].size(); j++)
              tempFaceCrease[j] = creases_[mappings[i][j]];
	// tempFaceCrease.push_back(creases_[i]); 
	 for (size_t j = 0; j < fs[i].size(); j++) 
	      tempFaceVertex[j] = vertices_[fs[i][j]]; 
         switch (typeIdx[i]) {
            case origamiSurface::FACEONEARC: 
	        faces_.push_back(new FaceOneArc(tempFaceVertex, 
                    tempFaceCrease, cen));
                break;  
            case origamiSurface::POLYGON: 
                faces_.push_back(new Polygon(tempFaceVertex, tempFaceCrease)); 
         }
    }
    for (auto it : vertices_) {
         if (!it->judgeNonVirtual()) continue; 

         Nvvertex* nv = static_cast<Nvvertex*>(it); 

         for (int i = 0; i < creases.size(); i++) 
              if (creases[i].first == nv->getIdx())
                  nv->insertCrease(creases_[i]);         
    }
    optAlgoType_ = optAlgoType; 
}


void OrigamiFold::ogmComputeNewPosition(SpringVertex* sv, std::vector<double>& new_crds)
{
    OgmPoint* op = static_cast<OgmPoint*>(sv); 

    new_crds[0] = op->getInitialCoords()[0]; 
    new_crds[1] = op->getInitialCoords()[1];
    new_crds[2] = op->getInitialCoords()[2];
    // note: for the point on the crease, 
    // we only move it once (rotate along creases of previous face or
    // current face. The two motion is equivalent since that point is on
    // the crease line)
    if (op->getOgmFaces().size() == 0) return; 

    Face* bface = op->getOgmFaces().back(); 

    bface->crsFoldCrds(new_crds); 
}

OrigamiFold::~OrigamiFold() {
    if (this->m_opt)
    {
        delete this->m_opt;
        this->m_opt = NULL;
    }
    for (auto v : vertices_)
         delete v; 
    for (auto c : creases_) 
         delete c; 
    for (auto f : faces_) 
         delete f; 
}

// class Crease
Crease::Crease(Vertex* v1, Vertex* v2, int idx1, int idx2, double rho) {
    v1_ = v1; 
    v2_ = v2; 
    index1 = idx1; 
    index2 = idx2; 
    rho_T = rho; 
    rot_matrix = arma::mat(4, 4, arma::fill::zeros);
    dir = arma::normalise(v2_->getCoords()-v1_->getCoords()); 
}

void Crease::calRotMatrix(double rho) {
    double a = dir(0), b = dir(1), c = dir(2);
    arma::vec v = v2_->getCoords(); 
    double x = v(0), y = v(1), z = v(2);
    arma::mat T(4, 4, arma::fill::eye);
    
    T(0, 3) = -x, T(1, 3) = -y, T(2, 3) = -z;

    arma::mat Ti(4, 4, arma::fill::eye);

    Ti(0, 3) = x, Ti(1, 3) = y, Ti(2, 3) = z;

    arma::mat Rz(4, 4, arma::fill::zeros);
    
    Rz(0, 0) = cos(rho), Rz(0, 1) = -sin(rho);
    Rz(1, 0) = sin(rho), Rz(1, 1) = cos(rho);
    Rz(2, 2) = Rz(3, 3) = 1.0;

    arma::mat Ry(4, 4, arma::fill::zeros);
    double d = sqrt(b*b+c*c);

    Ry(0, 0) = d, Ry(0, 2) = -a;
    Ry(1, 1) = 1.0;
    Ry(2, 0) = a, Ry(2, 2) = d;
    Ry(3, 3) = 1.0;

    arma::mat Ryi(arma::trans(Ry));
    arma::mat Rx(4, 4, arma::fill::eye);
    arma::mat Rxi(4, 4, arma::fill::eye);

    if (fabs(d) > 1.0e-12) {
        Rx(1, 1) = c/d, Rx(1, 2) = -b/d;
        Rx(2, 1) = b/d, Rx(2, 2) = c/d;
        Rxi = arma::trans(Rx);
    }

    rot_matrix = Rx * T;
    rot_matrix = Ry * rot_matrix; 
    rot_matrix = Rz * rot_matrix; 
    rot_matrix = Ryi * rot_matrix; 
    rot_matrix = Rxi * rot_matrix;
    rot_matrix = Ti * rot_matrix; 
}

void Crease::updateRotMatrix(double rho) {
    calRotMatrix(rho); 
}

// class Face
Face::Face(const std::vector<Vertex*>& vv, const std::vector<Crease*>& vc) {
    vertices_ = vv; 
    crsAlongPath = vc; 
    fd_matrix = arma::mat(4, 4, arma::fill::zeros);
}

void Face::crsFoldCrds(std::vector<double>& new_crds) {
    std::vector<double> coords(new_crds); 
    
    coords.push_back(1.0);

    arma::vec vcoords(coords);
    arma::vec ans(fd_matrix * vcoords); 

    for (int i = 0 ;i < 3; i++) 
         new_crds[i] = ans(i); 
}


void Face::updateFoldingMatrix() {
    arma::mat M(4, 4, arma::fill::eye);

    for (auto it : crsAlongPath) 
         M *= it->getRotMatrix(); 
    fd_matrix = M;
}

Face::~Face() {}

bool Face::leftOnStraightLine(arma::vec vp1, arma::vec vp2, arma::vec vq) {
    arma::vec va(vp2 - vp1);
    arma::vec vb(vq - vp1);
    arma::vec ans = arma::cross(va, vb);

    if (fabs(ans(2)) < 1.0e-12 || ans(2) > 0)
        return true;
    else return false;
}

bool Face::leftOnArc(arma::vec vp1, arma::vec vp2, 
        arma::vec vcen, arma::vec vq) {
    double radius = arma::norm(vp1 - vcen);
    double dist = arma::norm(vq - vcen);

    if (dist < radius || fabs(dist-radius) < 1.0e-4)
        return true;
    return false; 
}

bool Polygon::poInside(OgmPoint* op) {
    std::vector<Vertex*> vertices = getVertices();
    int size = vertices.size();
    arma::vec vOpIniCoords(op->getInitialCoords());

    for (size_t i = 0; i < size-1; i++) {
         arma::vec v1(vertices[i]->getCoords());
         arma::vec v2(vertices[i+1]->getCoords());

         if (leftOnStraightLine(v1, v2, vOpIniCoords)) continue;
         else return false;
    }

    arma::vec v1(vertices[size-1]->getCoords());
    arma::vec v2(vertices[0]->getCoords());

    if (!leftOnStraightLine(v1, v2, vOpIniCoords))
        return false;
    return true;
}

bool FaceOneArc::poInside(OgmPoint* op) {
    std::vector<Vertex*> vertices = getVertices(); 
    int size = vertices.size();
    arma::vec vOpIniCoords(op->getInitialCoords());

    for (size_t i = 0; i < size-1; i++) {
         arma::vec v1(vertices[i]->getCoords());
         arma::vec v2(vertices[i+1]->getCoords());
         if (leftOnStraightLine(v1, v2, vOpIniCoords)) continue;
         else return false;
    }

    arma::vec v1(vertices[size-1]->getCoords());
    arma::vec v2(vertices[0]->getCoords());

    if (!leftOnArc(v1, v2, center, vOpIniCoords))
        return false;
    return true;
}

bool FaceTwoArc::poInside(OgmPoint* op) {
    std::vector<Vertex*> vertices = getVertices(); 
    arma::vec vOpIniCoords(op->getInitialCoords());

    if (vertices.size() != 4) 
        throw std::runtime_error("FACETWOARC must have for edges!");
    for (size_t i = 0; i < 4; i++) {
         arma::vec v1(vertices[i]->getCoords());
         arma::vec v2(vertices[i+1]->getCoords());

         if (i & 1) 
             if (leftOnArc(v1, v2, center, vOpIniCoords)) continue; 
             else return false;
         else 
             if (leftOnStraightLine(v1, v2, vOpIniCoords)) continue; 
            else return false; 
    }
    return true; 
}

std::ostream& operator << (std::ostream& os, const std::vector<double>& v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(os, " "));
    return os;
}
