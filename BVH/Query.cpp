#include "BVH.h"
#include "CramersRule2d.h"

using NodePair = std::pair<BVH_Node*,BVH_Node*>;

static std::vector<NodePair> GetProximityCandidates(
        BVH_Node* const nodeA, BVH_Node* const nodeB);

static void ProcessProximityCandidates(std::vector<NodePair>& candidates);

static double HseToHseDistance(const Hse* A, const Hse* B);

static std::vector<double>
PointToClosestPointOfTriVec(POINT* p, std::vector<POINT*> triPts);

static std::vector<double>
PointToClosestPointOfEdgeVec(POINT* p, std::vector<POINT*> edgePts);

static double TriToTriDistance(const std::vector<POINT*>& ptsA,
                               const std::vector<POINT*>& ptsB);

static std::vector<double> Pt2Vec(const POINT* p);
static std::vector<double> Pts2Vec(const POINT* p1, const POINT* p2);

static std::vector<double> CrossVec(const std::vector<double>& u,
                                    const std::vector<double>& v);

static std::vector<double> ScalarVec(const std::vector<double>& u);

static std::vector<double> AddVec(const std::vector<double>& u,
                                  const std::vector<double>& v);

static std::vector<double> MinusVec(const std::vector<double>& u,
                                    const std::vector<double>& v);

static std::vector<double> NormalizeVec(std::vector<double>& u);

static double DotVec(const std::vector<double>& u,
                     const std::vector<double>& v);
static double MagVec(const std::vector<double>& v);

static double SignedParallelogramArea(const std::vector<double>& a,
                                      const std::vector<double>& b,
                                      const std::vector<double>& c);

static bool LeftTurn(const std::vector<double>& a,
                     const std::vector<double>& b,
                     const std::vector<double>& c);

static bool CollinearOrLeftTurn(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c);

static bool PointInTri(const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts);

static std::pair<std::vector<int>,int>
PointToTriTangencyDecomposition(const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts)

//NOTE: checkProximity() is the only externally callable function.

//TODO: Return type of this function is temporary for testing.
//      Need to actually perform distance computations for the
//      candidate Hse's corresponding to the nodes/BVs returned
//      in the stack by queryProximity().
const bool checkProximity(const BVH* A, const BVH* B)
{
    assert(A && B);
    BVH_Node* rootA = A->getRoot();
    BVH_Node* rootB = B->getRoot();
    assert(rootA && rootB);
    
    auto proximity_candidates = GetProximityCandidates(rootA,rootB);
    
    //TODO: implement this
    //ProcessProximityCandidates(proximity_candidates);
    
    //NOTE: at this point the CollisionSolver can apply the
    //      necessary repulsion forces to the nodes remaining
    //      in the stack i.e. the ones that are actually in
    //      proximity of each other.
    //      


    if( !proximity_candidates.empty() )
        return true;
    else
        return false;
}

///////////////////////////////////////////////
////    All functions below are static    ////
////       i.e. local to this file       ////
////////////////////////////////////////////

std::vector<NodePair> GetProximityCandidates(
        BVH_Node* const nodeA,
        BVH_Node* const nodeB)
{
    std::stack<NodePair> qstack;
    qstack.push(std::make_pair(nodeA,nodeB));
    std::vector<NodePair> candidates;

    while( !qstack.empty() )
    {
        auto A = qstack.top().first;
        auto B = qstack.top().second;
        qstack.pop();

        if( A->overlaps(B) )
        {
            if( A->isLeaf() && B->isLeaf() )
            {
                if( A->hasAdjacentHse(B) )
                    continue;
                candidates.push_back(std::make_pair(A,B));
            }
            else if( A->isLeaf() )
            {
                auto rc = B->getRightChild();
                qstack.push(std::make_pair(A,rc));
                auto lc = B->getLeftChild();
                qstack.push(std::make_pair(A,lc));
            }
            else if( B->isLeaf() )
            {
                auto rc = A->getRightChild();
                qstack.push(std::make_pair(rc,B));
                auto lc = A->getLeftChild();
                qstack.push(std::make_pair(lc,B));
            }
            else
            {
                if( A->volume() < B->volume() )
                {
                    auto rc = B->getRightChild();
                    qstack.push(std::make_pair(A,rc));
                    auto lc = B->getLeftChild();
                    qstack.push(std::make_pair(A,lc));
                }
                else
                {
                    auto rc = A->getRightChild();
                    qstack.push(std::make_pair(rc,B));
                    auto lc = A->getLeftChild();
                    qstack.push(std::make_pair(lc,B));
                }
            }
        }
    }

    return candidates;
}

void ProcessProximityCandidates(std::vector<NodePair>& candidates)
{
    for( auto& pair : candidates  )
    {
        auto hseA = pair.first->getHse();
        auto hseB = pair.second->getHse();
        
        double distAB = HseToHseDistance(hseA,hseB);

            //auto ptsA = hseA->getHsePoints();
            //auto ptsB = hseB->getHsePoints();
            //double distAB = TriToTriDistance(ptsA,ptsB);
    }

}

double HseToHseDistance(const Hse* A, const Hse* B)
{
    int nA = A->num_pts();
    int nB = B->num_pts();

    //TODO: Verify these are given in CCW order.
    auto ptsA = A->getHsePoints();
    auto ptsB = B->getHsePoints();

    if( nA == 2 )
    {
        if( nB == 2 )
            return EdgeToEdgeDistance(ptsA,ptsB);
        else
            return EdgeToTriDistance(ptsA,ptsB);
    }
    else
    {
        if( nB == 2 )
            return EdgeToTriDistance(ptsB,ptsA);
        else
            return TriToTriDistance(ptsA,ptsB);
    }
}

/*
double TriToTriDistance(
        const std::vector<POINT*>& ptsA,
        const std::vector<POINT*>& ptsB)
{
    //6 PointToTriDistance() checks
    //9 EdgeToEdgeDistance() checks
}
*/

//TODO: Should we return a distance and vector pair?
//         i.e. std::pair<double,std::vector<double>>
//      Or should we just return the closest point, and
//      compute the rest external to this function?
//
//TODO: Const Correctness
std::vector<double> PointToClosestPointOfTriVec(
        POINT* p, std::vector<POINT*> triPts)
{
    //TODO: TOL should be given as argument or set within
    //      the class, that this becomes a member function of.
    //      The CollisionSolver class or some component member
    //      class thereof.
    double TOL = 1.0e-06;

    auto x12 = Pts2Vec(triPts[0],triPts[1]);
    auto x13 = Pts2Vec(triPts[0],triPts[2]);
    auto x14 = Pts2Vec(triPts[0],p]);
    
    //TODO: Need to check if triangle is degenerate?
    auto ntri = CrossVec(x13,x23);
    auto unormal = NormalizeVec(ntri);
    double distToPlaneOfTri = DotVec(x14,unormal);
    
    //Check if point is close enough to the plane to
    //warrant computing the actual distance to the triangle.
    if( fabs(distTriPlane) > TOL )
    {
        return -1;
    }

    //Correct the normal vector to point in the direction of p
    //TODO: is this necessary?
    if( distToPlaneOfTri < 0.0 )
    {
        unormal = ScalarVec(-1.0,unormal);
    }

    CramersRule2d Lsq;

    Lsq.setA( DotVec(x12,x12), DotVec(x12,x13),
              DotVec(x12,x13), DotVec(x13,x13) );
    
    Lsq.setRHS( DotVec(x12,x14), DotVec(x13,x14) );

    auto W = Lsq.solve();
    W.push_front(1.0 - W[0] - W[1]);

    //TODO: determine geometric significance of this calculation (delta)
    double triArea = 0.5*MagVec(ntri);
    double charLength = std::sqrt(triArea);
    double delta = TOL/charLength;

    for( int i = 0; i < 3; ++i )
    {
        if( W[i] < -1.0*delta || W[i] > 1.0 + delta )
            return -1;
    }

    auto x1 = Pt2Vec(triPts[0]);
    auto x2 = Pt2Vec(triPts[1]);
    auto x3 = Pt2Vec(triPts[2]);
    std::vector<std::vector<double>> tpts = {x1,x2,x3};

    //Compute the closest point in the plane to p
    auto v = AddVec(ScalarVec(W[1],x12),ScalarVec(W[2],x13));
    auto projx4 = AddVec(x1,v);
    
    //Case 1: projx4 is in the triangle interior, then
    //        it is the closest point of the triangle.
    if( PointInTri(projx4,tpts) )
    {
        auto PointToTriVec = MinusVec(projx4,Pt2Vec(p));
        return PointToTriVec;
        //double distance = MagVec(triToPointVec);
    }

    //Compute tangent points of triangle from p
    auto tangencyDecomp = PointToTriTangencyDecomposition(projx4,tpts);
    auto tanIndices = tangencyDecomp.first;
    int tan0 = tanIndices[0];
    int tan1 = tanIndices[1];
    int nontan = tangencyDecomp.second;

    //Case 2: The nontangent point and p are on same side of
    //        the edge joining the two tangent points, then the
    //        nontangent point is the closest point of the triangle.
    if( LeftTurn(tpts[tan0],tpts[tan1],projx4) ==
            LeftTurn(tpts[tan0],tpts[tan1],tpts[nontan]) )
    {
        auto PointToTriVec = MinusVec(tpts[nontan],Pt2Vec(p));
        return PointToTriVec;
        //double distance = MagVec(triToPointVec);
    }

    //TODO: Now implement
    //Case 3: The nontangent point and p are on opposite
    //        sides of the edge joining the two tangent points.

    std::vector<POINT*> tanPts = {triPts[tan0],triPts[tan1]};
    return PointToClosestPointOfEdgeVec(p,tanPts);

}

/*
std::vector<double> PointToClosestPointOfEdgeVec(
        POINT* p, std::vector<POINT*> edgePts)
{
}
*/

std::vector<double> Pt2Vec(const POINT* p)
{
    std::vector<double> v(3,0.0);
    for( int i = 0; i < 3; ++i )
        v[i] = Coords(p)[i];
    return v;
}

std::vector<double> Pts2Vec(const POINT* p1, const POINT* p2)
{
    std::vector<double> v12(3,0.0);
    for( int i = 0; i < 3; ++i )
        v12[i] = Coords(p2)[i] - Coords(p1)[i];
    return v12;
}

double MagVec(const std::vector<double>& v)
{
    double sqrMag = DotVec(v,v);
    return std::sqrt(sqrMag);
}

double DotVec(const std::vector<double>& u,
              const std::vector<double>& v)
{
    double val = 0.0;
    for( int i = 0; i < 3; ++i )
        val += u[i]*v[i];
    return val;
}

std::vector<double> CrossVec(const std::vector<double>& u,
                             const std::vector<double>& v)
{
    double n1 = u[1]*v[2] - u[2]*v[1];
    double n2 = u[2]*v[0] - u[0]*v[2];
    double n3 = u[0]*v[1] - u[1]*v[0];
    return std::vector<double> {n1, n2, n3};
}

std::vector<double> ScalarVec(double c, const std::vector<double>& u)
{
    std::vector<double> cu(u);
    for( int i = 0; i < 3; ++i )
        cu[i] *= c;
    return cu;
}

std::vector<double> AddVec(const std::vector<double>& u,
                           const std::vector<double>& v)
{
    std::vector<double> w(u);
    for( int i = 0; i < 3; ++i )
        w[i] += v[i];
    return w;
}

std::vector<double> MinusVec(const std::vector<double>& u,
                             const std::vector<double>& v)
{
    return AddVec(u,ScalarVec(-1.0,v));
}

std::vector<double> NormalizeVec(std::vector<double>& u)
{
    return ScalarVec(1.0/MagVec(u),u);
}

double SignedParallelogramArea(const std::vector<double>& a,
                               const std::vector<double>& b,
                               const std::vector<double>& c)
{
    return (b[0] - a[0])*(c[1] - a[1])
            - (c[0] - a[0])*(b[1] - a[1]);
}

bool LeftTurn(const std::vector<double>& a,
              const std::vector<double>& b,
              const std::vector<double>& c)
{
    return SignedParallelogramArea(a,b,c) > 0.0;
}

bool CollinearOrLeftTurn(const std::vector<double>& a,
                         const std::vector<double>& b,
                         const std::vector<double>& c)
{
    return SignedParallelogramArea(a,b,c) >= 0.0;
}

//TODO: Verify triPts given in CCW order for below functions
bool PointInTri(const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts)
{
    for( int i = 0; i < 3; ++i )
    {
        int iplus1 = (i+1)%3;
        if( !LeftTurn(triPts[i],triPts[iplus1],p) )
            return false;
    }
    return true;
}

std::pair<std::vector<int>,int>
PointToTriTangencyDecomposition(const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts)
{
    std::vector<int> tangentPtsIndices;
    std::vector<int> nontangentPointIndex;

    for( int i = 0; i < 3; ++i )
    {
        int iminus1 = (i+2)%3;
        int iplus1 = (i+1)%3;

        if( CollinearOrLeftTurn(triPts[iminus1],triPts[i],p)
                != CollinearOrLeftTurn(triPts[i],triPts[iplus1],p) )
        {
            tangentPtsIndices.push_back(i);
        }
        else
        {
            nontangentPointIndex.push_back(i);
        }
    }

    assert( nontangentPointIndex.size() == 1
            && tangentPtsIndices.size() == 2 );
    
    tangentPtsIndices.shrink_to_fit();
    nontangentPointIndex.shrink_to_fit();
    return std::make_pair(tangentPtsIndices,nontangentPointIndex[0]);
}


