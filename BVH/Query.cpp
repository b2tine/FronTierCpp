#include "BVH.h"
#include "LSQ.h"


using NodePair = std::pair<BVH_Node*,BVH_Node*>;


static std::vector<NodePair> GetProximityCandidates(
        BVH_Node* const nodeA, BVH_Node* const nodeB);

static void ProcessProximityCandidates(std::vector<NodePair>& candidates);

static double HseToHseDistance(Hse* A, Hse* B);

static double TriToTriDistance(std::vector<POINT*> ptsA, std::vector<POINT*> ptsB);

static double PointToTriDistance(POINT* p, std::vector<POINT*> TriPoints);

static double TriToTriDistance(const std::vector<POINT*>& ptsA,
                                const std::vector<POINT*>& ptsB);

static std::vector<double> Pt2Vec(const POINT* p)
static std::vector<double> Pts2Vec(const POINT* p1, const POINT* p2);

static double DotVec(const std::vector<double>& u, const std::vector<double>& v);
static std::vector<double> CrossVec(const std::vector<double>& u, const std::vector<double>& v);
static double MagVec(const std::vector<double>& v);
static std::vector<double> ScalarVec(const std::vector<double>& u);
static std::vector<double> NormalizeVec(std::vector<double>& u);


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
        
        double distAB = HseToHseDistance(ptsA,ptsB);
        //auto ptsA = hseA->getHsePoints();
        //auto ptsB = hseB->getHsePoints();
        //double distAB = TriToTriDistance(ptsA,ptsB);
    }

}

/*
double HseToHseDistance(Hse* A, Hse* B)
{
    int nA = A->num_pts();
    int nB = B->num_pts();
    auto ptsA->getHsePoints();
    auto ptsB->getHsePoints();

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
*/

double TriToTriDistance(
        const std::vector<POINT*>& ptsA,
        const std::vector<POINT*>& ptsB)
{
    //6 PointToTriDistance() checks
    //9 EdgeToEdgeDistance() checks

}

double PointToTriDistance(POINT* p, std::vector<POINT*> triPts)
{
    double TOL = 1.0e-06;

    auto x13 = Pts2Vec(triPts[0],triPts[2]);
    auto x23 = Pts2Vec(triPts[1],triPts[2]);
    auto x43 = Pts2Vec(p,triPts[2]);
    
    auto ntri = CrossVec(x13,x23);
    auto unormal = Normalize(ntri);
    double distToPlaneOfTri = fabs(DotVec(x4,unormal));
    
    if( distTriPlane > TOL )
        return -1;

    LeastSquares2d Lsq;

    Lsq.setA( DotVec(x13,x13), DotVec(x13,x23),
              DotVec(x23,x13), DotVec(x23,x23) );
    
    Lsq.setRHS( DotVec(x13,x43), DotVec(x23,x43) );

    auto W = Lsq.solve();
    W.push_back(1.0 - W[0] - W[1]);

    double triArea = 0.5*MagVec(ntri);
    double charLength = std::sqrt(triArea);
    double delta = TOL/charLength;

    for( int i = 0; i < 3; ++i )
    {
        if( W[i] < -1.0*delta || W[i] > 1.0 + delta )
            return -1;
    }

    //TODO: convert from barycentric coords back to cartesian
    //      and compute the distance from p.

}

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
        cu[i] /= c;
    return cu;
}

std::vector<double> NormalizeVec(std::vector<double>& u)
{
    return ScalarVec(1.0/MagVec(u),u);
}

