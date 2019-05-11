#include "BVH.h"
#include "CramersRule2d.h"

using NodePair = std::pair<BVH_Node*,BVH_Node*>;

static std::vector<NodePair> GetProximityCandidates(
        BVH_Node* const nodeA, BVH_Node* const nodeB);

static void ProcessProximityCandidates(std::vector<NodePair>& candidates);

static double HseToHseDistance(const Hse* A, const Hse* B);

static std::vector<double> ClosestPointOfTriToPoint(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts,
        double TOL = HUGE);

static std::vector<double> ClosestPointOfEdgeToPoint(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& edgePts);

/*
static std::vector<double> PointToClosestPointOfEdgeVec(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& edgePts);
*/

static double TriToTriDistance(const std::vector<POINT*>& ptsA,
                               const std::vector<POINT*>& ptsB);


//TODO: Move these primitive functions into their own file
//      and make globally available.

static std::vector<double> Pt2Vec(const POINT* p);

static std::vector<double> Pts2Vec(const POINT* p1, const POINT* p2);

static std::vector<double> CrossVec(const std::vector<double>& u,
                                    const std::vector<double>& v);

static std::vector<double> ScalarVec(double c, const std::vector<double>& u);

static std::vector<double> AddVec(const std::vector<double>& u,
                                  const std::vector<double>& v);

static std::vector<double> MinusVec(const std::vector<double>& u,
                                    const std::vector<double>& v);

static std::vector<double> NormalizeVec(std::vector<double>& u);

static double DotVec(const std::vector<double>& u,
                     const std::vector<double>& v);

static double MagVec(const std::vector<double>& v);

static int LargestComponentIndexVec(const std::vector<double>& v);

static std::pair<std::vector<double>, std::vector<std::vector<double>> >
ProjectOutComponentPointAndTri(int index,
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts);

static bool PointInTri(const std::vector<double>& pp,
        const std::vector<std::vector<double>>& ptriPts);

static std::pair<std::vector<int>,int>
PointToTriTangencyDecomposition(const std::vector<double>& pp,
        const std::vector<std::vector<double>>& ptriPts);

static bool LeftTurn(const std::vector<double>& a,
                     const std::vector<double>& b,
                     const std::vector<double>& c);

static bool CollinearOrLeftTurn(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c);

static double SignedParallelogramArea(const std::vector<double>& a,
                                      const std::vector<double>& b,
                                      const std::vector<double>& c);

//NOTE: checkProximity() is the entry point for simulation runs

//TODO: Return type of this function is temporary for testing.
//      Need to actually perform distance computations for the
//      candidate Hse's corresponding to the nodes/BVs returned
//      in the stack by GetProximityCandidates().
const bool checkProximity(const BVH* A, const BVH* B)
{
    assert(A && B);
    BVH_Node* rootA = A->getRoot();
    BVH_Node* rootB = B->getRoot();
    assert(rootA && rootB);
    
    auto proximity_candidates = GetProximityCandidates(rootA,rootB);
    
    //TODO: Implement this
    //ProcessProximityCandidates(proximity_candidates);
    
    //NOTE: at this point the CollisionSolver can apply the
    //      necessary repulsion forces to the nodes remaining
    //      in the stack i.e. the ones that are actually within
    //      proximity of each other.

    //Temp return type for testing.
    if( !proximity_candidates.empty() )
        return true;
    else
        return false;
}

///////////////////////////////////////////////
////      Functions below are static      ////
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

/*
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
*/

/*
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
*/

/*
double TriToTriDistance(
        const std::vector<POINT*>& ptsA,
        const std::vector<POINT*>& ptsB)
{
    //6 PointToTriDistance() checks
    //9 EdgeToEdgeDistance() checks
    
    //TODO: convert to vector<double> equivalents and pass into functions.
    //      Something like this...
    auto x1 = Pt2Vec(triPts[0]);
    auto x2 = Pt2Vec(triPts[1]);
    auto x3 = Pt2Vec(triPts[2]);
    std::vector<std::vector<double>> tpts = {x1,x2,x3};

    and a for loop with Pt2Vec(otherTri[i]) and tpts passed as args
    to PointToClosestPointOfTriVec()
}
*/

//If the point and triangle are within proximity (function of TOL)
//of each other, the returned vector corresponds to the closest point
//of the triangle. If they are not within proximity of each other, the
//returned vector is a default constructed (i.e. empty) std::vector<double>.
//Use the  vector<double>::empty() method to determine which is the case.
std::vector<double> ClosestPointOfTriToPoint(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts,
        double TOL)
{
    //NOTE: TOL has default value of HUGE, and if not specified
    //      the closest point of the triangle will be returned regardless
    //      of how far it is from the point, p.

    auto x12 = MinusVec(triPts[1],triPts[0]);
    auto x13 = MinusVec(triPts[2],triPts[0]);
    auto x14 = MinusVec(p,triPts[0]);
    
    //TODO: Need to check if triangle is degenerate here,
    //      or is fine inside of NormalizeVec()?
    //      Maybe just do both for safety and immediate clarity
    auto ntri = CrossVec(x12,x13);
    auto unitnormal = NormalizeVec(ntri);
    double signedDistToPlaneOfTri = DotVec(x14,unitnormal);
    
    //Check if the point, p, is close enough to the plane of the
    //triangle for it to be possible for the point, p, and the triangle
    //to be within the prescribed proximity tolerance of each other.
    if( fabs(signedDistToPlaneOfTri) > TOL )
    {
        //Return default constructed std::vector<double> if not.
        return {};
    }

    //Compute barycentric coordinates of the closest point in the plane
    //of the triangle by solving the following linear system.
    CramersRule2d Lsq;

    Lsq.setA( DotVec(x12,x12), DotVec(x12,x13),
              DotVec(x12,x13), DotVec(x13,x13) );
    
    Lsq.setRHS( DotVec(x12,x14), DotVec(x13,x14) );

    auto W = Lsq.solve();
    W.insert(W.begin(),1.0 - W[0] - W[1]);

    //TODO: Determine geometric significance of this calculation (delta).
    //      See wikipedia page on barycentric coordinates,
    //      https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    double triArea = 0.5*MagVec(ntri);
    double charLength = std::sqrt(triArea);
    double delta = TOL/charLength;

    //Check if the projected point in the plane is close enough to the
    //triangle (via its barycentric coordinates W[i]) for it to be
    //possible for the point, p, and the triangle to be within the
    //prescribed proximity tolerance of each other.
    for( int i = 0; i < 3; ++i )
    {
        if( W[i] < -1.0*delta || W[i] > 1.0 + delta )
        {
            //Return default constructed std::vector<double> if not.
            return {};
        }
    }

    //Compute the closest point in the plane of the triangle to p, projx4.
    auto v = AddVec(ScalarVec(W[1],x12),ScalarVec(W[2],x13));
    auto projx4 = AddVec(triPts[0],v);
    
    //Project out the dimension corresponding to the largest component
    //of the triangle's normal vector to obtain Pprojx4 and PtriPts.
    //Prevents degeneracy that can result using a single projection plane.
    int proj_index = LargestComponentIndexVec(ntri);
    auto ProjectedPointTriPair =
        ProjectOutComponentPointAndTri(proj_index,projx4,triPts);
    auto Pprojx4 = ProjectedPointTriPair.first;    
    auto PtriPts = ProjectedPointTriPair.second;    

    std::cout << "proj_index = " << proj_index << "\n\n";
    for( int i = 0; i < 2; ++i )
    {
        std::cout << "Pprojx4[" << i << "] = " << Pprojx4[i] << "\n"; 
    }
    std::cout << "\n\n";
    for( int k = 0; k < 3; ++k )
    {
        for( int i = 0; i < 2; ++i )
        {
            std::cout << "PtriPts[" << k << "][" << i << "] = ";
            std::cout << PtriPts[k][i] << "\n"; 
        }
        std::cout << "\n\n";
    }

    //Case 1: Pprojx4 is in the projected triangle interior,
    //        then projx4 is the closest point of the triangle.
    if( PointInTri(Pprojx4,PtriPts) )
    {
        return projx4;
        //auto PointToTriVec = MinusVec(projx4,p);
        //return PointToTriVec;
    }

    //Compute indices of tangent and nontangent points of the
    //projected triangle from Pprojx4's line of sight.
    auto tangencyDecomp = PointToTriTangencyDecomposition(Pprojx4,PtriPts);
    auto tanIndices = tangencyDecomp.first;
    int tan0 = tanIndices[0];
    int tan1 = tanIndices[1];
    int nontan = tangencyDecomp.second;

    //Case 2: The nontangent point of the projected triangle and Pprojx4
    //        are on same side of the edge joining the two tangent points,
    //        then the nontangent point's unprojected preimage is the closest
    //        point of the triangle.
    if( LeftTurn(PtriPts[tan0],PtriPts[tan1],Pprojx4) ==
            LeftTurn(PtriPts[tan0],PtriPts[tan1],PtriPts[nontan]) )
    {
        return triPts[nontan];
        //auto PointToTriVec = MinusVec(triPts[nontan],p);
        //return PointToTriVec;
    }

    //Case 3: The nontangent point of the projected triangle and Pprojx4 are
    //        on opposite sides of the edge joining the two tangent points.
    //        Find closest point to projx4 on the unprojected edge corresponding
    //        to the tangent points
    std::vector<std::vector<double>> closestTriEdge = {triPts[tan0],triPts[tan1]};
    return ClosestPointOfEdgeToPoint(projx4,closestTriEdge);

    //auto PointToTriVec = MinusVec(ClosestPointOnEdge,p);
    //return PointToTriVec;
}

/*
//Not sure if we need this yet ...
std::vector<double> PointToClosestPointOfEdgeVec(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& edgePts)
{
    auto ClosestPointOnEdge = ClosestPointOfEdgeToPoint(p,edgePts);
}
*/

//For details of this implementation see
//http://geomalgorithms.com/a02-_lines.html#Distance-to-Ray-or-Segment
std::vector<double> ClosestPointOfEdgeToPoint(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& edgePts)
{
    auto u = MinusVec(p,edgePts[0]);
    auto v = MinusVec(edgePts[1],edgePts[0]);

    double uDotv = DotVec(u,v);
    if( uDotv <= 0.0 ) 
        return edgePts[0];

    double vsquared = DotVec(v,v);
    if( vsquared <= uDotv )
        return edgePts[1];

    auto projv_u = ScalarVec(uDotv/vsquared,v);
    auto point_of_proj = AddVec(edgePts[0],projv_u);
    return point_of_proj;
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

//Returns the vector u-v (i.e. the vector from v to u)
std::vector<double> MinusVec(const std::vector<double>& u,
                             const std::vector<double>& v)
{
    return AddVec(u,ScalarVec(-1.0,v));
}

std::vector<double> NormalizeVec(std::vector<double>& u)
{
    auto mag = MagVec(u);
    //TODO: Temporary assertion/tolerance for debugging.
    //      Unsure if this should be here or outside of function.
    //      Tolerance should be defined in a variable.
    assert(mag > 1.0e-12);
    return ScalarVec(1.0/mag,u);
}

int LargestComponentIndexVec(const std::vector<double>& v)
{
    int index = -1;
    double largest = 0.0;
    for( int i = 0; i < 3; ++i )
    {
        if( fabs(v[i]) > largest )
            index = i;
    }
    assert( index >= 0 );
    return index;
}

std::pair<std::vector<double>, std::vector<std::vector<double>> >
ProjectOutComponentPointAndTri(int index,
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts)
{
    std::vector<double> Pp;
    std::vector<std::vector<double>> PtriPts(3);
    for( int i = 0; i < 3; ++i )
    {
        if( i != index )
        {
            std::cout << "i != index \n\n";

            Pp.push_back(p[i]);
            for( int k = 0; k < 3; ++k )
            {
                PtriPts[k].push_back(triPts[k][i]);
            }
        }
    }
            
    std::cout << "Pp.size() = " << Pp.size() << "\n\n";

    Pp.shrink_to_fit();
    for( int k = 0; k < 3; ++k )
    {
        PtriPts[k].shrink_to_fit();
    }

    return std::make_pair(Pp,PtriPts);
}

//NOTE: Need to use ProjectOutComponentPointAndTri() before
//      calling any of these functions.
    
//TODO: Verify triPts given in CCW order for below functions
bool PointInTri(const std::vector<double>& pp,
        const std::vector<std::vector<double>>& ptriPts)
{
    for( int i = 0; i < 3; ++i )
    {
        int iplus1 = (i+1)%3;
        if( !LeftTurn(ptriPts[i],ptriPts[iplus1],pp) )
            return false;
    }
    return true;
}

std::pair<std::vector<int>,int>
PointToTriTangencyDecomposition(const std::vector<double>& pp,
        const std::vector<std::vector<double>>& ptriPts)
{
    std::vector<int> tangentPtsIndices;
    std::vector<int> nontangentPointIndex;

    for( int i = 0; i < 3; ++i )
    {
        int iminus1 = (i+2)%3;
        int iplus1 = (i+1)%3;

        if( CollinearOrLeftTurn(ptriPts[iminus1],ptriPts[i],pp) !=
                CollinearOrLeftTurn(ptriPts[i],ptriPts[iplus1],pp) )
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

//Magnitude of ab cross ac, where a, b and c are 2d points in the plane.
double SignedParallelogramArea(const std::vector<double>& a,
                               const std::vector<double>& b,
                               const std::vector<double>& c)
{
    return (b[0] - a[0])*(c[1] - a[1])
            - (c[0] - a[0])*(b[1] - a[1]);
}

//////////////////////////////////////////////////
//           Functions for Testing              //
//   Needed to test the above static methods    //
//////////////////////////////////////////////////

std::vector<double> TestPointToEdge(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& edgePts)
{
    return ClosestPointOfEdgeToPoint(p,edgePts);
}

std::vector<double> TestPointToTri(
        const std::vector<double>& p,
        const std::vector<std::vector<double>>& triPts)
{
    return ClosestPointOfTriToPoint(p,triPts);
}

