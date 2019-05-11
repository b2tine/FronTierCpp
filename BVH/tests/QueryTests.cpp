#include <gmock/gmock.h>
#include "QueryTests.h"


class QueryTestData : public testing::Test
{
    protected:

    std::vector<double> O = {0,0,0};
    std::vector<double> I = {1,0,0};
    std::vector<double> J = {0,1,0};
    std::vector<double> K = {0,0,1};
    
    std::vector<double> Kn = {0,0,-1};
    std::vector<double> U = {1,1,1};
    std::vector<double> Upert = {1.05,1,1};

    std::vector<std::vector<double>> edge_KnU = {Kn,U}; 
    std::vector<std::vector<double>> edge_JK = {J,K}; 
    std::vector<double> c = {0,3,0};
    std::vector<double> d = {0,4,0};
    std::vector<double> dpert = {0,4.1,0};

    std::vector<std::vector<double>> tri_OIJ = {O,I,J}; 
    std::vector<double> e = {0.25,0.25,1};

    std::vector<std::vector<double>> tri_OJK = {O,J,K}; 
    std::vector<std::vector<double>> tri_OUJ = {O,U,J}; 
    std::vector<std::vector<double>> tri_OUJpert = {O,Upert,J}; 
    std::vector<double> f = {5.0,0.5,0.25};
    std::vector<double> g = {2.0,0.5,0.25};
    std::vector<double> h = {-0.5,1.5,-0.5};


    QueryTestData() 
    {}

    virtual ~QueryTestData() override
    {}

};



class QueryTests : public QueryTestData
{
    protected:


    public:

    QueryTests()
    {}

    ~QueryTests() override
    {}

};



using DISABLED_QueryTests = QueryTests;

//TODO: Split into individual tests for each case.
TEST_F(QueryTests, ClosestPointOfTriToPointTest)
{
    //Project to xy plane, Ppoint in Ptriangle interior
    auto ClosestPointOIJ_ToE = TestPointToTri(e,tri_OIJ);
    ASSERT_DOUBLE_EQ(ClosestPointOIJ_ToE[0],e[0]);
    ASSERT_DOUBLE_EQ(ClosestPointOIJ_ToE[1],e[1]);
    ASSERT_DOUBLE_EQ(ClosestPointOIJ_ToE[2],0.0);

    //Project to yz plane, Ppoint in Ptriangle interior
    auto ClosestPointOJK_ToF = TestPointToTri(f,tri_OJK);
    ASSERT_DOUBLE_EQ(ClosestPointOJK_ToF[0],0.0);
    ASSERT_DOUBLE_EQ(ClosestPointOJK_ToF[1],f[1]);
    ASSERT_DOUBLE_EQ(ClosestPointOJK_ToF[2],f[2]);

    //Project to yz plane (tie goes to comp with lower index),
    //Ppoint is outside Ptriangle, on the opposite side of the
    //nontangent vertex, and outside of the near edge.
    auto ClosestPointOUJ_ToF = TestPointToTri(f,tri_OUJ);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToF[0],U[0]);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToF[1],U[1]);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToF[2],U[2]);

    //Project to xy plane, Ppoint is outside Ptriangle,
    //on the opposite side of the nontangent vertex, and
    //outside of the near edge.
    auto ClosestPointOUJpert_ToF = TestPointToTri(f,tri_OUJpert);
    ASSERT_DOUBLE_EQ(ClosestPointOUJpert_ToF[0],Upert[0]);
    ASSERT_DOUBLE_EQ(ClosestPointOUJpert_ToF[1],Upert[1]);
    ASSERT_DOUBLE_EQ(ClosestPointOUJpert_ToF[2],Upert[2]);

    //Project to yz plane (tie goes to comp with lower index),
    //Ppoint is outside Ptriangle, on the opposite side of the
    //nontangent vertex, and is an interior point of the near edge.
    auto ClosestPointOUJ_ToG = TestPointToTri(g,tri_OUJ);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToG[0],11.0/12.0);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToG[1],11.0/12.0);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToG[2],11.0/12.0);

    //Project to yz plane (tie goes to comp with lower index),
    //Ppoint is outside Ptriangle, and on the same side as the
    //nontangent vertex.
    auto ClosestPointOUJ_ToH = TestPointToTri(h,tri_OUJ);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToH[0],J[0]);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToH[1],J[1]);
    ASSERT_DOUBLE_EQ(ClosestPointOUJ_ToH[2],J[2]);

}

//TODO: Split into individual tests for each case.
TEST_F(QueryTests, ClosestPointOfEdgeToPointTest)
{
    //Projected point is an interior point of edge
    auto ClosestPointToOrigin = TestPointToEdge(O,edge_KnU);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[0],1.0/3.0);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[1],1.0/3.0);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[2],-1.0/3.0);

    //Projected point is an interior point of edge
    auto ClosestPointToC = TestPointToEdge(c,edge_KnU);
    ASSERT_DOUBLE_EQ(ClosestPointToC[0],5.0/6.0);
    ASSERT_DOUBLE_EQ(ClosestPointToC[1],5.0/6.0);
    ASSERT_DOUBLE_EQ(ClosestPointToC[2],2.0/3.0);

    //Projected point is an endpoint of edge
    auto ClosestPointToD = TestPointToEdge(d,edge_KnU);
    ASSERT_DOUBLE_EQ(ClosestPointToD[0],U[0]);
    ASSERT_DOUBLE_EQ(ClosestPointToD[1],U[1]);
    ASSERT_DOUBLE_EQ(ClosestPointToD[2],U[2]);

    //Projected point is outside of edge
    auto ClosestPointToDP = TestPointToEdge(dpert,edge_KnU);
    ASSERT_DOUBLE_EQ(ClosestPointToDP[0],U[0]);
    ASSERT_DOUBLE_EQ(ClosestPointToDP[1],U[1]);
    ASSERT_DOUBLE_EQ(ClosestPointToDP[2],U[2]);
}

TEST_F(QueryTests, ClosestPointPairLineToLineTest)
{
    auto ClosestPoints = TestLineToLine(edge_JK,edge_KnU);
    auto ClosestOnLineJK = ClosestPoints.first;
    auto ClosestOnLineKnU = ClosestPoints.second;
}
