#include <gmock/gmock.h>
#include "QueryTests.h"


class QueryTestData : public testing::Test
{
    protected:

    std::vector<double> O = {0,0,0};
    std::vector<double> I = {1,0,0};
    std::vector<double> J = {0,1,0};
    std::vector<double> K = {0,0,1};
    
    std::vector<double> p1 = {0,0,-1};
    std::vector<double> p2 = {1,1,1};

    std::vector<std::vector<double>> edge12 = {p1,p2}; 
    std::vector<double> c = {0,3,0};
    std::vector<double> d = {0,4,0};
    std::vector<double> dd = {0,4.1,0};

    std::vector<std::vector<double>> tri_OIJ = {O,I,J}; 
    std::vector<double> e = {0.25,0.25,1};

    std::vector<std::vector<double>> tri_OJK = {O,J,K}; 
    std::vector<std::vector<double>> tri_Op2J = {O,p2,J}; 
    std::vector<double> f = {5.0,0.5,0.25};
    std::vector<double> g = {2.0,0.5,0.25};


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


TEST_F(QueryTests, ClosestPointOfEdgeToPointTest)
{
    //Projected point is an interior point of edge
    auto ClosestPointToOrigin = TestPointToEdge(O,edge12);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[0],1.0/3.0);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[1],1.0/3.0);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[2],-1.0/3.0);

    //Projected point is an interior point of edge
    auto ClosestPointToC = TestPointToEdge(c,edge12);
    ASSERT_DOUBLE_EQ(ClosestPointToC[0],5.0/6.0);
    ASSERT_DOUBLE_EQ(ClosestPointToC[1],5.0/6.0);
    ASSERT_DOUBLE_EQ(ClosestPointToC[2],2.0/3.0);

    //Projected point is an endpoint of edge
    auto ClosestPointToD = TestPointToEdge(d,edge12);
    ASSERT_DOUBLE_EQ(ClosestPointToD[0],p2[0]);
    ASSERT_DOUBLE_EQ(ClosestPointToD[1],p2[1]);
    ASSERT_DOUBLE_EQ(ClosestPointToD[2],p2[2]);

    //Projected point is outside of edge
    auto ClosestPointToDD = TestPointToEdge(dd,edge12);
    ASSERT_DOUBLE_EQ(ClosestPointToDD[0],p2[0]);
    ASSERT_DOUBLE_EQ(ClosestPointToDD[1],p2[1]);
    ASSERT_DOUBLE_EQ(ClosestPointToDD[2],p2[2]);
}

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

    //Project to xy plane (tie goes to first one encountered)
    auto ClosestPointOp2J_ToF = TestPointToTri(f,tri_Op2J);
    ASSERT_DOUBLE_EQ(ClosestPointOp2J_ToF[0],p2[0]);
    ASSERT_DOUBLE_EQ(ClosestPointOp2J_ToF[1],p2[1]);
    ASSERT_DOUBLE_EQ(ClosestPointOp2J_ToF[2],p2[2]);

    //Project to xy plane (tie goes to first one encountered)
    auto ClosestPointOp2J_ToG = TestPointToTri(g,tri_Op2J);
    ASSERT_DOUBLE_EQ(ClosestPointOp2J_ToG[0],11.0/12.0);
    ASSERT_DOUBLE_EQ(ClosestPointOp2J_ToG[1],11.0/12.0);
    ASSERT_DOUBLE_EQ(ClosestPointOp2J_ToG[2],11.0/12.0);
}

