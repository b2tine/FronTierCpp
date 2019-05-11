#include <gmock/gmock.h>
#include "QueryTests.h"


class QueryTestData : public testing::Test
{
    protected:

    std::vector<double> origin = {0,0,0};
    
    std::vector<double> a = {0,0,-1};
    std::vector<double> b = {1,1,1};
    std::vector<std::vector<double>> edge_ab = {a,b}; 

    std::vector<double> c = {0,3,0};
    std::vector<double> d = {0,4,0};
    std::vector<double> e = {0,4.1,0};

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
    auto ClosestPointToOrigin = TestPointToEdgePt(origin,edge_ab);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[0],1.0/3.0);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[1],1.0/3.0);
    ASSERT_DOUBLE_EQ(ClosestPointToOrigin[2],-1.0/3.0);

    //Projected point is an interior point of edge
    auto ClosestPointToC = TestPointToEdgePt(c,edge_ab);
    ASSERT_DOUBLE_EQ(ClosestPointToC[0],5.0/6.0);
    ASSERT_DOUBLE_EQ(ClosestPointToC[1],5.0/6.0);
    ASSERT_DOUBLE_EQ(ClosestPointToC[2],2.0/3.0);

    //Projected point is an endpoint of edge
    auto ClosestPointToD = TestPointToEdgePt(d,edge_ab);
    ASSERT_DOUBLE_EQ(ClosestPointToD[0],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToD[1],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToD[2],1.0);

    //Projected point is not on the edge
    auto ClosestPointToE = TestPointToEdgePt(e,edge_ab);
    ASSERT_DOUBLE_EQ(ClosestPointToE[0],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToE[1],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToE[2],1.0);
}

TEST_F(DISABLED_QueryTests, PointToClosestPointOfTriVecTest)
{
    //auto TestPointToTriVec();
}

