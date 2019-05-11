#include <gmock/gmock.h>
#include "QueryTests.h"


class QueryTestData : public testing::Test
{
    protected:

    std::vector<double> O = {0,0,0};
    std::vector<double> a = {1,0,0};
    std::vector<double> b = {0,1,0};
    
    std::vector<double> p1 = {0,0,-1};
    std::vector<double> p2 = {1,1,1};

    std::vector<std::vector<double>> edge12 = {p1,p2}; 
    std::vector<double> c = {0,3,0};
    std::vector<double> d = {0,4,0};
    std::vector<double> dd = {0,4.1,0};

    std::vector<std::vector<double>> tri_Obc = {O,a,b}; 
    std::vector<double> e = {0.5,0.5,1};


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
    ASSERT_DOUBLE_EQ(ClosestPointToD[0],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToD[1],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToD[2],1.0);

    //Projected point is outside of edge
    auto ClosestPointToDD = TestPointToEdge(dd,edge12);
    ASSERT_DOUBLE_EQ(ClosestPointToDD[0],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToDD[1],1.0);
    ASSERT_DOUBLE_EQ(ClosestPointToDD[2],1.0);
}

TEST_F(DISABLED_QueryTests, ClosestPointOfTriToPointTest)
{
    //auto TestPointToTriVec();
}

