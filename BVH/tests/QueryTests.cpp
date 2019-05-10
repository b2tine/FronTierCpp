#include <gmock/gmock.h>


class QueryTestData : public testing::Test
{
    protected:

    std::vector<double> a = {0,0,0};
    std::sdt::vector<vector<double>> edge = { {1,0,0}, {0,1,0} };

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
    auto TestPointToEdgePt(a,edge);
}


