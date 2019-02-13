#include <gmock/gmock.h>
#include <CGAL_Point.h>


class CGAL_PointTests :  public ::testing::Test
{
    protected:

        CGAL_Point defaultPoint;
        CGAL_Point xyzTriplePoint;

        CGAL_Point p1, p2, p3, p4;

    public:

        CGAL_PointTests()
            : xyzTriplePoint(1,2,3),
            p1(0,0,0), p2(1,2,3),
            p3(0,2,3), p4(1,2,0)
        {}

        ~CGAL_PointTests() = default;

};

using DISABLED_CGAL_PointTests = CGAL_PointTests;


TEST_F(CGAL_PointTests, ZeroVectorByDefault)
{
    ASSERT_DOUBLE_EQ(defaultPoint[0],0.0);
    ASSERT_DOUBLE_EQ(defaultPoint[1],0.0);
    ASSERT_DOUBLE_EQ(defaultPoint[2],0.0);
}

TEST_F(CGAL_PointTests, xyzTripleConstructor)
{
    ASSERT_DOUBLE_EQ(xyzTriplePoint[0],1.0);
    ASSERT_DOUBLE_EQ(xyzTriplePoint[1],2.0);
    ASSERT_DOUBLE_EQ(xyzTriplePoint[2],3.0);
}

TEST_F(CGAL_PointTests, OutOfRangeDeathTest)
{
    ASSERT_DEATH(xyzTriplePoint[-1],"");
    ASSERT_DEATH(xyzTriplePoint[3],"");
}

TEST_F(CGAL_PointTests, xCoordsNotEqualLessThan)
{
    ASSERT_TRUE(p1 < p2);
}

TEST_F(CGAL_PointTests, xCoordsEqualLessThan)
{
    ASSERT_TRUE(p1 < p3);
}

TEST_F(CGAL_PointTests, xyCoordsBothEqualLessThan)
{
    ASSERT_TRUE(p4 < p2);
}

TEST_F(CGAL_PointTests, PointsAreEqualLessThan)
{
    ASSERT_FALSE(defaultPoint < p1);
}

TEST_F(CGAL_PointTests, IsValidKeyForStdMap)
{
    std::map<CGAL_Point,int> imap;
    ASSERT_TRUE(imap.empty());

    imap[p1] = 1;
    ASSERT_FALSE(imap.empty());
    ASSERT_EQ(imap[p1],1);
}

//TODO: Start test driving HilbertSortingTraits


