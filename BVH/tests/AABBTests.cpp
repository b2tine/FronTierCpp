#include <gmock/gmock.h>
#include <BoundingVolume.h>


class AABBTests : public testing::Test
{
    protected:

    static TRI *t1, *t2, *t3, *t4;
    static BOND *s1, *s2, *s3;
    static POINT *a, *b, *c, *d, *e, *f, *g, *h;

    static HsTri *T1, *T2, *T3, *T4;
    static HsBond *B1, *B2, *B3;
    static HsPoint* P;

    AABB bbT1, bbT2, bbT3, bbT4, bbB1, bbB2, bbB3;//, bbBVPts;

    static void SetUpTestCase()
    {
        a = new POINT;         b = new POINT;
        Coords(a)[0] = 2.0;    Coords(b)[0] = 1.0;
        Coords(a)[1] = 2.0;    Coords(b)[1] = 1.0;
        Coords(a)[2] = 3.0;    Coords(b)[2] = 1.0;

        c = new POINT;          d = new POINT;
        Coords(c)[0] = 9.0;     Coords(d)[0] = 10.0;
        Coords(c)[1] = 0.5;     Coords(d)[1] = 0.0;
        Coords(c)[2] = 2.0;     Coords(d)[2] = 0.0;

        e = new POINT;          f = new POINT;
        Coords(e)[0] = 0.0;     Coords(f)[0] = 0.0;
        Coords(e)[1] = 10.0;    Coords(f)[1] = 0.0;
        Coords(e)[2] = 0.0;     Coords(f)[2] = 10.0;

        g = new POINT;          h = new POINT;
        Coords(g)[0] = 0.0;     Coords(h)[0] = -1.0;
        Coords(g)[1] = 0.0;     Coords(h)[1] = -1.0;
        Coords(g)[2] = -10.0;   Coords(h)[2] = -1.0;

        t1 = new TRI;               t2 = new TRI;
        Point_of_tri(t1)[0] = a;    Point_of_tri(t2)[0] = d;
        Point_of_tri(t1)[1] = b;    Point_of_tri(t2)[1] = e;
        Point_of_tri(t1)[2] = c;    Point_of_tri(t2)[2] = f;

        t3 = new TRI;               t4 = new TRI;
        Point_of_tri(t3)[0] = d;    Point_of_tri(t4)[0] = d;
        Point_of_tri(t3)[1] = e;    Point_of_tri(t4)[1] = e;
        Point_of_tri(t3)[2] = g;    Point_of_tri(t4)[2] = b;

        s1 = new BOND;  s2 = new BOND;  s3 = new BOND;
        s1->start = a;  s2->start = c;  s3->start = h;
        s1->end = b;    s2->end = f;    s3->end = b;

        T1 = new HsTri(t1);     T2 = new HsTri(t2);
        T3 = new HsTri(t3);     T4 = new HsTri(t4);

        B1 = new HsBond(s1);   B2 = new HsBond(s2);
        B3 = new HsBond(s3);
        
        P = new HsPoint(a);
    }

    static void TearDownTestCase()
    {
        delete a; delete b; delete c;
        delete d; delete e; delete f;
        delete g; delete h;
        delete s1; delete s2; delete s3;
        delete t1; delete t2; delete t3;
        delete t4;
        delete B1; delete B2; delete B3;
        delete T1; delete T2; delete T3;
        delete T4;
        delete P;
    }

    void SetUp() override
    {
        bbT1 = AABB(T1);
        bbT2 = AABB(T2);
        bbT3 = AABB(T3);
        bbB1 = AABB(B1);
        bbB2 = AABB(B2);
        bbB3 = AABB(B3);
    }
  
    ~AABBTests() = default;

};

TRI* AABBTests::t1 = nullptr;
TRI* AABBTests::t2 = nullptr;
TRI* AABBTests::t3 = nullptr;
TRI* AABBTests::t4 = nullptr;
BOND* AABBTests::s1 = nullptr;
BOND* AABBTests::s2 = nullptr;
BOND* AABBTests::s3 = nullptr;
POINT* AABBTests::a = nullptr;
POINT* AABBTests::b = nullptr;
POINT* AABBTests::c = nullptr;
POINT* AABBTests::d = nullptr;
POINT* AABBTests::e = nullptr;
POINT* AABBTests::f = nullptr;
POINT* AABBTests::g = nullptr;
POINT* AABBTests::h = nullptr;
HsTri* AABBTests::T1 = nullptr;
HsTri* AABBTests::T2 = nullptr;
HsTri* AABBTests::T3 = nullptr;
HsTri* AABBTests::T4 = nullptr;
HsBond* AABBTests::B1 = nullptr;
HsBond* AABBTests::B2 = nullptr;
HsBond* AABBTests::B3 = nullptr;
HsPoint* AABBTests::P = nullptr;


using DISABLED_AABBTests = AABBTests;


TEST_F(DISABLED_AABBTests, BoxesOverlapVsContain)
{
    //contains is supposed to mean strictly contained
    EXPECT_TRUE(bbT2.contains(bbT1));
    EXPECT_TRUE(bbT1.volume() < bbT2.volume());

    //shared surface should count as overlap
    EXPECT_TRUE(bbT2.overlaps(bbT3));
    //but not containment
    EXPECT_FALSE(bbT2.contains(bbT4));
    //TODO: This may not be an issue if preconditions
    //      can be checked at node level.
}

TEST_F(AABBTests, ConstructorTwoAABBs)
{
    AABB parentbox(bbT2,bbT3);

    BV_Point lower = parentbox.lower;
    BV_Point upper = parentbox.upper;

    ASSERT_DOUBLE_EQ(lower[2],-10.0);
    ASSERT_DOUBLE_EQ(upper[2],10.0);
}

TEST_F(DISABLED_AABBTests, ConstructorTwoBV_Points)
{
    /*
    CGAL_Point centroid = bbBVPts.Centroid();
    ASSERT_DOUBLE_EQ(centroid.x(),0.0);
    ASSERT_DOUBLE_EQ(centroid.y(),0.0);
    ASSERT_DOUBLE_EQ(centroid.z(),0.0);
    */
}

TEST_F(AABBTests, ConstructorOneHse)
{
    ASSERT_DOUBLE_EQ(bbT1.upper[0],9.0);
    ASSERT_DOUBLE_EQ(bbT1.lower[2],1.0);
}



