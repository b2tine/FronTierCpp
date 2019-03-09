#include <gmock/gmock.h>
#include <HyperSurfElement.h>

class HseTestData : public testing::Test
{
    protected:

    POINT *a, *b, *c, *d, *e;
    BOND *s0, *s1, *s2;
    TRI *t0;

    HseTestData() 
        : a{new POINT}, b{new POINT}, c{new POINT},
        d{new POINT}, e{new POINT},
        s0{new BOND}, s1{new BOND}, s2{new BOND},
        t0{new TRI}
    {
    
    Coords(a)[0] = 1.0;     Coords(b)[0] = 0.0;     Coords(c)[0] = 0.0;
    Coords(a)[1] = 0.0;     Coords(b)[1] = 1.0;     Coords(c)[1] = 0.0;
    Coords(a)[2] = 0.0;     Coords(b)[2] = 0.0;     Coords(c)[2] = 1.0;

    Coords(d)[0] = 1.0;     Coords(e)[0] = 1.0;     //Coords(f)[0] = 0.0;
    Coords(d)[1] = 1.0;     Coords(e)[1] = 1.0;     //Coords(f)[1] = 0.0;
    Coords(d)[2] = 1.0;     Coords(e)[2] = 0.0;     //Coords(f)[2] = 1.0;

    Point_of_tri(t0)[0] = a;
    Point_of_tri(t0)[1] = b;
    Point_of_tri(t0)[2] = c;

    s0->start = b;      s1->start = d;      s2->start = e;
    s0->end = d;        s1->end = e;        s2->end = a;

    }

    virtual ~HseTestData() override
    {
        delete a;   delete b;   delete c;
        delete d;   delete e;   //delete f;
        delete s0;  delete s1;  delete s2;
        delete t0;
    }

};



class HseTests : public HseTestData
{
    protected:

    Hse *B0, *B1, *B2;
    Hse *T0;

    public:

    HseTests()
        : B0{new HsBond(s0)}, B1{new HsBond(s1)},
        B2{new HsBond(s2)}, T0{new HsTri(t0)}
    {}

    ~HseTests() override
    {
        delete B0;  delete B1;  delete B2;
        delete T0;
    }

};


TEST_F(HseTests, HsBondOutOfRangeDeathTest)
{
    ASSERT_DEATH(B0->Point_of_hse(2),"");
    ASSERT_DEATH(B0->Point_of_hse(-1),"");
}

TEST_F(HseTests, HsTriOutOfRangeDeathTest)
{
    ASSERT_DEATH(T0->Point_of_hse(3),"");
    ASSERT_DEATH(T0->Point_of_hse(-1),"");
}

TEST_F(HseTests, AdjacencyTest)
{
    ASSERT_TRUE(areAdjacentHse(T0,B0));
    ASSERT_TRUE(areAdjacentHse(B0,B1));
    ASSERT_FALSE(areAdjacentHse(T0,B1));
    ASSERT_TRUE(areAdjacentHse(B1,B2));
    ASSERT_TRUE(areAdjacentHse(B2,T0));
    ASSERT_FALSE(areAdjacentHse(B0,B2));
    ASSERT_TRUE(areAdjacentHse(T0,T0));
}





