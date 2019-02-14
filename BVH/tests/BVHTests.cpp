#include <gmock/gmock.h>
#include <BVH.h>


class BVH_Tests : public testing::Test
{
    protected:

    static TRI *t1, *t2, *t3, *t4,*t5, *t6, *t7, *t8;
    static POINT *a, *b, *c, *d, *e, *f, *g, *h, *i, *j, *k, *l;

    static HsTri *T1, *T2, *T3, *T4, *T5, *T6, *T7, *T8;

    BVH bvh1;
    BVH bvh2;
    std::vector<Hse*> hseList1;
    std::vector<Hse*> hseList2;

    static void SetUpTestCase()
    {
        a = new POINT;         b = new POINT;         c = new POINT;
        Coords(a)[0] = 0.0;    Coords(b)[0] = 1.0;    Coords(b)[0] = 0.0;
        Coords(a)[1] = 0.0;    Coords(b)[1] = 0.0;    Coords(b)[0] = 1.0;
        Coords(a)[2] = 0.0;    Coords(b)[2] = 0.0;    Coords(b)[0] = 1.0;

        d = new POINT;         e = new POINT;         f = new POINT;
        Coords(c)[0] = -1.0;   Coords(d)[0] = -0.5;   Coords(f)[0] = -0.25;
        Coords(c)[1] = 1.0;    Coords(d)[1] = 2.0;    Coords(f)[1] = 2.5;
        Coords(c)[2] = 0.0;    Coords(d)[2] = 1.0;    Coords(f)[2] = 0.5;


        g = new POINT;         h = new POINT;         i = new POINT;
        Coords(g)[0] = 1.0;    Coords(h)[0] = 0.0;    Coords(i)[0] = 2.0;
        Coords(g)[1] = 1.0;    Coords(h)[1] = 1.0;    Coords(i)[1] = 0.0;
        Coords(g)[2] = 0.0;    Coords(h)[2] = 0.75;   Coords(i)[2] = 1.0;

        j = new POINT;         k = new POINT;         l = new POINT;
        Coords(j)[0] = 1.5;    Coords(k)[0] = 2.5;    Coords(l)[0] = 2.75;
        Coords(j)[1] = 1.5;    Coords(k)[1] = 1.0;    Coords(l)[1] = 5.0;
        Coords(j)[2] = 1.0;    Coords(k)[2] = 0.0;    Coords(l)[2] = 0.25;


        t1 = new TRI;                t2 = new TRI;              t3 = new TRI;
        Point_of_tri(t1)[0] = a;     Point_of_tri(t2)[0] = a;   Point_of_tri(t3)[0] = d;
        Point_of_tri(t1)[1] = b;     Point_of_tri(t2)[1] = c;   Point_of_tri(t3)[1] = c;
        Point_of_tri(t1)[2] = c;     Point_of_tri(t2)[2] = d;   Point_of_tri(t3)[2] = e;

        t4 = new TRI;                t5 = new TRI;              t6 = new TRI;
        Point_of_tri(t4)[0] = d;     Point_of_tri(t5)[0] = h;   Point_of_tri(t6)[0] = i;
        Point_of_tri(t4)[1] = c;     Point_of_tri(t5)[1] = g;   Point_of_tri(t6)[1] = j;
        Point_of_tri(t4)[2] = f;     Point_of_tri(t5)[2] = i;   Point_of_tri(t6)[2] = h;

        t7 = new TRI;                t8 = new TRI;              //t9 = new TRI;
        Point_of_tri(t7)[0] = i;     Point_of_tri(t8)[0] = k;   //Point_of_tri(t9)[0] = i;
        Point_of_tri(t7)[1] = k;     Point_of_tri(t8)[1] = l;   //Point_of_tri(t9)[1] = j;
        Point_of_tri(t7)[2] = j;     Point_of_tri(t8)[2] = j;   //Point_of_tri(t9)[2] = h;


        T1 = new HsTri(t1);  T2 = new HsTri(t2);  T3 = new HsTri(t3);
        T4 = new HsTri(t4);  T5 = new HsTri(t5);  T6 = new HsTri(t6);
        T7 = new HsTri(t7);  T8 = new HsTri(t8);  //T9 = new HsTri(t9);
    }

    static void TearDownTestCase()
    {
        delete t1; delete t2; delete t3;
        delete t4; delete t5; delete t6;
        delete t7; delete t8; //delete t9;
        delete T1; delete T2; delete T3;
        delete T4; delete T5; delete T6;
        delete T7; delete T8; //delete T9;
        delete a; delete b; delete c; delete d;
        delete e; delete f; delete g; delete h;
        delete i; delete j; delete k; delete l;
    }

    void SetUp() override
    {
        hseList1.push_back(T1);
        hseList1.push_back(T2);
        hseList1.push_back(T3);
        hseList1.push_back(T4);
        bvh1.buildTester(hseList1);
        
        hseList2.push_back(T5);
        hseList2.push_back(T6);
        hseList2.push_back(T7);
        hseList2.push_back(T8);
        bvh2.buildTester(hseList2);
    }

    void TearDown() override
    {

    }

    ~BVH_Tests() = default;
};


TRI* BVH_Tests::t1 = nullptr;
TRI* BVH_Tests::t2 = nullptr;
TRI* BVH_Tests::t3 = nullptr;
TRI* BVH_Tests::t4 = nullptr;
TRI* BVH_Tests::t5 = nullptr;
TRI* BVH_Tests::t6 = nullptr;
TRI* BVH_Tests::t7 = nullptr;
TRI* BVH_Tests::t8 = nullptr;
//TRI* BVH_Tests::t9 = nullptr;

POINT* BVH_Tests::a = nullptr;
POINT* BVH_Tests::b = nullptr;
POINT* BVH_Tests::c = nullptr;
POINT* BVH_Tests::d = nullptr;
POINT* BVH_Tests::e = nullptr;
POINT* BVH_Tests::f = nullptr;
POINT* BVH_Tests::g = nullptr;
POINT* BVH_Tests::h = nullptr;
POINT* BVH_Tests::i = nullptr;
POINT* BVH_Tests::j = nullptr;
POINT* BVH_Tests::k = nullptr;
POINT* BVH_Tests::l = nullptr;

HsTri* BVH_Tests::T1 = nullptr;
HsTri* BVH_Tests::T2 = nullptr;
HsTri* BVH_Tests::T3 = nullptr;
HsTri* BVH_Tests::T4 = nullptr;
HsTri* BVH_Tests::T5 = nullptr;
HsTri* BVH_Tests::T6 = nullptr;
HsTri* BVH_Tests::T7 = nullptr;
HsTri* BVH_Tests::T8 = nullptr;
//HsTri* BVH_Tests::T9 = nullptr;



using DISABLED_BVH_Tests = BVH_Tests;


TEST_F(DISABLED_BVH_Tests, AABBSelfProximity)
{

}


TEST_F(DISABLED_BVH_Tests, inProximityTest)
{
    //auto nodestack = queryProximity(&bvh1,&bvh2);
    //ASSERT_FALSE(nodestack.empty());
}


TEST_F(BVH_Tests, BuildTest)
{
    ASSERT_NE(bvh1.getRoot().lock(),nullptr);
}


