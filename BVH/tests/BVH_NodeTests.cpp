#include <gmock/gmock.h>
#include <BVH_Node.h>

class BVH_NodeTests : public testing::Test
{
    protected:

    static POINT *a, *b, *c, *d, *e, *f, *g;
    static TRI *t1, *t2, *t3, *t4,*t5;
    static Hse *T1, *T2, *T3, *T4, *T5;

    BVH_Node *l1, *l2, *l3, *l4, *l5;

    static void SetUpTestCase()
    {
        a = new POINT;         b = new POINT;
        Coords(a)[0] = 0.0;    Coords(b)[0] = 1.0;
        Coords(a)[1] = 0.0;    Coords(b)[1] = 0.0;
        Coords(a)[2] = 0.0;    Coords(b)[2] = 0.0;

        c = new POINT;          d = new POINT;
        Coords(c)[0] = 0.0;     Coords(d)[0] = 0.0;
        Coords(c)[1] = 1.0;     Coords(d)[1] = 0.0;
        Coords(c)[2] = 0.0;     Coords(d)[2] = 1.0;

        e = new POINT;          f = new POINT;
        Coords(e)[0] = -1.0;     Coords(f)[0] = 0.0;
        Coords(e)[1] = 0.0;     Coords(f)[1] = -1.0;
        Coords(e)[2] = 0.0;     Coords(f)[2] = 0.0;

        g = new POINT;          //f = new POINT;
        Coords(g)[0] = 0.0;     //Coords(f)[0] = 1.0;
        Coords(g)[1] = 0.0;     //Coords(f)[1] = 1.0;
        Coords(g)[2] = -1.0;    //Coords(f)[2] = 0.0;

        t1 = new TRI;                t2 = new TRI;
        Point_of_tri(t1)[0] = a;     Point_of_tri(t2)[0] = b;
        Point_of_tri(t1)[1] = b;     Point_of_tri(t2)[1] = c;
        Point_of_tri(t1)[2] = c;     Point_of_tri(t2)[2] = d;

        t3 = new TRI;                t4 = new TRI;
        Point_of_tri(t3)[0] = c;     Point_of_tri(t4)[0] = d;
        Point_of_tri(t3)[1] = d;     Point_of_tri(t4)[1] = e;
        Point_of_tri(t3)[2] = e;     Point_of_tri(t4)[2] = f;

        t5 = new TRI;                //t6 = new TRI;
        Point_of_tri(t5)[0] = e;     //Point_of_tri(t6)[0] = f;
        Point_of_tri(t5)[1] = f;     //Point_of_tri(t6)[1] = g;
        Point_of_tri(t5)[2] = g;     //Point_of_tri(t6)[2] = h;

        T1 = new HsTri(t1);    T2 = new HsTri(t2);
        T3 = new HsTri(t3);    T4 = new HsTri(t4);
        T5 = new HsTri(t5);    //T6 = new HsTri(t6);
    }

    static void TearDownTestCase()
    {
        delete t1; delete t2; delete t3;
        delete t4; delete t5;
        delete T1; delete T2; delete T3;
        delete T4; delete T5;
        delete a; delete b; delete c;
        delete d; delete e; delete f;
        delete g;
    }

    void SetUp() override
    {
        l1 = new LeafNode(T1);
        l2 = new LeafNode(T2);
        l3 = new LeafNode(T3);
        l4 = new LeafNode(T4);
        l5 = new LeafNode(T5);
    }

    void TearDown() override
    {
        delete l1;   delete l2;
        delete l3;   delete l4;
        delete l5;

    }

    ~BVH_NodeTests() = default;
};

TRI* BVH_NodeTests::t1 = nullptr;
TRI* BVH_NodeTests::t2 = nullptr;
TRI* BVH_NodeTests::t3 = nullptr;
TRI* BVH_NodeTests::t4 = nullptr;
TRI* BVH_NodeTests::t5 = nullptr;
POINT* BVH_NodeTests::a = nullptr;
POINT* BVH_NodeTests::b = nullptr;
POINT* BVH_NodeTests::c = nullptr;
POINT* BVH_NodeTests::d = nullptr;
POINT* BVH_NodeTests::e = nullptr;
POINT* BVH_NodeTests::f = nullptr;
POINT* BVH_NodeTests::g = nullptr;
Hse* BVH_NodeTests::T1 = nullptr;
Hse* BVH_NodeTests::T2 = nullptr;
Hse* BVH_NodeTests::T3 = nullptr;
Hse* BVH_NodeTests::T4 = nullptr;
Hse* BVH_NodeTests::T5 = nullptr;



using DISABLED_BVH_NodeTests = BVH_NodeTests;



TEST_F(BVH_NodeTests, CheckIfLeafNodeIsSelfOrHasAdjacentHse)
{
    ASSERT_TRUE(l1->hasAdjacentHse(l1));
    ASSERT_TRUE(l1->hasAdjacentHse(l2));
    ASSERT_FALSE(l2->hasAdjacentHse(l5));
}

TEST_F(BVH_NodeTests, InternalNodeCtorDeathTest)
{
    BVH_Node* l6;
    ASSERT_DEATH(auto p = new InternalNode(l1,l6,1),"");
}


TEST_F(BVH_NodeTests, ConstructorLeafNode)
{
    ASSERT_TRUE(l5->isLeaf());
    ASSERT_NE(l5->getHse(),nullptr);
    ASSERT_EQ(l5->getParent(),nullptr);

    BoundingVolume bv5 = l5->getBV();
    ASSERT_DOUBLE_EQ(bv5.upper[0],0.0);
    ASSERT_DOUBLE_EQ(bv5.lower[1],-1.0);
}


