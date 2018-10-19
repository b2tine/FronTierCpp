#include "FronTier.h"
#include "dummyfront.h"
#include "gtest/gtest.h"


TEST(DummyFrontTests, CompileTest)
{
    cpp::Front front;
    front.testfunction();
}
