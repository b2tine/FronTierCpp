#include <FronTier.h>
#include <front/dummyfront.hpp>
#include <gtest/gtest.h>

TEST(FrontTests, CompileTest)
{
    Front front;
}

TEST(cppFrontTests, CompileTest)
{
    cpp::Front front;
    front.testfunction();
}
