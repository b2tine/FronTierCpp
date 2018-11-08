#include <FronTier.h>
#include <front/dummyfront.hpp>
#include <gtest/gtest.h>


TEST(DummyFrontTests, CompileTest)
{
    cpp::Front front;
    front.testfunction();
}
