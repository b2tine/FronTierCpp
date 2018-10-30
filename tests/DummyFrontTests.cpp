#include <FronTier.h>
#include <dummyfront.hpp>
#include <gtest/gtest.h>


TEST(DummyFrontTests, CompileTest)
{
    cpp::Front front;
    front.testfunction();
}
