# Build rules for tests.
# TESTS_ENVIRONMENT: Set environment variables for the test run
# TESTS: Define programs run automatically by "make check"
# UNIT_TESTS: Define programs run automatically by "make check-unit-tests"
# check_PROGRAMS: Define programs built by "make check" but not necessarily run

TEST_ENVIRONMENT=
TESTS =
UNIT_TESTS=
check_PROGRAMS = 

UNIT_TESTS += dummytest
check_PROGRAMS += dummytest

dummytest_SOURCES = dummytest.cpp
dummytest_CPPFLAGS = $(GTEST_CPPFLAGS)
dummytest_CXXFLAGS = $(GTEST_CXXFLAGS)
dummytest_LDFLAGS = $(GTEST_LDFLAGS)
dummytest_LDADD = $(GTEST_LIBS) \
					 $(top_builddir)/lib/gmock/lib/libgmock_main.la \
					 $(top_builddir)/lib/gmock/lib/libgmock.la
