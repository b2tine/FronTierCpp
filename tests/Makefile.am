AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(abs_top_srcdir)/include \
			 -I$(abs_top_srcdir)/include/front \
			 -I$(abs_top_srcdir)/include/intfc \
			 -I$(abs_top_srcdir)/include/util @HDF_CFLAGS@

AM_CXXFLAGS = $(GTEST_CXXFLAGS)

AM_LDFLAGS = $(GTEST_LDFLAGS) $(GTEST_LIBS) -lgmock_main -lgmock \
			 -L$(abs_top_builddir)/lib -lfront -lintfc -lutil @HDF_LIBS@


# Build rules for tests.
# TESTS_ENVIRONMENT: Set environment variables for the test run
# TESTS: Define programs run automatically by "make check"
# check_PROGRAMS: Define programs built by "make check" but not necessarily run

TESTS_ENVIRONMENT = 

TESTS =
check_PROGRAMS =

TESTS += DummyFrontTests
check_PROGRAMS += DummyFrontTests

DummyFrontTests_SOURCES = DummyFrontTests.cpp
