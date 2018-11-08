AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(abs_top_srcdir)/include \
			 -I$(abs_top_srcdir)/include/front \
			 -I$(abs_top_srcdir)/include/intfc \
			 -I$(abs_top_srcdir)/include/util \
			 @MPI_CFLAGS@ @HDF_CFLAGS@

AM_CXXFLAGS = $(GTEST_CXXFLAGS)


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

DummyFrontTests_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) -lgmock_main -lgmock \
		   				-L$(top_builddir)/lib -lFronTier @MPI_LIBS@ @HDF_LIBS@

