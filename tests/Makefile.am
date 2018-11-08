AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(top_srcdir)/include \
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
		   				$(top_builddir)/lib/libFronTier.la \
						@MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

