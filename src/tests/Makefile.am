AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(top_srcdir)/include \
			 @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_CXXFLAGS = $(GTEST_CXXFLAGS)


# Build rules for tests.
# TESTS_ENVIRONMENT: Set environment variables for the test run
# TESTS: Define programs run automatically by "make check"
# check_PROGRAMS: Define programs built by "make check" but not necessarily run

TESTS_ENVIRONMENT = 

TESTS =
check_PROGRAMS =

TESTS += libFronTierTests
check_PROGRAMS += libFronTierTests

libFronTierTests_SOURCES = libFronTierTests.cpp

libFronTierTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
						 $(top_builddir)/lib/libFronTier.la


$(libFronTierTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check
	@cd $(top_builddir)/src && $(MAKE) all

FORCE:




tagsfile:
	ctags *.cpp $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
