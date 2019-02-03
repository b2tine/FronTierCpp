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

TESTS += utilTests
check_PROGRAMS += utilTests

utilTests_SOURCES = utilTests.cpp

utilTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
				  $(top_builddir)/src/util/objs/libutil.la


$(utilTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check
	@cd $(top_builddir)/src/util/objs && $(MAKE) all

FORCE:


tagsfile:
	ctags $(top_srcdir)/src/util/*.h $(top_srcdir)/src/util/*.c
