AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(top_srcdir)/airfoil \
			 -I$(top_srcdir)/Collision \
			 @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			 @ARMA_CFLAGS@ \
			 -I$(top_srcdir)/iFluid \
			 -I$(top_srcdir)/solver \
			 @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			 -I$(top_srcdir)/include \
			 @HDF_CFLAGS@ @GD_CFLAGS@

AM_CXXFLAGS = $(GTEST_CXXFLAGS)


# Build rules for tests.
# TESTS_ENVIRONMENT: Set environment variables for the test run
# TESTS: Define programs run automatically by "make check"
# check_PROGRAMS: Define programs built by "make check" but not necessarily run

TESTS_ENVIRONMENT = 

TESTS =
check_PROGRAMS =

TESTS += airfoilTests
check_PROGRAMS += airfoilTests

airfoilTests_SOURCES = airfoilTests.cpp

airfoilTests_LDADD = $(top_builddir)/airfoil/objs/libairfoil.la \
					 $(top_builddir)/lib/gmock/lib/libgmock_main.la

$(airfoilTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/airfoil/objs && $(MAKE) libairfoil.la
	@cd $(top_builddir)/lib/gmock && $(MAKE) check

FORCE:


tagsfile:
	ctags $(top_srcdir)/airfoil/*.h $(top_srcdir)/airfoil/*.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
