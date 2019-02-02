AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
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

TESTS += CollisionTests
check_PROGRAMS += CollisionTests

CollisionTests_SOURCES = CollisionTests.cpp

CollisionTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
						$(top_builddir)/lib/libCollision.la

$(CollisionTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check

FORCE:


tagsfile:
	ctags $(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
