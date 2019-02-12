AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(top_srcdir)/BVH \
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

TESTS += HSETests
check_PROGRAMS += HSETests

HSETests_SOURCES = HSETests.cpp \
				   ../HyperSurfElement.h \
				   ../HyperSurfElement.cpp

HSETests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la

$(HSETests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check


TESTS += BV_PointTests
check_PROGRAMS += BV_PointTests

BV_PointTests_SOURCES = BV_PointTests.cpp \
						../BV_Point.h


BV_PointTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la

$(BV_PointTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check


TESTS += AABBTests
check_PROGRAMS += AABBTests

AABBTests_SOURCES = AABBTests.cpp

AABBTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
				  $(top_builddir)/lib/libBVH.la

$(AABBTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check
	@cd $(top_builddir)/BVH/objs && $(MAKE) all


TESTS += BVH_NodeTests
check_PROGRAMS += BVH_NodeTests

BVH_NodeTests_SOURCES = BVH_NodeTests.cpp

BVH_NodeTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
					  $(top_builddir)/lib/libBVH.la


$(BVH_NodeTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check
	@cd $(top_builddir)/BVH/objs && $(MAKE) all


TESTS += BVHTests
check_PROGRAMS += BVHTests

BVHTests_SOURCES = BVHTests.cpp

BVHTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
					   $(top_builddir)/lib/libBVH.la


$(BVHTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check
	@cd $(top_builddir)/BVH/objs && $(MAKE) all













FORCE:


tagsfile:
	ctags $(top_srcdir)/BVH/*.h $(top_srcdir)/BVH/*.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
