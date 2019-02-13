AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(top_srcdir)/BVH \
			 @CGAL_CFLAGS@ \
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
	@cd $(top_builddir)/BVH/objs && $(MAKE) all


TESTS += CGAL_PointTests
check_PROGRAMS += CGAL_PointTests

CGAL_PointTests_SOURCES = CGAL_PointTests.cpp \
						  ../CGAL_Point.h \
						  ../CGAL_Point.cpp

CGAL_PointTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la



TESTS += AABBTests
check_PROGRAMS += AABBTests

AABBTests_SOURCES = AABBTests.cpp

AABBTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
				  $(top_builddir)/lib/libBVH.la



TESTS += BVH_NodeTests
check_PROGRAMS += BVH_NodeTests

BVH_NodeTests_SOURCES = BVH_NodeTests.cpp

BVH_NodeTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
					  $(top_builddir)/lib/libBVH.la




TESTS += BVHTests
check_PROGRAMS += BVHTests

BVHTests_SOURCES = BVHTests.cpp

BVHTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
					   $(top_builddir)/lib/libBVH.la





FORCE:


tagsfile:
	ctags $(top_srcdir)/BVH/*.h $(top_srcdir)/BVH/*.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c


