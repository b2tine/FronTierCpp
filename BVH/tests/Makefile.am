#SUBDIRS = ../objs ../../lib/gmock

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

TEST_LIBS = $(top_builddir)/BVH/objs/libBVH.la \
			$(top_builddir)/lib/gmock/lib/libgmock_main.la


# TESTS_ENVIRONMENT: Set environment variables for the test run
TESTS_ENVIRONMENT = 

# TESTS: Define programs run automatically by "make check"
TESTS =

# check_PROGRAMS: Define programs built by "make check" but not necessarily run
check_PROGRAMS =



TESTS += HseTests
check_PROGRAMS += HseTests

HseTests_SOURCES = HseTests.cpp

HseTests_LDADD = $(TEST_LIBS)


#NOTE: CGAL_Point.h and CGAL_Point.cpp no longer used

#TESTS += CGAL_PointTests
#check_PROGRAMS += CGAL_PointTests

#CGAL_PointTests_SOURCES = CGAL_PointTests.cpp \
#						  ../CGAL_Point.h \
#						  ../CGAL_Point.cpp

#CGAL_PointTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la



TESTS += AABBTests
check_PROGRAMS += AABBTests

AABBTests_SOURCES = AABBTests.cpp

AABBTests_LDADD = $(TEST_LIBS)


TESTS += BVH_NodeTests
check_PROGRAMS += BVH_NodeTests

BVH_NodeTests_SOURCES = BVH_NodeTests.cpp

BVH_NodeTests_LDADD = $(TEST_LIBS)


TESTS += BVHTests
check_PROGRAMS += BVHTests

BVHTests_SOURCES = BVHTests.cpp \
				   BVH_Tester.h

BVHTests_LDADD = $(TEST_LIBS)


TESTS += QueryTests
check_PROGRAMS += QueryTests

QueryTests_SOURCES = QueryTests.cpp \
					 QueryTests.h

QueryTests_LDADD = $(TEST_LIBS)


$(TEST_LIBS): libraries

libraries:
	@$(MAKE) -C $(top_builddir)/BVH/objs libBVH.la
	@$(MAKE) -C $(top_builddir)/lib/gmock check



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/BVH/*.h $(top_srcdir)/BVH/*.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c


