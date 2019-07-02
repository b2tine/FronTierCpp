AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(top_srcdir)/FrontWrapper \
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

TESTS += FrontWrapperTests
check_PROGRAMS += FrontWrapperTests

FrontWrapperTests_SOURCES = FrontWrapperTests.cpp

FrontWrapperTests_LDADD = $(top_builddir)/FrontWrapper/objs/libFronTierCpp.la \
						  $(top_builddir)/lib/gmock/lib/libgmock_main.la

$(FrontWrapperTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/FrontWrapper/objs && $(MAKE) libFronTierCpp.la
	@cd $(top_builddir)/lib/gmock && $(MAKE) check

FORCE:


tagsfile:
	ctags $(top_srcdir)/FrontWrapper/*.cpp $(top_srcdir)/FrontWrapper/*.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
