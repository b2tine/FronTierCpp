AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(top_srcdir)/cFluid \
			 -I$(top_srcdir)/include \
			 @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_CXXFLAGS = $(GTEST_CXXFLAGS)


# Build rules for tests.
# TESTS_ENVIRONMENT: Set environment variables for the test run
# TESTS: Define programs run automatically by "make check"
# check_PROGRAMS: Define programs built by "make check" but not necessarily run

TESTS_ENVIRONMENT = 

TESTS =
UNIT_TESTS = 
check_PROGRAMS =

UNIT_TESTS += cFluidTests
check_PROGRAMS += cFluidTests

cFluidTests_SOURCES = cFluidTests.cpp

cFluidTests_LDADD = $(top_builddir)/cFluid/objs/libcFluid.la \
					$(top_builddir)/lib/gmock/lib/libgmock_main.la

TESTS += $(UNIT_TESTS)
TESTS += integrationtests.sh

.PHONY: check-unit-tests
check-unit-tests:
	$(MAKE) $(AM_MAKEFLAGS) check TESTS=$(UNIT_TESTS)


$(cFluidTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/cFluid/objs && $(MAKE) libcFluid.la
	@cd $(top_builddir)/cFluid && $(MAKE) cFluid
	@cd $(top_builddir)/lib/gmock && $(MAKE) check

FORCE:


tagsfile:
	ctags $(top_srcdir)/cFluid/*.h $(top_srcdir)/cFluid/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
