AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
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
UNIT_TESTS =
check_PROGRAMS =

UNIT_TESTS += iFluidTests
check_PROGRAMS += iFluidTests

iFluidTests_SOURCES = iFluidTests.cpp

iFluidTests_LDADD = $(top_builddir)/iFluid/objs/libiFluid.la \
				   	$(top_builddir)/lib/gmock/lib/libgmock_main.la

TESTS += $(UNIT_TESTS)
TESTS += integrationtests.sh

.PHONY: check-unit-tests
check-unit-tests:
	$(MAKE) $(AM_MAKEFLAGS) check TESTS=$(UNIT_TESTS)


$(iFluidTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/iFluid/objs && $(MAKE) libiFluid.la
	@cd $(top_builddir)/iFluid && $(MAKE) iFluid
	@cd $(top_builddir)/lib/gmock && $(MAKE) check

FORCE:


tagsfile:
	ctags $(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
