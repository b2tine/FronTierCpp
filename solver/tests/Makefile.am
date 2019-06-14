AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
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

TESTS += solverTests
check_PROGRAMS += solverTests

solverTests_SOURCES = solverTests.cpp

solverTests_LDADD = $(top_builddir)/solver/objs/libsolver.la \
					$(top_builddir)/lib/gmock/lib/libgmock_main.la

$(solverTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/solver/objs && $(MAKE) libsolver.la
	@cd $(top_builddir)/lib/gmock && $(MAKE) check

FORCE:


tagsfile:
	ctags $(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
