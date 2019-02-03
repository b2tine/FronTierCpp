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

TESTS += intfcTests
check_PROGRAMS += intfcTests

intfcTests_SOURCES = intfcTests.cpp

intfcTests_LDADD = $(top_builddir)/lib/gmock/lib/libgmock_main.la \
				  $(top_builddir)/src/intfc/objs/libintfc.la


$(intfcTests_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib/gmock && $(MAKE) check
	@cd $(top_builddir)/src/intfc/objs && $(MAKE) all

FORCE:


tagsfile:
	ctags $(top_srcdir)/src/intfc/*.h $(top_srcdir)/src/intfc/*.c
