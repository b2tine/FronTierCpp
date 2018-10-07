AM_CPPFLAGS= $(GTEST_CPPFLAGS) \
			 -I$(abs_top_srcdir)/include \
			 -I$(abs_top_srcdir)/include/front \
			 -I$(abs_top_srcdir)/include/intfc \
			 -I$(abs_top_srcdir)/include/util @HDF_CFLAGS@

AM_CXXFLAGS = $(GTEST_CXXFLAGS)

AM_LDFLAGS = $(GTEST_LDFLAGS) $(GTEST_LIBS) -lgmock_main -lgmock \
			 -L$(abs_top_builddir)/lib -lfront -lintfc -lutil @HDF_LIBS@

#$(abs_top_builddir)/lib/gmock/lib/libgmock_main.la \
#$(abs_top_builddir)/lib/gmock/lib/libgmock.la \
#$(abs_top_builddir)/lib/libfront.la \
#$(abs_top_builddir)/lib/libintfc.la \
#$(abs_top_builddir)/lib/libutil.la

# Build rules for tests.
# TESTS_ENVIRONMENT: Set environment variables for the test run
# TESTS: Define programs run automatically by "make check"
# UNIT_TESTS: Define programs run automatically by "make check-unit-tests"
# check_PROGRAMS: Define programs built by "make check" but not necessarily run

TESTS_ENVIRONMENT = 

TESTS =
check_PROGRAMS =

TESTS += RectGridTests
check_PROGRAMS += RectGridTests

RectGridTests_SOURCES = RectGridTests.cpp
