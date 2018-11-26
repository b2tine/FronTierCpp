AM_CPPFLAGS = $(GTEST_CPPFLAGS)
AM_CXXFLAGS = $(GTEST_CXXFLAGS)
AM_LDFLAGS = $(GTEST_LDFLAGS) $(GTEST_LIBS)


check_HEADERS = gmock/gmock.h \
				gtest/gtest.h


check_LTLIBRARIES = lib/libgmock.la \
					lib/libgmock_main.la


lib_libgmock_la_SOURCES = src/gmock-gtest-all.cc


lib_libgmock_main_la_SOURCES = src/gmock_main.cc
lib_libgmock_main_la_LIBADD = lib/libgmock.la
