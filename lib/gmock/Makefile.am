check_LTLIBRARIES = lib/libgmock.la \
					lib/libgmock_main.la

lib_libgmock_la_SOURCES = \
	gmock/gmock.h gmock/gtest.h gmock-gtest-all.cc

lib_libgmock_la_CPPFLAGS = \
	-I$(srcdir)/gmock -I$(srcdir)/gtest $(GTEST_CPPFLAGS)
lib_libgmock_la_CXXFLAGS = $(GTEST_CXXFLAGS)
lib_libgmock_la_LDFLAGS = $(GTEST_LDFLAGS)
lib_libgmock_la_LIBADD = $(GTEST_LIBS)

lib_libgmock_main_la_SOURCES = \
	gmock/gmock.h gmock/gtest.h gmock_main.cc

lib_libgmock_main_la_CPPFLAGS = \
	-I$(srcdir)/gmock -I$(srcdir)/gtest $(GTEST_CPPFLAGS)
lib_libgmock_main_la_CXXFLAGS = $(GTEST_CXXFLAGS)
lib_libgmock_main_la_LDFLAGS = $(GTEST_LDFLAGS)
lib_libgmock_main_la_LIBADD = $(GTEST_LIBS) $(srcdir)/lib/libgmock.la
