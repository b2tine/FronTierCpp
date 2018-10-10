AM_CPPFLAGS = -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util @HDF_CFLAGS@

AM_LDFLAGS = -L$(abs_top_builddir)/lib -lfront -lintfc -lutil @HDF_LIBS@

noinst_PROGRAMS = example00

example00_SOURCES = example00.cpp
