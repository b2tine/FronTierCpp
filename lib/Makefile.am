SUBDIRS = . gmock tests

AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

FT_LIBS = $(top_builddir)/src/front/objs/libfront.la \
		  $(top_builddir)/src/intfc/objs/libintfc.la \
		  $(top_builddir)/src/util/objs/libutil.la


lib_LTLIBRARIES = libFronTier.la

libFronTier_la_SOURCES = ../include/FronTier.h \
						 dummy.cpp

libFronTier_la_LIBADD = $(FT_LIBS)



$(FT_LIBS): libraries

libraries:
	@$(MAKE) -C $(top_builddir)/src libraries
