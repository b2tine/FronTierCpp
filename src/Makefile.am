SUBDIRS = util intfc front


AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


noinst_LTLIBRARIES = ../lib/libFronTier.la


___lib_libFronTier_la_SOURCES = $(top_srcdir)/include/FronTier.h


___lib_libFronTier_la_LIBADD = $(builddir)/front/objs/libfront.la \
							   $(builddir)/intfc/objs/libintfc.la \
							   $(builddir)/util/objs/libutil.la



