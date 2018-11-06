AM_CPPFLAGS = -I$(top_srcdir)/include \
			  -I$(top_srcdir)/include/front \
			  -I$(top_srcdir)/include/intfc \
			  -I$(top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -static


noinst_LTLIBRARIES = libutil.la

include_HEADERS = ../../include/util/cdecs.h \
				 ../../include/util/fnamedebug.h \
				 ../../include/util/navdecs.h \
				 ../../include/util/plotdecs.h \
				 ../../include/util/uapi.h \
				 ../../include/util/uprotos.h \
				 ../../include/util/vmalloc.h

libutil_la_SOURCES = cleanup.c \
					 debug.c \
					 error.c \
					 fft.c \
					 fgetstrin.c \
					 fsort.c \
					 machine.c \
					 matrix.c \
					 other.c \
					 output.c \
					 ppsub.c \
					 quad.c \
					 random.c \
					 roots.c \
					 runga.c \
					 screen.c \
					 sigplot.c \
					 simpleio.c \
					 sphhar.c \
					 times.c \
					 uinit.c \
					 umap.c \
					 vectormalloc.c \
					 vtk.c \
					 dummyutil.cpp

libutil_la_LIBADD = @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

