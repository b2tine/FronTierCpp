AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -static
AM_LIBTOOLFLAGS = --preserve-dup-deps


noinst_HEADERS = \
../../../include/util/cdecs.h \
../../../include/util/fnamedebug.h \
../../../include/util/navdecs.h \
../../../include/util/plotdecs.h \
../../../include/util/uapi.h \
../../../include/util/uprotos.h \
../../../include/util/vmalloc.h


noinst_LTLIBRARIES = libutil.la

libutil_la_SOURCES = $(noinst_HEADERS)

libutil_la_SOURCES += \
../cleanup.c \
../debug.c \
../error.c \
../fft.c \
../fgetstrin.c \
../fsort.c \
../machine.c \
../matrix.c \
../other.c \
../output.c \
../ppsub.c \
../quad.c \
../random.c \
../roots.c \
../runga.c \
../screen.c \
../sigplot.c \
../simpleio.c \
../sphhar.c \
../times.c \
../uinit.c \
../umap.c \
../vectormalloc.c \
../vtk.c

#Add a dummy c++ source file so the linker behaves properly.
libutil_la_SOURCES += ../dummyutil.cpp

