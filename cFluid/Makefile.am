SUBDIRS = objs


input_directories = input-1d input-2d input-3d
loose_infiles = in-test in-rt2d in-rt3d0 in-rt3d1 in-rt3d2

EXTRA_DIST = $(input_directories) \
			 $(loose_infiles) \
			 cFweno-save.cpp



AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = cFluid


cFluid_SOURCES = cFluid.cpp


cFluid_LDADD = $(top_builddir)/lib/libcFluid.la



tagsfile:
	ctags *.h *.cpp $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

