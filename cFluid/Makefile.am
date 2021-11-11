SUBDIRS = objs . tests

input_directories = input-1d input-2d input-3d
loose_infiles = in-test in-rt2d in-rt3d0 in-rt3d1 in-rt3d2

EXTRA_DIST = $(input_directories) \
			 $(loose_infiles) \
			 cFweno-save.cpp


AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@



noinst_PROGRAMS = cFluid


cFluid_SOURCES = cFluid.cpp


cFluid_LDADD = $(builddir)/objs/libcFluid.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/CGAL_Interface/*.h $(top_srcdir)/CGAL_Interface/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

