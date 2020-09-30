SUBDIRS = objs . tests

input_directories = input-bubble \
				   	input-cauchy \
					input-drop \
					input-eddy \
					input-obst \
					input-others \
					input-rt2d \
					input-rt3d \
					input-TG

EXTRA_DIST = $(input_directories)


AM_CPPFLAGS = -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_PROGRAMS = iFluid


iFluid_SOURCES = iFluid.cpp


iFluid_LDADD = $(builddir)/objs/libiFluid.la




tagsfile:
	ctags *.h *.cpp $(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/CGAL_Interface/*.h $(top_srcdir)/CGAL_Interface/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c


