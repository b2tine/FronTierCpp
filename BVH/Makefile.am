SUBDIRS = objs tests

input_directories = input-OFF\
					input-collid-new

EXTRA_DIST = $(input_directories) \
			 CGAL_Point.h \
			 CGAL_Point.cpp
			 

AM_CPPFLAGS = -I$(top_srcdir)/BVH \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = OFF2FTsurf

				  
OFF2FTsurf_SOURCES = OFF2FTsurf.cpp

OFF2FTsurf_LDADD = $(top_builddir)/lib/libBVH.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



