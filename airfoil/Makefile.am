SUBDIRS = objs tests


input_directories = input-2d \
					input-7ft \
					input-C9 \
					input-cloth \
					input-collsn \
					input-Cross \
					input-drum \
					input-G11 \
					input-G12D \
					input-intruder \
					input-others \
					input-T10 \
					input-young3d

EXTRA_DIST = $(input_directories)



AM_CPPFLAGS = -I$(top_srcdir)/airfoil \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = airfoil

				  
airfoil_SOURCES = airfoil.cpp


airfoil_LDADD =	$(builddir)/objs/libairfoil.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



