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

EXTRA_DIST = $(input_directories) airfoil_gpu.cuh airfoil_gpu.cu



AM_CPPFLAGS = -I$(top_srcdir)/airfoil \
			  -I$(top_srcdir)/folding \
			  @NLOPT_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@ \
			 @PETSC_LIBS@ \
			 @CVODE_LIBS@

if FOLDING_OPT
AM_LDFLAGS += @NLOPT_LIBS@
endif


noinst_PROGRAMS = airfoil

				  
airfoil_SOURCES = airfoil.cpp


airfoil_LDADD =	$(builddir)/objs/libairfoil.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/folding/*.h $(top_srcdir)/folding/*.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/iFcartsn3d.cpp \
		$(top_srcdir)/iFluid/iFinit.cpp $(top_srcdir)/iFluid/iFbasic.cpp \
		$(top_srcdir)/iFluid/iFsub.cpp $(top_srcdir)/iFluid/rigidbody.cpp \
		$(top_srcdir)/iFluid/iFturb.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/CGAL_Interface/*.h $(top_srcdir)/CGAL_Interface/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



