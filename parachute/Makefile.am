input_directories = input-C9 \
					input-cross \
					input-G11 \
					input-intruder \
					input-mchutes \
					input-others \
					input-T10 \
					input-T11

EXTRA_DIST = $(input_directories)


AM_CPPFLAGS = -I$(top_srcdir)/parachute \
			  -I$(top_srcdir)/airfoil \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@



noinst_PROGRAMS = parachute



parachute_SOURCES = modules.cpp \
					parachute.cpp


parachute_LDADD = $(top_builddir)/airfoil/objs/libairfoil.la



$(parachute_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/airfoil/objs libairfoil.la

FORCE:



tagsfile:
	ctags parachute.cpp modules.cpp \
		$(top_srcdir)/airfoil/*.h $(top_srcdir)/airfoil/*.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/iFcartsn3d.cpp \
		$(top_srcdir)/iFluid/iFinit.cpp $(top_srcdir)/iFluid/iFbasic.cpp \
		$(top_srcdir)/iFluid/iFsub.cpp $(top_srcdir)/iFluid/keps.cpp \
		$(top_srcdir)/iFluid/rigidbody.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/CGAL_Interface/*.h $(top_srcdir)/CGAL_Interface/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c


