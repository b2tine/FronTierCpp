AM_CPPFLAGS = -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_LTLIBRARIES = libiFluid.la



libiFluid_la_SOURCES = \
../iFluid.h \
../ifluid_state.h \
../iFbasic.cpp \
../iFcartsn2d.cpp \
../iFcartsn3d.cpp \
../iFinit.cpp \
../iFsub.cpp \
../iFturb.h \
../iFturb.cpp \
../keps.h \
../keps.cpp \
../rigidbody.h \
../rigidbody.cpp \
../iFvisualize.cpp
								

libiFluid_la_LIBADD = $(top_builddir)/solver/objs/libsolver.la


$(libiFluid_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/solver/objs libsolver.la

FORCE:


