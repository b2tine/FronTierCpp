AM_CPPFLAGS = -I$(top_srcdir)/AngStiff_sprModel \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@




noinst_LTLIBRARIES = libangstiff.la



libangstiff_la_SOURCES = \
../airfoil.h \
../airfoil_sv.h \
../afcnpy.cpp \
../afdata.cpp \
../afprop.cpp \
../afsetd.cpp \
../aftest.cpp \
../afvelo.cpp \
../modules.cpp



libangstiff_la_LIBADD = $(top_builddir)/Collision/objs/libCollision.la



$(libangstiff_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/Collision/objs libCollision.la

FORCE:

