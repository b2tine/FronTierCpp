AM_CPPFLAGS = -I$(top_srcdir)/airfoil \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_CPPFLAGS += -DDAMPING_FORCE -DCOLLISION_DETECTION



noinst_LTLIBRARIES = ../../lib/libairfoil.la



______lib_libairfoil_la_SOURCES = \
../airfoil.h \
../airfoil_sv.h \
../afcnpy.cpp \
../afdata.cpp \
../afinit3d.cpp \
../afinit.cpp \
../afprop.cpp \
../afsetd.cpp \
../aftest.cpp \
../afvelo.cpp \
../cgal.cpp



______lib_libairfoil_la_LIBADD = $(top_builddir)/lib/libCollision.la




$(______lib_libairfoil_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/Collision/objs && $(MAKE) all

FORCE:



clean-libtool:
