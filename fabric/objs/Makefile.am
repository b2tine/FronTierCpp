AM_CPPFLAGS = -I$(top_srcdir)/fabric \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_CPPFLAGS += -DDAMPING_FORCE



noinst_LTLIBRARIES = libairfoil.la



libairfoil_la_SOURCES = \
../airfoil.h \
../airfoil_sv.h \
../afprop.cpp \
../afcnpy.cpp \
../afdata.cpp \
../afinit3d.cpp \
../afinit.cpp \
../afsetd.cpp \
../aftest.cpp \
../afvelo.cpp \
../cgal.cpp \
../fluidprop.cpp



libairfoil_la_LIBADD = $(top_builddir)/Collision/objs/libCollision.la



$(libairfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/Collision/objs libCollision.la

FORCE:

