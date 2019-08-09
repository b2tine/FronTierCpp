AM_CPPFLAGS = -I$(top_srcdir)/fabric \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_LTLIBRARIES = libfabric.la



libfabric_la_SOURCES = \
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
../collid.h \
../dcollid.cpp \
../dcollid3d.cpp \
../AABB.h \
../AABB.cpp \
../fluidinit.cpp \
../fluidprop.cpp \
../rigidbody.h \
../rigidbody.cpp \
../modules.cpp \
../fabric_api_internal.cpp \
../fabric_api.cpp



libfabric_la_LIBADD = $(top_builddir)/lib/libFronTier.la



$(libfabric_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:



