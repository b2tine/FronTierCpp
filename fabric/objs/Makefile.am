AM_CPPFLAGS = -I$(top_srcdir)/fabric \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_LTLIBRARIES = libfabric.la



libfabric_la_SOURCES = \
../fabric.h \
../fabric_sv.h \
../state.h \
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
../CD_HSE.h \
../CD_HSE.cpp \
../AABB.h \
../AABB.cpp \
../fluid.h \
../fluidinit.cpp \
../fluidprop.cpp \
../rigidbody.h \
../rigidbody.cpp \
../modules.cpp \
../fabric_api.cpp



if CGAL_INTFC_OPT
libfabric_la_LIBADD = $(top_builddir)/CGAL_Interface/objs/libCgalIntfc.la
else
libfabric_la_LIBADD = $(top_builddir)/lib/libFronTier.la
endif



if CGAL_INTFC_OPT
$(libfabric_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/CGAL_Interface/objs libCgalIntfc.la
else
$(libfabric_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all
endif


FORCE:



