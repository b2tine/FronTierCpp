AM_CPPFLAGS = -I$(top_srcdir)/cfabric \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_LTLIBRARIES = libcfabric.la



libcfabric_la_SOURCES = \
../fabric.h \
../fabric_sv.h \
../state.h \
../state.cpp \
../afprop.cpp \
../afcnpy.cpp \
../afdata.cpp \
../afinit3d.cpp \
../afinit.cpp \
../afsetd.cpp \
../aftest.cpp \
../afvelo.cpp \
../cgal_intfc.h \
../cgal_intfc.cpp \
../cgal.cpp \
../collid.h \
../dcollid.cpp \
../dcollid3d.cpp \
../CD_HSE.h \
../CD_HSE.cpp \
../AABB.h \
../AABB.cpp \
../cFluid.h \
../cFcartsn.cpp \
../cFphys.cpp \
../cFbasic.cpp \
../cFbase.cpp \
../cFriem.cpp \
../cFweno.cpp \
../cFtvd.cpp \
../cFinit.cpp \
../cFsub.cpp \
../cFeos.cpp \
../fluid.h \
../fluidinit.cpp \
../iFprop.cpp \
../rigidbody.h \
../rigidbody.cpp \
../modules.cpp \
../fabric_api.cpp



libcfabric_la_LIBADD = $(top_builddir)/lib/libFronTier.la



$(libcfabric_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:



