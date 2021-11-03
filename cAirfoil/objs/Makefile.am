AM_CPPFLAGS = -I$(top_srcdir)/cAirfoil \
			  -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@




noinst_LTLIBRARIES = libcAirfoil.la



libcAirfoil_la_SOURCES = \
../fabric.h \
../fabric_sv.h \
../afprop.cpp \
../afcnpy.cpp \
../afdata.cpp \
../afinit3d.cpp \
../afinit.cpp \
../afsetd.cpp \
../aftest.cpp \
../afvelo.cpp \
../cgal.cpp \
../modules.cpp \
../fabric_api.cpp \
../cFabric.h \
../cFabric.cpp \
../cFbasic.cpp

libcAirfoil_la_SOURCES += \
../CD_HSE.h \
../CD_HSE.cpp \
../AABB.h \
../AABB.cpp \
../collid.h \
../dcollid.cpp \
../dcollid3d.cpp \
../dummy.cpp


libcAirfoil_la_LIBADD = $(top_builddir)/cFluid/objs/libcFluid.la


$(libcAirfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/cFluid/objs libcFluid.la


FORCE:



