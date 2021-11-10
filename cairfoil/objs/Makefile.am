AM_CPPFLAGS = -I$(top_srcdir)/cairfoil \
			  -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@




noinst_LTLIBRARIES = libcairfoil.la



libcairfoil_la_SOURCES = \
../FabricManager.h \
../FabricManager.cpp \
../airfoil.h \
../airfoil_sv.h \
../afcnpy.cpp \
../afdata.cpp \
../afinit3d.cpp \
../afinit.cpp \
../afprop.cpp \
../afsetd.cpp \
../aftest.cpp \
../afstat.cpp \
../afvelo.cpp \
../bending.h \
../bending.cpp \
../cgal.cpp \
../modules.cpp \
../cFabric.h \
../cFabric.cpp

libcairfoil_la_SOURCES += \
../CD_HSE.h \
../CD_HSE.cpp \
../AABB.h \
../AABB.cpp \
../collid.h \
../dcollid.cpp \
../dcollid3d.cpp \
../vtk.h \
../vtk.cpp


libcairfoil_la_LIBADD = $(top_builddir)/cFluid/objs/libcFluid.la


$(libcairfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/cFluid/objs libcFluid.la


FORCE:



