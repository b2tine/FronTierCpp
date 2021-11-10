AM_CPPFLAGS = -I$(top_srcdir)/cairfoil \
			  -I$(top_srcdir)/folding \
			  @NLOPT_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/CGAL_Interface \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@ \
			 @CVODE_LIBS@


if FOLDING_OPT
AM_LDFLAGS += @NLOPT_LIBS@
endif



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



if FOLDING_OPT
libcairfoil_la_LIBADD = $(top_builddir)/folding/objs/libfolding.la
else
libcairfoil_la_LIBADD = $(top_builddir)/Collision/objs/libCollision.la
endif


if FOLDING_OPT
$(libcairfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/folding/objs libfolding.la
else
$(libcairfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/Collision/objs libCollision.la
endif


FORCE:



