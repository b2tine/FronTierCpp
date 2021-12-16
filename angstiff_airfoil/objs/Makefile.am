AM_CPPFLAGS = -I$(top_srcdir)/angstiff_airfoil \
			  -I$(top_srcdir)/folding \
			  @NLOPT_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@ \
			 @PETSC_LIBS@ \
			 @CVODE_LIBS@


if FOLDING_OPT
AM_LDFLAGS += @NLOPT_LIBS@
endif



noinst_LTLIBRARIES = libangstiffairfoil.la



libangstiffairfoil_la_SOURCES = \
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
../cgal.cpp


if FOLDING_OPT
libangstiffairfoil_la_LIBADD = $(top_builddir)/folding/objs/libfolding.la
else
libangstiffairfoil_la_LIBADD = $(top_builddir)/Collision/objs/libCollision.la
endif


if FOLDING_OPT
$(libangstiffairfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/folding/objs libfolding.la
else
$(libangstiffairfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/Collision/objs libCollision.la
endif


FORCE:



