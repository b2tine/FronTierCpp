AM_CPPFLAGS = -I$(top_srcdir)/BVH \
			  -I$(top_srcdir)/airfoil \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_LTLIBRARIES = ../../lib/libBVH.la



______lib_libBVH_la_SOURCES = \
../HyperSurfElement.h \
../HyperSurfElement.cpp \
../BoundingVolume.h \
../BoundingVolume.cpp \
../BVH_Node.h \
../BVH_Node.cpp \
../BVH.h \
../BVH.cpp \
../BVH_Visualization.cpp \
../Query.cpp \
../CramersRule2d.h \
../BVH_util.cpp \
../BVH_util.h


______lib_libBVH_la_LIBADD = $(top_builddir)/lib/libairfoil.la



$(______lib_libBVH_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/airfoil/objs && $(MAKE) all


FORCE:


clean-libtool:
