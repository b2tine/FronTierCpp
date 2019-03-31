AM_CPPFLAGS = -I$(top_srcdir)/BVH \
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
../Query.cpp


______lib_libBVH_la_LIBADD = $(top_builddir)/lib/libiFluid.la



$(______lib_libBVH_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/iFluid/objs && $(MAKE) all


FORCE:


clean-libtool:
