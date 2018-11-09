AM_CPPFLAGS = -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_LTLIBRARIES = ../../lib/libCollision.la



______lib_libCollision_la_SOURCES = \
../collid.h \
../dcollid.cpp \
../dcollid3d.cpp \
../cdinit.cpp \
../AABB.cpp \
../AABB.h


______lib_libCollision_la_LIBADD = @CGAL_LIBS@ @ARMA_LIBS@ \
								   $(top_builddir)/lib/libiFluid.la



$(______lib_libCollision_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/iFluid/objs && $(MAKE) all

FORCE:



