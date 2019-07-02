AM_CPPFLAGS = -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_LTLIBRARIES = libCollision.la



libCollision_la_SOURCES = \
../collid.h \
../dcollid.cpp \
../dcollid3d.cpp \
../cdinit.cpp \
../AABB.cpp \
../AABB.h


libCollision_la_LIBADD = $(top_builddir)/iFluid/objs/libiFluid.la



$(libCollision_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/iFluid/objs libiFluid.la

FORCE:

