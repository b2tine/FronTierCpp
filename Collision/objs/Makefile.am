AM_CPPFLAGS = -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


#NOTE: libCollision.la is a nonstatic libtool convenience library

#AM_LDFLAGS = -static
AM_LIBTOOLFLAGS = --preserve-dup-deps


noinst_LTLIBRARIES = ../../lib/libCollision.la



______lib_libCollision_la_SOURCES = \
../ifluid_state.h \
../collid.h \
../dcollid.cpp \
../dcollid3d.cpp


______lib_libCollision_la_LIBADD = @CGAL_LIBS@



EXTRA______lib_libCollision_la_DEPENDENCIES = $(top_builddir)/lib/libFronTier.la



$(EXTRA______lib_libCollision_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/src && $(MAKE) all

FORCE:



