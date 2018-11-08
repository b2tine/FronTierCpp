AM_CPPFLAGS = -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -static
AM_LIBTOOLFLAGS = --preserve-dup-deps


noinst_LTLIBRARIES = ../../lib/libiFluid.la



______lib_libiFluid_la_SOURCES = \
../iFluid.h \
../ifluid_state.h \
../iFbasic.cpp \
../iFcartsn2d.cpp \
../iFcartsn3d.cpp \
../iFinit.cpp \
../iFsub.cpp
								

______lib_libiFluid_la_LIBADD = $(top_builddir)/lib/libsolver.la



$(______lib_libiFluid_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/solver && $(MAKE) all

FORCE:



