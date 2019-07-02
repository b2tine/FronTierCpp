AM_CPPFLAGS = -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

LIB_SOLVER = $(top_builddir)/solver/objs/libsolver.la


noinst_LTLIBRARIES = libiFluid.la



libiFluid_la_SOURCES = \
../iFluid.h \
../ifluid_state.h \
../iFbasic.cpp \
../iFcartsn2d.cpp \
../iFcartsn3d.cpp \
../iFinit.cpp \
../iFsub.cpp
								

libiFluid_la_LIBADD = $(LIB_SOLVER)


$(libiFluid_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/solver/objs libsolver.la

FORCE:


