AM_CPPFLAGS = -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @PETSC_LIBS@


noinst_LTLIBRARIES = ../../lib/libsolver.la


______lib_libsolver_la_SOURCES = \
../solver.h \
../cim_ellip.cpp \
../cim_parab.cpp \
../dual_ellip.cpp \
../ellip.cpp \
../hyperb.cpp \
../parab.cpp \
../solver.cpp


______lib_libsolver_la_LIBADD = $(top_builddir)/lib/libFronTier.la
								



$(______lib_libsolver_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/src && $(MAKE) all

FORCE:


clean-libtool:
