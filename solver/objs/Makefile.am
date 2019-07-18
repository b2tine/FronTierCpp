AM_CPPFLAGS = -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @PETSC_LIBS@


LIB_FT = $(top_builddir)/lib/libFronTier.la


noinst_LTLIBRARIES = libsolver.la


libsolver_la_SOURCES = \
../solver.h \
../ellip.cpp \
../double_ellip.cpp \
../hyperb.cpp \
../parab.cpp \
../solver.cpp


libsolver_la_LIBADD = $(LIB_FT)
								

$(libsolver_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:


