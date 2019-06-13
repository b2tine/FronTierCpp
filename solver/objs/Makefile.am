AM_CPPFLAGS = -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @PETSC_LIBS@


LIB_FT = $(top_builddir)/lib/libFronTier.la


noinst_LTLIBRARIES = libsolver.la


libsolver_la_SOURCES = \
../solver.h \
../cim_ellip.cpp \
../cim_parab.cpp \
../dual_ellip.cpp \
../ellip.cpp \
../hyperb.cpp \
../parab.cpp \
../solver.cpp


libsolver_la_LIBADD = $(LIB_FT)
								

$(LIB_FT): libraries

libraries:
	@$(MAKE) -C $(top_builddir)/lib libFronTier.la


