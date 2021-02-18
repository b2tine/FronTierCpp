AM_CPPFLAGS = -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @PETSC_LIBS@




noinst_LTLIBRARIES = libsolver.la


libsolver_la_SOURCES = \
../solver.h \
../ellip.cpp \
../double_ellip.cpp \
../hyperb.cpp \
../parab.cpp \
../solver.cpp


if CGAL_INTFC_OPT
libsolver_la_LIBADD = $(top_builddir)/CGAL_Interface/objs/libCgalIntfc.la
else
libsolver_la_LIBADD = $(top_builddir)/lib/libFronTier.la
endif


if CGAL_INTFC_OPT
$(libsolver_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/CGAL_Interface/objs libCgalIntfc.la
else
$(libsolver_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all
endif

FORCE:


