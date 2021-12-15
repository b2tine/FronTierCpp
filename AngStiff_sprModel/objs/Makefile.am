AM_CPPFLAGS = -I$(top_srcdir)/AngStiff_sprModel \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@



noinst_LTLIBRARIES = libangstiff.la



libangstiff_la_SOURCES = \
../state.h \
../airfoil.h \
../airfoil_sv.h \
../afcnpy.cpp \
../afdata.cpp \
../afprop.cpp \
../afsetd.cpp \
../aftest.cpp \
../afvelo.cpp \
../modules.cpp



libangstiff_la_LIBADD = $(top_builddir)/lib/libFronTier.la



$(libangstiff_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:

