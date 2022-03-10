AM_CPPFLAGS = -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@


noinst_LTLIBRARIES = libCgalIntfc.la



libCgalIntfc_la_SOURCES = \
../cgal_intfc.h \
../cgal_intfc.cpp \
../cgal_surf.h



libCgalIntfc_la_LIBADD = $(top_builddir)/lib/libFronTier.la



$(libCgalIntfc_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:



