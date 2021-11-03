AM_CPPFLAGS = -I$(top_srcdir)/cAirfoil \
			  @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@




noinst_LTLIBRARIES = libcAirfoil.la



libcAirfoil_la_SOURCES = \
../dummy.cpp


libcAirfoil_la_LIBADD = $(top_builddir)/cFluid/objs/libcFluid.la


$(libcAirfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/cFluid/objs libcFluid.la


FORCE:



