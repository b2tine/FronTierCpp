SUBDIRS = objs


AM_CPPFLAGS = -I$(top_srcdir)/cAirfoil \
			  -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@



noinst_PROGRAMS = cfabric_Cplus


cfabric_Cplus_SOURCES = cfabric_Cplus.cpp
cfabric_Cplus_LDADD = $(builddir)/objs/libcAirfoil.la



tagsfile:
	ctags \
		$(srcdir)/*.h $(srcdir)/*.cpp \
		$(top_srcdir)/cFluid/*.h $(top_srcdir)/cFluid/*.cpp \
		$(top_srcdir)/CGAL_Interface/*.h $(top_srcdir)/CGAL_Interface/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

