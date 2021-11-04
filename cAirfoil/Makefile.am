SUBDIRS = objs


AM_CPPFLAGS = -I$(top_srcdir)/cAirfoil \
			  -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@



noinst_PROGRAMS = cairfoil \
				  cFluid_test



cairfoil_SOURCES = fabric_Cplus.cpp

cairfoil_LDADD = $(builddir)/objs/libcAirfoil.la



cFluid_test_SOURCES = cFluid_test.cpp

cFluid_test_LDADD = $(builddir)/objs/libcAirfoil.la



tagsfile:
	ctags \
		$(srcdir)/*.h $(srcdir)/*.cpp \
		$(top_srcdir)/cFluid/*.h $(top_srcdir)/cFluid/*.cpp \
		$(top_srcdir)/CGAL_Interface/*.h $(top_srcdir)/CGAL_Interface/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

