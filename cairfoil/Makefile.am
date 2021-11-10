SUBDIRS = objs


AM_CPPFLAGS = -I$(top_srcdir)/cairfoil \
			  -I$(top_srcdir)/folding \
			  @NLOPT_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/CGAL_Interface \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@ \
			 @CVODE_LIBS@


if FOLDING_OPT
AM_LDFLAGS += @NLOPT_LIBS@
endif


noinst_PROGRAMS = cairfoil

cairfoil_SOURCES = cairfoil.cpp

cairfoil_LDADD = $(builddir)/objs/libcairfoil.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/folding/*.h $(top_srcdir)/folding/*.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/cFluid/*.h $(top_srcdir)/cFluid/*.cpp \
		$(top_srcdir)/CGAL_Interface/*.h $(top_srcdir)/CGAL_Interface/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



