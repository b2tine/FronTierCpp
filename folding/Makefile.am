SUBDIRS = objs


EXTRA_DIST = input


#TODO: update nlopt inc/lib
AM_CPPFLAGS = -I$(top_srcdir)/folding \
			  -I$(top_srcdir)/folding/nlopt \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = -L$(top_srcdir)/folding/nlopt -lnlopt

noinst_PROGRAMS = test

				  
test_SOURCES = test.cpp


test_LDADD = $(top_builddir)/lib/libfolding.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



