SUBDIRS = objs tests


EXTRA_DIST = input


AM_CPPFLAGS = -I$(top_srcdir)/folding \
			  @NLOPT_CFLAGS@ \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/Collision \
			  -I$(top_srcdir)/iFluid \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CVODE_LIBS@ \
			 @NLOPT_LIBS@ \
			 @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@



noinst_PROGRAMS = test

				  
test_SOURCES = test.cpp


test_LDADD = $(builddir)/objs/libfolding.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



