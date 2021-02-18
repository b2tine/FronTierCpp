SUBDIRS = objs tests

input_directories = input

EXTRA_DIST = $(input_directories) \
			 filecheck.h \
			 filecheck.cpp \
			 vtk.cpp
			 

AM_CPPFLAGS = -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_PROGRAMS = test

				  
test_SOURCES = test.cpp \
			   cdinit.cpp

test_LDADD = $(builddir)/objs/libCollision.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



