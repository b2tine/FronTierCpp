SUBDIRS = objs

input_directories = input

EXTRA_DIST = $(input_directories) \
			 filecheck.h \
			 filecheck.cpp \
			 vtk.cpp
			 

AM_CPPFLAGS = -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = test

				  
test_SOURCES = test.cpp

test_LDADD = $(top_builddir)/lib/libCollision.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/iFluid/ifluid_state.h \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



