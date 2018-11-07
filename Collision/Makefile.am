input_directories = input

EXTRA_DIST = $(input_directories) \
			 filecheck.h \
			 filecheck.cpp \
			 vtk.cpp
			 

AM_CPPFLAGS = -I$(abs_top_srcdir)/Collision \
			  @CGAL_CFLAGS@ \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = test

				  
test_SOURCES = collid.h \
			   cdinit.cpp \
			   test.cpp

test_LDADD = $(abs_top_builddir)/lib/libCollision.la \
			 @CGAL_LIBS@ \
			 $(abs_top_builddir)/lib/libFronTier.la \
			 @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@



tagsfile:
	ctags *.h *.cpp ../Collision/*.h ../Collision/*.cpp ../iFluid/*.h ../iFluid/*.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c



