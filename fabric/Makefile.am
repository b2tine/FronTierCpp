SUBDIRS = objs
#SUBDIRS += tests


input_directories = input-2d \
					input-7ft \
					input-C9 \
					input-cloth \
					input-collsn \
					input-Cross \
					input-drum \
					input-G11 \
					input-G12D \
					input-intruder \
					input-others \
					input-T10 \
					input-young3d

EXTRA_DIST = $(input_directories) airfoil_gpu.cuh  airfoil_gpu.cu



AM_CPPFLAGS = -I$(top_srcdir)/fabric \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = fabric \
				  test1 \
				  fortrantest

				  
fabric_SOURCES = fabric.cpp
fabric_LDADD = $(builddir)/objs/libfabric.la


test1_SOURCES = test1.cpp
test1_LDADD = $(builddir)/objs/libfabric.la


fortrantest_SOURCES = fabric.F90
#fortrantest_LDADD =



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



