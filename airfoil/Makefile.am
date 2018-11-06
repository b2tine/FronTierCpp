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

EXTRA_DIST = $(input_directories) airfoil_gpu.cuh airfoil_gpu.cu

AM_CPPFLAGS = -I$(abs_top_srcdir)/airfoil \
			  -I$(abs_top_srcdir)/Collision \
			  @CGAL_CFLAGS@ \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_CPPFLAGS += -DDAMPING_FORCE -DCOLLISION_DETECTION


noinst_PROGRAMS = airfoil

				  
airfoil_SOURCES = airfoil.h \
				  airfoil_sv.h \
				  airfoil.cpp

airfoil_LDADD =	$(abs_top_builddir)/lib/libairfoil.la \
				$(abs_top_builddir)/lib/libCollision.la \
				$(abs_top_builddir)/lib/libiFluid.la \
				$(abs_top_builddir)/lib/libsolver.la \
				@CGAL_LIBS@	@PETSC_LIBS@ @MPI_LIBS@ \
				$(abs_top_builddir)/lib/libFronTier.la \
				@HDF_LIBS@ @GD_LIBS@



tagsfile:
	ctags *.h *.cpp ../Collision/*.h ../Collision/*.cpp ../iFluid/*.h ../iFluid/*.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c



