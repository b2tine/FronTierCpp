input_directories = input-1d input-2d input-3d

EXTRA_DIST = $(input_directories)

AM_CPPFLAGS = -I$(abs_top_srcdir)/crystal \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @HDF_CFLAGS@ @GD_CFLAGS@

noinst_PROGRAMS = crystal

crystal_SOURCES = crystal.h \
				  crysub.cpp \
				  ccartsn.cpp \
				  crystal.cpp

crystal_LDADD = $(abs_top_builddir)/lib/libiFluid.la \
				$(abs_top_builddir)/lib/libsolver.la \
				@CGAL_LIBS@	@PETSC_LIBS@ @MPI_LIBS@ \
				$(abs_top_builddir)/lib/libFronTier.la \
				@HDF_LIBS@ @GD_LIBS@

tagsfile:
	ctags *.h *.cpp ../solver/*.h ../solver/*.cpp ../iFluid/*.h ../iFluid/*cpp ../include/*/*.h ../src/*/*.c


