EXTRA_DIST = input-2d

AM_CPPFLAGS = -I$(abs_top_srcdir)/subsurf \
			  -I$(abs_top_srcdir)/crystal \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = subsurf


subsurf_SOURCES = subsurf.h \
				  subsurf.cpp \
				  subsub.cpp

subsurf_LDADD =	$(abs_top_builddir)/crystal/ccartsn.o \
				$(abs_top_builddir)/crystal/crysub.o \
				$(abs_top_builddir)/lib/libiFluid.la \
				$(abs_top_builddir)/lib/libsolver.la \
				@CGAL_LIBS@	@PETSC_LIBS@ @MPI_LIBS@ \
				$(abs_top_builddir)/lib/libFronTier.la \
				@HDF_LIBS@ @GD_LIBS@

tagsfile:
	ctags *.h *.cpp ../crystal/*.h ../crystal/*.cpp ../iFluid/*.h ../iFluid/*.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c


