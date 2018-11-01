AM_CPPFLAGS = -I$(abs_top_srcdir)/subsurf \
			  -I$(abs_top_srcdir)/crystal \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver @PETSC_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = subsurf

subsurf_LDADD =	$(abs_top_builddir)/crystal/ccartsn.o \
				$(abs_top_builddir)/crystal/crysub.o \
				-L$(abs_top_builddir)/lib -liFluid \
				@PETSC_LIBS@ @MPI_LIBS@	@HDF_LIBS@ @GD_LIBS@

EXTRA_DIST = input-2d

subsurf_SOURCES = subsurf.h \
				  subsurf.cpp \
				  subsub.cpp


tagsfile:
	ctags *.h *.cpp ../crystal/*.h ../crystal/*.cpp ../iFluid/*.h ../iFluid/*.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c
