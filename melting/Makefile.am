input_directories = input

EXTRA_DIST = $(input_directories) mrun1d mrun2d mrun3d


AM_CPPFLAGS = -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -all-static


noinst_PROGRAMS = melting

melting_SOURCES = iFluid.h \
 				  melting.h \
 				  msub.cpp \
 				  minit.cpp \
 				  mcartsn.cpp \
 				  mprop.cpp \
 				  melting.cpp

melting_LDADD = $(abs_top_builddir)/lib/libiFluid.la \
 				$(abs_top_builddir)/lib/libsolver.la \
 				@PETSC_LIBS@ @MPI_LIBS@ \
 				$(abs_top_builddir)/lib/libFronTier.la \
 				@HDF_LIBS@ @GD_LIBS@

EXTRA_melting_DEPENDENCIES = $(abs_top_builddir)/lib/libiFluid.la \
 							 $(abs_top_builddir)/lib/libsolver.la \
 							 $(abs_top_builddir)/lib/libFronTier.la


tagsfile:
	ctags *.h *.cpp ../iFluid/*.h ../iFluid/*.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c


