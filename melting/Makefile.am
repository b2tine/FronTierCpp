input_directories = input

EXTRA_DIST = $(input_directories) mrun1d mrun2d mrun3d


AM_CPPFLAGS = -I$(top_srcdir)/melting \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = melting


melting_SOURCES = melting.h \
 				  msub.cpp \
 				  minit.cpp \
 				  mcartsn.cpp \
 				  mprop.cpp \
 				  melting.cpp


melting_LDADD = $(top_builddir)/iFluid/objs/libiFluid.la



$(melting_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/iFluid/objs libiFluid.la

FORCE:



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c


