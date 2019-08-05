input_directories = input-1d input-2d input-3d

EXTRA_DIST = $(input_directories)

AM_CPPFLAGS = -I$(top_srcdir)/crystal \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = crystal


crystal_SOURCES = crystal.h \
				  crysub.cpp \
				  ccartsn.cpp \
				  crystal.cpp


crystal_LDADD = $(top_builddir)/iFluid/objs/libiFluid.la




$(crystal_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/iFluid/objs libiFluid.la
	

FORCE:



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c


