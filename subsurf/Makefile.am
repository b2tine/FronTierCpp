EXTRA_DIST = input-2d

AM_CPPFLAGS = -I$(top_srcdir)/subsurf \
			  -I$(top_srcdir)/crystal \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = subsurf


subsurf_SOURCES = subsurf.h \
				  subsub.cpp \
				  subsurf.cpp

subsurf_LDADD =	$(top_builddir)/crystal/ccartsn.o \
				$(top_builddir)/crystal/crysub.o \
				$(top_builddir)/iFluid/objs/libiFluid.la


$(subsurf_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/crystal all
	@$(MAKE) -C $(top_builddir)/iFluid/objs libiFluid.la

FORCE:



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/crystal/*.h $(top_srcdir)/crystal/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

