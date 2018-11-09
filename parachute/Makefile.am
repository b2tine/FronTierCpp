input_directories = input-C9 \
					input-cross \
					input-G11 \
					input-intruder \
					input-mchutes \
					input-others \
					input-T10 \
					input-T11

EXTRA_DIST = $(input_directories)


AM_CPPFLAGS = -I$(top_srcdir)/parachute \
			  -I$(top_srcdir)/airfoil \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = parachute



parachute_SOURCES = modules.cpp \
					parachute.cpp


parachute_LDADD = $(top_builddir)/lib/libairfoil.la



$(parachute_DEPENDENCIES): FORCE
	@cd $(top_builddir)/airfoil/objs && $(MAKE) all

FORCE:



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c


