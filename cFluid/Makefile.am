AM_CPPFLAGS = -I$(abs_top_srcdir)/cFluid \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


noinst_PROGRAMS = cFluid

cFluid_LDADD = -L$(abs_top_builddir)/lib -lcFluid \
			   @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

cFluid_SOURCES = cFluid.cpp


tagsfile:
	ctags *.h *.cpp ../include/*/*.h ../src/*/*.c


input_directories = input-1d input-2d input-3d

EXTRA_DIST = $(input_directories)
