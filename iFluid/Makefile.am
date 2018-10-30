AM_CPPFLAGS = -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

noinst_PROGRAMS = iFluid

iFluid_LDADD = -L$(abs_top_builddir)/lib -liFluid \
			   @PETSC_LIBS@ @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

iFluid_SOURCES = iFluid.cpp


tagsfile:
	ctags *.h *.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c
