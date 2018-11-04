AM_CPPFLAGS = -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

noinst_PROGRAMS = iFluid

iFluid_DEPENDENCIES = $(abs_top_builddir)/lib/libiFluid.la \
					  $(abs_top_builddir)/lib/libFronTier.la

noinst_HEADERS = iFluid.h \
				 ../solver/solver.h \
				 ../include/FronTier.h \
				 ifluid_state.h

iFluid_SOURCES = iFluid.cpp

iFluid_LDADD = -L$(abs_top_builddir)/lib -liFluid \
			   @PETSC_LIBS@ @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


tagsfile:
	ctags *.h *.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c


input_directories = input-bubble \
				   	input-cauchy \
					input-drop \
					input-eddy \
					input-obst \
					input-others \
					input-rt2d \
					input-rt3d \
					input-TG

EXTRA_DIST = $(input_directories)
