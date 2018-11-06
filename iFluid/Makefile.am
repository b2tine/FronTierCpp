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


AM_CPPFLAGS = -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -all-static


noinst_PROGRAMS = iFluid

iFluid_SOURCES = iFluid.h \
				 ifluid_state.h \
				 iFluid.cpp

iFluid_LDADD = $(abs_top_builddir)/lib/libiFluid.la \
			   $(abs_top_builddir)/lib/libsolver.la \
			   @PETSC_LIBS@ @MPI_LIBS@ \
			   $(abs_top_builddir)/lib/libFronTier.la \
			   @HDF_LIBS@ @GD_LIBS@

EXTRA_iFluid_DEPENDENCIES = $(abs_top_builddir)/lib/libiFluid.la \
							$(abs_top_builddir)/lib/libsolver.la \
			 				$(abs_top_builddir)/lib/libFronTier.la


tagsfile:
	ctags *.h *.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c


