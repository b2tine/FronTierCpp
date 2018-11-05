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
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = iFluid

iFluid_SOURCES = iFluid.h \
				 ifluid_state.h \
				 iFluid.cpp

iFluid_LDADD = $(abs_top_builddir)/lib/libiFluid.la \
			   $(abs_top_builddir)/lib/libsolver.la \
			   @PETSC_LIBS@ @MPI_LIBS@ \
			   $(abs_top_builddir)/lib/libFronTier.la \
			   @HDF_LIBS@ @GD_LIBS@
			   


tagsfile:
	ctags *.h *.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c


clean-local:
	-rm -rf tags
	-rm -rf TAGS
	-rm -rf CTAGS
	-rm -rf Makefile

