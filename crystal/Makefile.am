AM_CPPFLAGS = -I$(abs_top_srcdir)/crystal \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

noinst_PROGRAMS = crystal

crystal_LDADD = -L$(abs_top_builddir)/lib -lsolver \
			   @PETSC_LIBS@ @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

crystal_SOURCES = crystal.h \
				  crystal_state.h \
				  ../iFluid/ifluid_state.h \
				  crysub.cpp \
				  ccartsn.cpp \
				  crystal.cpp

tagsfile:
	ctags *.h *.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c
