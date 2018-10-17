AM_CPPFLAGS = -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  -I$(abs_top_srcdir)/solver \
			  @HDF_CFLAGS@ @PETSC_CFLAGS@ @MPI_CFLAGS@

AM_LDFLAGS = -L$(abs_top_builddir)/lib -lfront -lintfc -lutil \
			 @HDF_LIBS@ @PETSC_LIBS@ @MPI_LIBS@

noinst_PROGRAMS = iFluid

iFluid_SOURCES = iFluid.h \
				 ifluid_state.h \
				 iFluid.cpp \
				 iFbasic.cpp \
				 iFcartsn2d.cpp \
				 iFcartsn3d.cpp \
				 iFsub.cpp \
				 iFinit.cpp \
				 ../solver/solver.cpp \
				 ../solver/ellip.cpp \
				 ../solver/cim_ellip.cpp \
				 ../solver/dual_ellip.cpp \
				 ../solver/hyperb.cpp \
				 ../solver/parab.cpp

