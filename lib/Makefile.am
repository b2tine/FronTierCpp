SUBDIRS = gmock

AM_CPPFLAGS = -I$(top_srcdir)/include \
			  -I$(top_srcdir)/include/front \
			  -I$(top_srcdir)/include/intfc \
			  -I$(top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


noinst_HEADERS = ../include/FronTier.h

noinst_LTLIBRARIES = libFronTier.la

libFronTier_la_SOURCES =
libFronTier_la_LIBADD = $(abs_top_builddir)/src/front/libfront.la \
  						$(abs_top_builddir)/src/intfc/libintfc.la \
  						$(abs_top_builddir)/src/util/libFTutil.la

if PETSC_OPT
noinst_LTLIBRARIES += libsolver.la libiFluid.la
endif


libsolver_la_CPPFLAGS = -I$(top_srcdir)/solver @PETSC_CFLAGS@ $(AM_CPPFLAGS)

libsolver_la_LIBADD = @PETSC_LIBS@ $(abs_top_builddir)/lib/libFronTier.la

libsolver_la_SOURCES = ../solver/solver.h \
					   ../solver/solver.cpp \
					   ../solver/ellip.cpp \
					   ../solver/dual_ellip.cpp \
					   ../solver/cim_ellip.cpp \
					   ../solver/hyperb.cpp \
					   ../solver/parab.cpp \
					   ../solver/cim_parab.cpp


libiFluid_la_CPPFLAGS = -I$(top_srcdir)/iFluid -I$(top_srcdir)/solver \
						@PETSC_CFLAGS@ $(AM_CPPFLAGS)

libiFluid_la_LIBADD = $(abs_top_builddir)/lib/libsolver.la

libiFluid_la_SOURCES = ../iFluid/iFluid.h \
					   ../iFluid/ifluid_state.h \
					   ../iFluid/iFbasic.cpp \
					   ../iFluid/iFcartsn2d.cpp \
					   ../iFluid/iFcartsn3d.cpp \
					   ../iFluid/iFsub.cpp \
					   ../iFluid/iFinit.cpp

