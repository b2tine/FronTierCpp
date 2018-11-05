SUBDIRS = gmock

AM_CFLAGS = $(CFLAGS)
AM_CXXFLAGS = $(CXXFLAGS)

AM_CPPFLAGS = -I$(top_srcdir)/include \
			  -I$(top_srcdir)/include/front \
			  -I$(top_srcdir)/include/intfc \
			  -I$(top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


noinst_HEADERS = ../include/FronTier.h

noinst_LTLIBRARIES = libFronTier.la

libFronTier_la_SOURCES =
libFronTier_la_LIBADD = $(abs_top_builddir)/src/front/libfront.la \
  						$(abs_top_builddir)/src/intfc/libintfc.la \
  						$(abs_top_builddir)/src/util/libutil.la \
						@MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

if PETSC_OPT
noinst_LTLIBRARIES += libsolver.la libiFluid.la
endif


libsolver_la_CPPFLAGS = -I$(top_srcdir)/solver @PETSC_CFLAGS@ $(AM_CPPFLAGS)

libsolver_la_SOURCES = ../solver/solver.h \
					   ../solver/solver.cpp \
					   ../solver/ellip.cpp \
					   ../solver/dual_ellip.cpp \
					   ../solver/cim_ellip.cpp \
					   ../solver/hyperb.cpp \
					   ../solver/parab.cpp \
					   ../solver/cim_parab.cpp

libsolver_la_LIBADD = @PETSC_LIBS@ @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


libiFluid_la_CPPFLAGS = -I$(top_srcdir)/iFluid \
						-I$(top_srcdir)/solver \
						@PETSC_CFLAGS@ $(AM_CPPFLAGS)

libiFluid_la_SOURCES = ../iFluid/iFluid.h \
					   ../iFluid/iFbasic.cpp \
					   ../iFluid/iFcartsn2d.cpp \
					   ../iFluid/iFcartsn3d.cpp \
					   ../iFluid/iFinit.cpp \
					   ../iFluid/iFsub.cpp

libiFluid_la_LIBADD = @PETSC_LIBS@ @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


if CGAL_OPT
noinst_LTLIBRARIES += libCollision.la
if PETSC_OPT
noinst_LTLIBRARIES += libairfoil.la
endif
endif


libCollision_la_CPPFLAGS = -I$(top_srcdir)/Collision @CGAL_CFLAGS@ \
						   -I$(top_srcdir)/iFluid $(AM_CPPFLAGS)

libCollision_la_SOURCES = ../Collision/collid.h \
						  ../Collision/dcollid.cpp \
						  ../Collision/dcollid3d.cpp
						  
libCollision_la_LIBADD = @CGAL_LIBS@ 


libairfoil_la_CPPFLAGS = -I$(top_srcdir)/airfoil \
						 -I$(top_srcdir)/Collision @CGAL_CFLAGS@ \
						 -I$(top_srcdir)/iFluid \
						 -I$(top_srcdir)/solver @PETSC_CFLAGS@ \
						 $(AM_CPPFLAGS)

libairfoil_la_SOURCES = ../airfoil/airfoil.h \
						../airfoil/airfoil_sv.h \
						../airfoil/afcnpy.cpp \
						../airfoil/afdata.cpp \
						../airfoil/afinit3d.cpp \
						../airfoil/afinit.cpp \
						../airfoil/afprop.cpp \
						../airfoil/afsetd.cpp \
						../airfoil/aftest.cpp \
						../airfoil/afvelo.cpp \
						../airfoil/cgal.cpp

libairfoil_la_LIBADD = @CGAL_LIBS@ @PETSC_LIBS@ @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@



