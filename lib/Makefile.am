SUBDIRS = gmock

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
noinst_LTLIBRARIES += libsolver.la
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


if CGAL_OPT
noinst_LTLIBRARIES += libCollision.la
endif

libCollision_la_CPPFLAGS = -I$(top_srcdir)/Collision @CGAL_CFLAGS@ \
						   -I$(top_srcdir)/iFluid $(AM_CPPFLAGS)

libCollision_la_SOURCES = ../Collision/collid.h \
						  ../Collision/dcollid.cpp \
						  ../Collision/dcollid3d.cpp
						  
libCollision_la_LIBADD = @CGAL_LIBS@ 

