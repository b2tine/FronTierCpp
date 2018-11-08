SUBDIRS = ../src


AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -static
AM_LIBTOOLFLAGS = --preserve-dup-deps



noinst_LTLIBRARIES = 


if CGAL_OPT
noinst_LTLIBRARIES += libCollision.la
if PETSC_OPT
noinst_LTLIBRARIES += libairfoil.la
endif
endif

SUBDIRS += gmock

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



