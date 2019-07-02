AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -static
AM_LIBTOOLFLAGS = --preserve-dup-deps

noinst_HEADERS = \
../../../include/intfc/array.h \
../../../include/intfc/geom.h \
../../../include/intfc/iapi.h \
../../../include/intfc/iBase.h \
../../../include/intfc/iGeom.h \
../../../include/intfc/iloc.h \
../../../include/intfc/ilocprotos.h \
../../../include/intfc/iMesh.h \
../../../include/intfc/iMeshP.h \
../../../include/intfc/int.h \
../../../include/intfc/iprotos.h \
../../../include/intfc/iRel.h \
../../../include/intfc/iTaps.h \
../../../include/intfc/table.h \
../../../include/intfc/triangledefs.h \
../../../include/intfc/userint.h


noinst_LTLIBRARIES = libintfc.la

libintfc_la_SOURCES = $(noinst_HEADERS)

libintfc_la_SOURCES += \
../arrayutils.c \
../comp1d.c \
../comp2d.c \
../comp3d.c \
../comp.c \
../cross2d.c \
../geomutils.c \
../iblkb.c \
../iblkc2.c \
../iblkc3.c \
../icheck3d.c \
../idiagnostic.c \
../iecomps.c \
../ifourier.c \
../igrid.c \
../igview.c \
../imap.c \
../iMesh.c \
../iMeshP.c \
../imkcurve.c \
../imksurf.c \
../int3d.c \
../intfc.c \
../iprt3d.c \
../iredist.c \
../iredist_o2.c \
../irefl.c \
../iscatter.c \
../isect2d.c \
../isect3d.c \
../isub.c \
../isurgery.c \
../itopo.c \
../itrisset.c \
../iuserintfc.c \
../iwallsurf.c \
../ixgraph.c \
../map.c \
../ppcopy.c \
../setb1d.c \
../setb2d.c \
../setb3d.c \
../shift.c \
../top.c \
../triangle.c \
../trisurf.c \
../userhooks.c \
../zoom.c

#Add a dummy c++ source file so the linker behaves properly.
libintfc_la_SOURCES += ../dummyintfc.cpp

