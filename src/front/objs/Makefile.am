AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = -static
AM_LIBTOOLFLAGS = --preserve-dup-deps


noinst_HEADERS = \
../../../include/front/fapi.h \
../../../include/front/fcell.h \
../../../include/front/fdecs.h \
../../../include/front/fpatrecon.h \
../../../include/front/fprotos.h \
../../../include/front/frp.h \
../../../include/front/fuserint.h \
../../../include/front/fvelo.h


noinst_LTLIBRARIES = libfront.la

libfront_la_SOURCES = $(noinst_HEADERS)

libfront_la_SOURCES += \
../fadv3d.c \
../fadv.c \
../fbdry1.c \
../fbdry2.c \
../fbdry3.c \
../fbdry4.c \
../fcell.c \
../fcheck3d.c \
../fcollsn3d.c \
../fcorrspnd.c \
../fcrosscur.c \
../fcrossext.c \
../fcrstatus.c \
../fdiagnostic.c \
../fgb2d.c \
../fgb3comp.c \
../fgb3d.c \
../fgb3dutil.c \
../fglobal.c \
../fgrid.c \
../finit.c \
../fint.c \
../fmap3.c \
../fmap2.c \
../fmap.c \
../fnode.c \
../fnodesub.c \
../fpatch3d.c \
../fprint.c \
../fprop2d.c \
../fprop3d.c \
../fredist1d.c \
../fredist2d.c \
../fredist3d.c \
../fredist.c \
../frp1.c \
../frp2.c \
../fscat1d.c \
../fscat2d.c \
../fscat3d1.c \
../fscat3d2.c \
../fscat3d3.c \
../fscatter.c \
../fstate2d.c \
../fstate.c \
../fsub.c \
../ftop.c \
../funtan2d.c \
../funtan3d.c \
../fuserhooks.c \
../fuserintfc.c \
../fvelo.c

#Add a dummy c++ source file so the linker behaves properly.
libfront_la_SOURCES += ../dummyfront.cpp


