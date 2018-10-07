SUBDIRS = gmock

AM_CPPFLAGS = -I$(top_srcdir)/include/front \
			  -I$(top_srcdir)/include/intfc \
			  -I$(top_srcdir)/include/util @HDF_CFLAGS@

AM_LDFLAGS = @HDF_LIBS@

noinst_LTLIBRARIES = libutil.la \
					 libintfc.la \
					 libfront.la

noinst_HEADERS = ../include/util/cdecs.h \
				 ../include/util/fnamedebug.h \
				 ../include/util/navdecs.h \
				 ../include/util/plotdecs.h \
				 ../include/util/uapi.h \
				 ../include/util/uprotos.h \
				 ../include/util/vmalloc.h \
				 ../include/intfc/array.h \
				 ../include/intfc/geom.h \
				 ../include/intfc/iapi.h \
				 ../include/intfc/iBase.h \
				 ../include/intfc/iGeom.h \
				 ../include/intfc/iloc.h \
				 ../include/intfc/ilocprotos.h \
				 ../include/intfc/iMesh.h \
				 ../include/intfc/iMeshP.h \
				 ../include/intfc/int.h \
				 ../include/intfc/iprotos.h \
				 ../include/intfc/iRel.h \
				 ../include/intfc/iTaps.h \
				 ../include/intfc/table.h \
				 ../include/intfc/triangledefs.h \
				 ../include/intfc/userint.h \
				 ../include/front/fapi.h \
				 ../include/front/fdecs.h \
				 ../include/front/fpatrecon.h \
				 ../include/front/fprotos.h \
				 ../include/front/frp.h \
				 ../include/front/fuserint.h \
				 ../include/front/fvelo.h 

libutil_la_SOURCES = ../src/util/cleanup.c \
					 ../src/util/debug.c \
					 ../src/util/error.c \
					 ../src/util/fft.c \
					 ../src/util/fgetstrin.c \
					 ../src/util/fsort.c \
					 ../src/util/machine.c \
					 ../src/util/matrix.c \
					 ../src/util/other.c \
					 ../src/util/output.c \
					 ../src/util/ppsub.c \
					 ../src/util/quad.c \
					 ../src/util/random.c \
					 ../src/util/roots.c \
					 ../src/util/runga.c \
					 ../src/util/screen.c \
					 ../src/util/sigplot.c \
					 ../src/util/simpleio.c \
					 ../src/util/sphhar.c \
					 ../src/util/times.c \
					 ../src/util/uinit.c \
					 ../src/util/umap.c \
					 ../src/util/vectormalloc.c \
					 ../src/util/vtk.c

libintfc_la_SOURCES = ../src/intfc/arrayutils.c \
					  ../src/intfc/comp1d.c \
					  ../src/intfc/comp2d.c \
					  ../src/intfc/comp3d.c \
					  ../src/intfc/comp.c \
					  ../src/intfc/cross2d.c \
					  ../src/intfc/geomutils.c \
					  ../src/intfc/iblkb.c \
					  ../src/intfc/iblkc2.c \
					  ../src/intfc/iblkc3.c \
					  ../src/intfc/icheck3d.c \
					  ../src/intfc/idiagnostic.c \
					  ../src/intfc/iecomps.c \
					  ../src/intfc/ifourier.c \
					  ../src/intfc/igrid.c \
					  ../src/intfc/igview.c \
					  ../src/intfc/imap.c \
					  ../src/intfc/iMesh.c \
					  ../src/intfc/iMeshP.c \
					  ../src/intfc/imkcurve.c \
					  ../src/intfc/imksurf.c \
					  ../src/intfc/int3d.c \
					  ../src/intfc/intfc.c \
					  ../src/intfc/iprt3d.c \
					  ../src/intfc/iredist.c \
					  ../src/intfc/iredist_o2.c \
					  ../src/intfc/irefl.c \
					  ../src/intfc/iscatter.c \
					  ../src/intfc/isect2d.c \
					  ../src/intfc/isect3d.c \
					  ../src/intfc/isub.c \
					  ../src/intfc/isurgery.c \
					  ../src/intfc/itopo.c \
					  ../src/intfc/itrisset.c \
					  ../src/intfc/iuserintfc.c \
					  ../src/intfc/iwallsurf.c \
					  ../src/intfc/ixgraph.c \
					  ../src/intfc/map.c \
					  ../src/intfc/ppcopy.c \
					  ../src/intfc/setb1d.c \
					  ../src/intfc/setb2d.c \
					  ../src/intfc/setb3d.c \
					  ../src/intfc/shift.c \
					  ../src/intfc/top.c \
					  ../src/intfc/triangle.c \
					  ../src/intfc/trisurf.c \
					  ../src/intfc/userhooks.c \
					  ../src/intfc/zoom.c

libfront_la_SOURCES = ../src/front/fadv3d.c \
					  ../src/front/fadv.c \
					  ../src/front/fbdry1.c \
					  ../src/front/fbdry2.c \
					  ../src/front/fbdry3.c \
					  ../src/front/fbdry4.c \
					  ../src/front/fcheck3d.c \
					  ../src/front/fcollsn3d.c \
					  ../src/front/fcorrspnd.c \
					  ../src/front/fcrosscur.c \
					  ../src/front/fcrossext.c \
					  ../src/front/fcrstatus.c \
					  ../src/front/fdiagnostic.c \
					  ../src/front/fgb2d.c \
					  ../src/front/fgb3comp.c \
					  ../src/front/fgb3d.c \
					  ../src/front/fgb3dutil.c \
					  ../src/front/fglobal.c \
					  ../src/front/fgrid.c \
					  ../src/front/finit.c \
					  ../src/front/fint.c \
					  ../src/front/fmap2.c \
					  ../src/front/fmap.c \
					  ../src/front/fnode.c \
					  ../src/front/fnodesub.c \
					  ../src/front/fpatch3d.c \
					  ../src/front/fprint.c \
					  ../src/front/fprop2d.c \
					  ../src/front/fprop3d.c \
					  ../src/front/fredist1d.c \
					  ../src/front/fredist2d.c \
					  ../src/front/fredist3d.c \
					  ../src/front/fredist.c \
					  ../src/front/frp1.c \
					  ../src/front/frp2.c \
					  ../src/front/fscat1d.c \
					  ../src/front/fscat2d.c \
					  ../src/front/fscat3d1.c \
					  ../src/front/fscat3d2.c \
					  ../src/front/fscat3d3.c \
					  ../src/front/fscatter.c \
					  ../src/front/fstate2d.c \
					  ../src/front/fstate.c \
					  ../src/front/fsub.c \
					  ../src/front/ftop.c \
					  ../src/front/funtan2d.c \
					  ../src/front/funtan3d.c \
					  ../src/front/fuserhooks.c \
					  ../src/front/fuserintfc.c \
					  ../src/front/fvelo.c

