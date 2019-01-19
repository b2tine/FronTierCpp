AM_CPPFLAGS = -I$(top_srcdir)/folding \
			  -I$(top_srcdir)/folding/nlopt \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_CPPFLAGS += -DDAMPING_FORCE -DCOLLISION_DETECTION


AM_LDFLAGS = -L$(top_srcdir)/folding/nlopt -lnlopt



noinst_LTLIBRARIES = ../../lib/libfolding.la



______lib_libfolding_la_SOURCES = \
../bending.cpp \
../bending.h \
../cgal.cpp \
../cgal.h \
../di_spring_solver.cpp \
../drag.cpp \
../drag.h \
../drag_proto.cpp \
../drag_proto.h \
../ex_spring_solver.cpp \
../folding.cpp \
../folding.h \
../folding_helper.cpp \
../im_spring_solver.cpp \
../origami.cpp \
../origami.h \
../singleton.h \
../spring_solver.cpp \
../spring_solver.h



______lib_libfolding_la_LIBADD = $(top_builddir)/lib/libCollision.la




$(______lib_libfolding_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/Collision/objs && $(MAKE) all

FORCE:



clean-libtool:
