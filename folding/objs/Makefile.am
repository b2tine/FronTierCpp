AM_CPPFLAGS = -I$(top_srcdir)/folding \
			  @NLOPT_CFLAGS@ \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/Collision \
			  -I$(top_srcdir)/iFluid \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_CPPFLAGS += -DDAMPING_FORCE


AM_LDFLAGS = @CVODE_LIBS@ \
			 @NLOPT_LIBS@ \
			 @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@


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
