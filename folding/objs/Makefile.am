AM_CPPFLAGS = -I$(top_srcdir)/folding \
			  @NLOPT_CFLAGS@ @CVODE_CFLAGS@ \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CVODE_LIBS@ \
			 @NLOPT_LIBS@ \
			 @CGAL_LIBS@ \
			 @BOOST_LIBS@ \
			 @ARMA_LIBS@


noinst_LTLIBRARIES = libfolding.la



libfolding_la_SOURCES = \
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



libfolding_la_LIBADD = $(top_builddir)/Collision/objs/libCollision.la



$(libfolding_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/Collision/objs libCollision.la

FORCE:


