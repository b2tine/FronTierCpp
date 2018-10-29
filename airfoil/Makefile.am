AM_CPPFLAGS = -I$(abs_top_srcdir)/airfoil \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @CGAL_CFLAGS@ @PETSC_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

noinst_PROGRAMS = airfoil

airfoil_LDADD = @CGAL_LIBS@ -L$(abs_top_builddir)/lib -liFluid \
			   @PETSC_LIBS@ @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

airfoil_SOURCES = airfoil.h \
				  airfoil.cpp \
				  afcnpy.cpp \
				  afdata.cpp \
				  afinit.cpp \
				  afinit3d.cpp \
				  afprop.cpp \
				  afsetd.cpp \
				  aftest.cpp \
				  afvelo.cpp \
				  cgal.cpp
