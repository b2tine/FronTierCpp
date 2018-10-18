AM_CPPFLAGS = -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  -I$(abs_top_srcdir)/solver \
			  @HDF_CFLAGS@ @MPI_CFLAGS@ @PETSC_CFLAGS@

AM_LDFLAGS = -L$(abs_top_builddir)/lib -lfront -lintfc -lutil \
			 @HDF_LIBS@ @MPI_LIBS@ @PETSC_LIBS@

noinst_PROGRAMS = airfoil

airfoil_SOURCES =
