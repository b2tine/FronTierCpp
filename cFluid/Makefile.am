AM_CPPFLAGS = -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @HDF_CFLAGS@ @MPI_CFLAGS@

AM_LDFLAGS = -L$(abs_top_builddir)/lib -lfront -lintfc -lutil -lcFluid\
			 @HDF_LIBS@ @MPI_LIBS@

noinst_PROGRAMS = cFluid

cFluid_SOURCES = cFluid.cpp

