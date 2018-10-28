AM_CPPFLAGS = -I$(abs_top_srcdir)/cFluid \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @HDF_CFLAGS@ @MPI_CFLAGS@


noinst_PROGRAMS = cFluid

cFluid_LDADD = -L$(abs_top_builddir)/lib -lcFluid \
  			   @HDF_LIBS@ @MPI_LIBS@

cFluid_SOURCES = cFluid.cpp

