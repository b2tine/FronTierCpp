AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@



noinst_PROGRAMS = timeSeries

timeSeries_SOURCES = timeSeries.h timeSeries.cpp

timeSeries_LDADD = $(top_builddir)/lib/libFronTier.la



$(timeSeries_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:



tagsfile:
	ctags *.h *.cpp $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
