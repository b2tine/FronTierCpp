AM_CPPFLAGS = -I$(top_srcdir)/FrontWrapper \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@



noinst_LTLIBRARIES = ../../lib/libFronTierCpp.la



______lib_libFronTierCpp_la_SOURCES = \
../FrontWrapper.h



______lib_libFronTierCpp_la_LIBADD = $(top_builddir)/lib/libFronTier.la




$(______lib_libFronTierCpp_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/src && $(MAKE) all

FORCE:


clean-libtool:
