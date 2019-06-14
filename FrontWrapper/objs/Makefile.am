AM_CPPFLAGS = -I$(top_srcdir)/FrontWrapper \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@



noinst_LTLIBRARIES = libFronTierCpp.la



libFronTierCpp_la_SOURCES = \
../FrontWrapper.h \
../FrontWrapper.cpp



libFronTierCpp_la_LIBADD = $(top_builddir)/lib/libFronTier.la




$(libFronTierCpp_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib libFronTier.la

FORCE:

