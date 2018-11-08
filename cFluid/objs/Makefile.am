AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = -static
AM_LIBTOOLFLAGS = --preserve-dup-deps


noinst_LTLIBRARIES = ../../lib/libcFluid.la



______lib_libcFluid_la_SOURCES = \
../cFluid.h \
../cFbase.cpp \
../cFcartsn.cpp \
../cFeos.cpp \
../cFinit.cpp \
../cFphys.cpp \
../cFriem.cpp \
../cFsub.cpp \
../cFtvd.cpp \
../cFweno.cpp
								

______lib_libcFluid_la_LIBADD = $(top_builddir)/lib/libFronTier.la



$(______lib_libcFluid_la_DEPENDENCIES): FORCE
	@cd $(top_builddir)/src && $(MAKE) all

FORCE:
