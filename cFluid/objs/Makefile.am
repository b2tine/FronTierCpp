AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


LIB_FT = $(top_builddir)/lib/libFronTier.la


noinst_LTLIBRARIES = libcFluid.la


libcFluid_la_SOURCES = \
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
								

libcFluid_la_LIBADD = $(LIB_FT)


$(LIB_FT): libraries

libraries:
	@$(MAKE) -C $(top_builddir)/lib libFronTier.la


