AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  @MPI_CFLAGS@ \
			  -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@



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
../cFweno.cpp \
../cFvisc.h \
../cFvisc.cpp \
../cFturb.h \
../cFturb.cpp \
../rigidbody.h \
../rigidbody.cpp
								

if CGAL_INTFC_OPT
libcFluid_la_LIBADD = $(top_builddir)/CGAL_Interface/objs/libCgalIntfc.la
else
libcFluid_la_LIBADD = $(top_builddir)/lib/libFronTier.la
endif


if CGAL_INTFC_OPT
$(libcFluid_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/CGAL_Interface/objs libCgalIntfc.la
else
$(libcFluid_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all
endif

FORCE:


