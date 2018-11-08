AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


noinst_PROGRAMS = ../cFluid

___cFluid_SOURCES = \
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
../cFluid.cpp


___cFluid_LDADD = $(top_builddir)/lib/libFronTier.la \
   				  @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


$(___cFluid_DEPENDENCIES): FORCE
	@cd $(top_builddir)/lib && $(MAKE) all

FORCE:


