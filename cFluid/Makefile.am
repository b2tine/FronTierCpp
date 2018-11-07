input_directories = input-1d input-2d input-3d
loose_infiles = in-test in-rt2d in-rt3d0 in-rt3d1 in-rt3d2

EXTRA_DIST = $(input_directories) $(loose_infiles) cFweno-save.cpp


AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


noinst_PROGRAMS = cFluid

cFluid_SOURCES = cFluid.h \
				 cFbase.cpp \
				 cFcartsn.cpp \
				 cFeos.cpp \
				 cFinit.cpp \
				 cFphys.cpp \
				 cFriem.cpp \
				 cFsub.cpp \
				 cFtvd.cpp \
				 cFweno.cpp \
				 cFluid.cpp

cFluid_LDADD = $(top_builddir)/lib/libFronTier.la \
			   @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


all-local:
	@cd $(top_builddir)/lib && $(MAKE) --silent all


tagsfile:
	ctags *.h *.cpp ../include/*/*.h ../src/*/*.c



