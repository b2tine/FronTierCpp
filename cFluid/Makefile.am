input_directories = input-1d input-2d input-3d
EXTRA_DIST = $(input_directories)

AM_CPPFLAGS = -I$(abs_top_srcdir)/cFluid \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

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

cFluid_LDADD = $(abs_top_builddir)/lib/libFronTier.la

tagsfile:
	ctags *.h *.cpp ../include/*/*.h ../src/*/*.c


