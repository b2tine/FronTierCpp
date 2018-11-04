AM_CPPFLAGS = -I$(abs_top_srcdir)/cFluid \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@

noinst_LTLIBRARIES = ../lib/libcFluid.la

___lib_libcFluid_la_CPPFLAGS = -I$(top_srcdir)/cFluid $(AM_CPPFLAGS)

___lib_libcFluid_la_SOURCES = cFluid.h \
							  cFbase.cpp \
							  cFcartsn.cpp \
							  cFeos.cpp \
							  cFinit.cpp \
							  cFluid.cpp \
							  cFphys.cpp \
							  cFriem.cpp \
							  cFsub.cpp \
							  cFtvd.cpp \
							  cFweno.cpp \
							  cFweno-save.cpp

___lib_libcFluid_la_LIBADD = $(abs_top_builddir)/lib/libFronTier.la


noinst_PROGRAMS = cFluid

cFluid_DEPENDENCIES = $(abs_top_builddir)/lib/libcFluid.la \
					  $(abs_top_builddir)/lib/libFronTier.la

cFluid_SOURCES = cFluid.h \
				 cFluid.cpp

cFluid_LDADD = -L$(abs_top_builddir)/lib -lcFluid \
			   $(AM_LDFLAGS)
			   


tagsfile:
	ctags *.h *.cpp ../include/*/*.h ../src/*/*.c


input_directories = input-1d input-2d input-3d

EXTRA_DIST = $(input_directories)
