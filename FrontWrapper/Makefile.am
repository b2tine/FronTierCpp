SUBDIRS = objs tests


AM_CPPFLAGS = -I$(top_srcdir)/FrontWrapper \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


noinst_PROGRAMS = example00

example00_SOURCES = example00.cpp

example00_LDADD = $(top_builddir)/lib/libFronTierCpp.la




tagsfile:
	ctags *.cpp *.h $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
