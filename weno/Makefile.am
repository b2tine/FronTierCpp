AM_CPPFLAGS = -I$(top_srcdir)/weno \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@ @MPI_CFLAGS@



noinst_PROGRAMS = test \
				  matdata


test_SOURCES = test.cpp \
			   weno.cpp \
			   weno.h

test_LDADD = $(top_builddir)/lib/libFronTier.la


matdata_SOURCES = matdata.cpp \
				  weno.cpp \
   				  weno.h

matdata_LDADD = $(top_builddir)/lib/libFronTier.la


$(matdata_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:


tagsfile:
	ctags *.h *.cpp $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

