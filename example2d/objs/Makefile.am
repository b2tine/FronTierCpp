AM_CPPFLAGS = -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@


noinst_PROGRAMS = \
../example00 \
../example01 \
../example02 \
../example03 \
../example04 \
../example05 \
../example06 \
../example07 \
../example08 \
../example09 \
../example10 \
../example11 \
../example12 \
../example13 \
../example14


___example00_SOURCES = ../example00.cpp

___example01_SOURCES = ../example01.cpp

___example02_SOURCES = ../example02.cpp

___example03_SOURCES = ../example03.cpp

___example04_SOURCES = ../example04.cpp

___example05_SOURCES = ../example05.cpp

___example06_SOURCES = ../example06.cpp

___example07_SOURCES = ../example07.cpp

___example08_SOURCES = ../example08.cpp

___example09_SOURCES = ../example09.cpp

___example10_SOURCES = ../example10.cpp

___example11_SOURCES = ../example11.cpp

___example12_SOURCES = ../example12.cpp

___example13_SOURCES = ../example13.cpp

___example14_SOURCES = ../example14.cpp


AM_LDFLAGS = $(top_builddir)/lib/libFronTier.la \
			 @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@



$(PROGRAMS): $(top_builddir)/lib/libFronTier.la


$(top_builddir)/lib/libFronTier.la: FORCE
	@cd $(top_builddir)/lib && $(MAKE) libFronTier.la

FORCE:



