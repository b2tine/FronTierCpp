broken_files = example06.cpp \
			   example12.cpp \
			   example13.cpp \
			   example14.cpp \
			   example_sphere_redistribute.cpp

EXTRA_DIST = input $(broken_files)



AM_CPPFLAGS = -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

noinst_PROGRAMS = example00 \
				  example01 \
				  example02 \
				  example03 \
				  example04 \
				  example05 \
				  example07 \
				  example08 \
				  example09 \
				  example10 \
				  example11


example00_SOURCES = example00.cpp

example01_SOURCES = example01.cpp

example02_SOURCES = example02.cpp

example03_SOURCES = example03.cpp

example04_SOURCES = example04.cpp

example05_SOURCES = example05.cpp

example07_SOURCES = example07.cpp

example08_SOURCES = example08.cpp

example09_SOURCES = example09.cpp

example10_SOURCES = example10.cpp

example11_SOURCES = example11.cpp


AM_LDFLAGS = $(abs_top_builddir)/lib/libFronTier.la \
			 @MPI_LIBS@ @HDF_LIBS@ @GD_LIBS@


tagsfile:
	ctags *.cpp ../include/*/*.h ../src/*/*.c



