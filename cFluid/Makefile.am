input_directories = input-1d input-2d input-3d
loose_infiles = in-test in-rt2d in-rt3d0 in-rt3d1 in-rt3d2

EXTRA_DIST = $(input_directories) $(loose_infiles) cFweno-save.cpp


SUBDIRS = objs


tagsfile:
	ctags *.h *.cpp $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

