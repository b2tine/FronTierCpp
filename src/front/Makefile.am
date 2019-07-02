SUBDIRS = objs tests


EXTRA_DIST = input \
			 testfront.c \
 			 front_gpu.cu \
 			 front_gpu.cuh


libraries:
	@$(MAKE) -C objs libfront.la



tagsfile:
	ctags $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
