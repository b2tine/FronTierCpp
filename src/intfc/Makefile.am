SUBDIRS = objs tests


EXTRA_DIST = testintfc.c \
			 intfc_gpu.cuh \
			 intfc_gpu.cu


libraries:
	@$(MAKE) -C objs libintfc.la


tagsfile:
	ctags $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
