SUBDIRS = objs tests


libraries:
	@$(MAKE) -C objs libutil.la


tagsfile:
	ctags $(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
