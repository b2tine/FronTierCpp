SUBDIRS = objs

AM_CPPFLAGS = -I$(top_srcdir)/CGAL_Interface \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@


noinst_PROGRAMS = cgal_intfc_test \
				  cgal_capsule_intfc_test \
				  cgal_cylindrical_shell

cgal_intfc_test_SOURCES = cgal_intfc_test.cpp
cgal_intfc_test_LDADD = $(builddir)/objs/libCgalIntfc.la

cgal_capsule_intfc_test_SOURCES = cgal_capsule_intfc_test.cpp
cgal_capsule_intfc_test_LDADD = $(builddir)/objs/libCgalIntfc.la

cgal_cylindrical_shell_SOURCES = cgal_cylindrical_shell.cpp
cgal_cylindrical_shell_LDADD = $(builddir)/objs/libCgalIntfc.la

tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
