SUBDIRS = objs


AM_CPPFLAGS = -I$(top_srcdir)/AngStiff_sprModel \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_PROGRAMS = spring

				  
spring_SOURCES = spring.cpp


spring_LDADD = $(builddir)/objs/libangstiff.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



