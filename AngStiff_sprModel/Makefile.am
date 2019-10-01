SUBDIRS = objs tests


AM_CPPFLAGS = -I$(top_srcdir)/AngStiff_sprModel \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = spring

				  
spring_SOURCES = spring.cpp


spring_LDADD =	$(builddir)/objs/libairfoil.la



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/Collision/*.h $(top_srcdir)/Collision/*.cpp \
		$(top_srcdir)/iFluid/*.h $(top_srcdir)/iFluid/*.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c



