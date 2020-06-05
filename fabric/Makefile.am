SUBDIRS = objs
#SUBDIRS += tests


input_directories = input-2d \
					input-7ft \
					input-C9 \
					input-cloth \
					input-collsn \
					input-Cross \
					input-drum \
					input-G11 \
					input-G12D \
					input-intruder \
					input-others \
					input-T10 \
					input-young3d

EXTRA_DIST = $(input_directories) airfoil_gpu.cuh  airfoil_gpu.cu



AM_CPPFLAGS = -I$(top_srcdir)/fabric \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@


AM_LDFLAGS = @CGAL_LIBS@ @BOOST_LIBS@ @ARMA_LIBS@


noinst_PROGRAMS = fabric_Cplus \
				  cgal_intfc_test \
				  spring \
				  fabric_Fortn

				  
fabric_Cplus_SOURCES = fabric_Cplus.cpp
fabric_Cplus_LDADD = $(builddir)/objs/libfabric.la

cgal_intfc_test_SOURCES = cgal_intfc_test.cpp
cgal_intfc_test_LDADD = $(builddir)/objs/libfabric.la

spring_SOURCES = spring.cpp
spring_LDADD = $(builddir)/objs/libfabric.la



fabric_Fortn_SOURCES = fortran_api.F90 \
					   fabric_Fortn.F90

fabric_Fortn_LDADD = $(builddir)/objs/libfabric.la -lstdc++

%.o %.mod: %.F90
	$(AM_V_PPFC)$(PPFCCOMPILE) -c $<

fabric_Fortn.o: fortran_api.mod

clean-local:
	$(RM) *.mod



tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c
