AM_CPPFLAGS = -I$(top_srcdir)/weno \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@



noinst_PROGRAMS = test \
				  matdata \
				  convection_diffusion


test_SOURCES = test.cpp \
			   weno.cpp \
			   weno.h

test_LDADD = $(top_builddir)/lib/libFronTier.la

$(test_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all


matdata_SOURCES = matdata.cpp \
				  weno.cpp \
   				  weno.h

matdata_LDADD = $(top_builddir)/lib/libFronTier.la



convection_diffusion_SOURCES = convection_diffusion.cpp \
			 				   implicit.cpp \
			 				   implicit.h \
			 				   weno.cpp \
			 				   weno.h

convection_diffusion_LDADD = $(top_builddir)/solver/objs/libsolver.la


FORCE:


tagsfile:
	ctags *.h *.cpp \
		$(top_srcdir)/solver/*.h $(top_srcdir)/solver/*.cpp \
		$(top_srcdir)/include/*/*.h $(top_srcdir)/src/*/*.c

