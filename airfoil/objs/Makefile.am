AM_CPPFLAGS = -I$(top_srcdir)/airfoil \
			  -I$(top_srcdir)/Collision \
			  @CGAL_CFLAGS@ @BOOST_CFLAGS@ \
			  @ARMA_CFLAGS@ \
			  -I$(top_srcdir)/iFluid \
			  -I$(top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(top_srcdir)/include \
			  @HDF_CFLAGS@ @GD_CFLAGS@

if CUDA_OPT

NVCC = @NVCC@
AM_CPPFLAGS += @CUDA_CFLAGS@ 
AM_LDFLAGS = @CUDA_LIBS@

endif

#AM_CPPFLAGS += -DDAMPING_FORCE



noinst_LTLIBRARIES = libairfoil.la

libairfoil_la_SOURCES = \
../airfoil.h \
../airfoil_sv.h \
../afcnpy.cpp \
../afdata.cpp \
../afinit3d.cpp \
../afinit.cpp \
../afprop.cpp \
../afsetd.cpp \
../aftest.cpp \
../afvelo.cpp \
../cgal.cpp

libairfoil_la_LIBADD = $(top_builddir)/Collision/objs/libCollision.la


$(libairfoil_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/Collision/objs libCollision.la

FORCE:


if CUDA_OPT

noinst_LTLIBRARIES += libaircuda.la

libaircuda_la_SOURCES = ../airfoil_gpu.cu \
					   ../airfoil_gpu.cuh

.cu.o:
	$(NVCC) -c $<

.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(NVCC) -c $<


libairfoil_la_LIBADD += $(top_builddir)/airfoil/objs/libaircuda.a

endif




