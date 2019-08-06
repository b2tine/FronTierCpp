AM_CPPFLAGS = -I$(top_srcdir)/cFluid \
			  -I$(top_srcdir)/include \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@



if CUDA_OPT

NVCC = @NVCC@
AM_CPPFLAGS += @CUDA_CFLAGS@ 
AM_LDFLAGS = @CUDA_LIBS@

endif



LIB_FT = $(top_builddir)/lib/libFronTier.la


noinst_LTLIBRARIES = libcFluid.la


libcFluid_la_SOURCES = \
../cFluid.h \
../cFbase.cpp \
../cFcartsn.cpp \
../cFeos.cpp \
../cFinit.cpp \
../cFphys.cpp \
../cFriem.cpp \
../cFsub.cpp \
../cFtvd.cpp \
../cFweno.cpp
								

libcFluid_la_LIBADD = $(LIB_FT)


$(libcFluid_la_DEPENDENCIES): FORCE
	@$(MAKE) -C $(top_builddir)/lib all

FORCE:



if CUDA_OPT

noinst_LIBRARIES = libwenocuda.a

libwenocuda_a_SOURCES = ../cFweno_gpu.cu \
					   ../cFweno_gpu.cuh

.cu.o:
	$(NVCC) -c $<


libcFluid_la_LIBADD += $(builddir)/libwenocuda.a

endif




