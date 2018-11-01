AM_CPPFLAGS = -I$(abs_top_srcdir)/airfoil \
			  -I$(abs_top_srcdir)/Collision \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @CGAL_CFLAGS@ @PETSC_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @MPI_CFLAGS@ @HDF_CFLAGS@ @GD_CFLAGS@

AM_CXXFLAGS = -DDAMPING_FORCE -DCOLLISION_DETECTION


noinst_PROGRAMS = airfoil

airfoil_LDADD = -L$(abs_top_builddir)/lib -liFluid \
				@CGAL_LIBS@ @PETSC_LIBS@ @MPI_LIBS@ \
				@HDF_LIBS@ @GD_LIBS@

noinst_HEADERS = airfoil.h \
				 airfoil_sv.h \
				 airfoil_gpu.cuh \
				 ../Collision/collid.h

airfoil_SOURCES = airfoil.cpp \
				  afcnpy.cpp \
				  afdata.cpp \
				  afinit.cpp \
				  afinit3d.cpp \
				  afprop.cpp \
				  afsetd.cpp \
				  aftest.cpp \
				  afvelo.cpp \
				  cgal.cpp \
				  airfoil_gpu.cu \
				  ../Collision/dcollid.cpp \
				  ../Collision/dcollid3d.cpp


tagsfile:
	ctags *.h *.cpp ../Collision/*.h ../Collision/*.cpp ../iFluid/*.h ../iFluid/*.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c


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


EXTRA_DIST = $(input_directories)
