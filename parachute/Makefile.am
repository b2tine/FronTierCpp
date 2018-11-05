input_directories = input-C9 \
					input-cross \
					input-G11 \
					input-intruder \
					input-mchutes \
					input-others \
					input-T10 \
					input-T11

EXTRA_DIST = $(input_directories)

AM_CPPFLAGS = -I$(abs_top_srcdir)/parachute \
			  -I$(abs_top_srcdir)/airfoil \
			  -I$(abs_top_srcdir)/Collision \
			  @CGAL_CFLAGS@ \
			  -I$(abs_top_srcdir)/iFluid \
			  -I$(abs_top_srcdir)/solver \
			  @PETSC_CFLAGS@ @MPI_CFLAGS@ \
			  -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @HDF_CFLAGS@ @GD_CFLAGS@

AM_CPPFLAGS += -DDAMPING_FORCE -DCOLLISION_DETECTION


noinst_PROGRAMS = parachute

parachute_SOURCES = modules.cpp \
					parachute.cpp
				  
parachute_LDADD = $(abs_top_builddir)/lib/libairfoil.la \
				  $(abs_top_builddir)/lib/libCollision.la \
  				  $(abs_top_builddir)/lib/libiFluid.la \
  				  $(abs_top_builddir)/lib/libsolver.la \
  				  @CGAL_LIBS@ @PETSC_LIBS@ @MPI_LIBS@ \
				  $(abs_top_builddir)/lib/libFronTier.la \
  				  @HDF_LIBS@ @GD_LIBS@



tagsfile:
	ctags *.h *.cpp ../Collision/*.h ../Collision/*.cpp ../iFluid/*.h ../iFluid/*.cpp ../solver/*.h ../solver/*.cpp ../include/*/*.h ../src/*/*.c


clean-local:
	-rm -rf tags
	-rm -rf TAGS
	-rm -rf CTAGS
	-rm -rf Makefile

