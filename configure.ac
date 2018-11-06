AC_PREREQ([2.69])
#Initialize autoconf with project name and version number
AC_INIT([FronTier],[1.1])

#Include macros from the autoconf-archive
AC_CONFIG_MACRO_DIR([m4])
#Directory for generated build tools
AC_CONFIG_AUX_DIR([build-aux])
#Header file with generated #DEFINE macros.
#   included at the top of include/util/cdecs.h
AC_CONFIG_HEADERS([config.h])

#Initialize automake and specify foreign (non-GNU) standards.
AM_INIT_AUTOMAKE([foreign subdir-objects])

#Initialize Libtool
AM_PROG_AR
LT_INIT([disable-shared])
LT_LANG([C++])

#Check for compilers
AC_PROG_CXX([mpic++ mpicxx g++])
AC_LANG([C++])

#Check for pthreads
AX_PTHREAD([have_pthread=yes],[have_pthread=no])

GTEST_VERSION="1.8.0"
AS_IF([test "x${have_pthread}" = xyes],
      [GTEST_CPPFLAGS="-I\$(abs_top_srcdir)/lib/gmock -DGTEST_HAS_PTHREADS=1"
       GTEST_CXXFLAGS="$PTHREAD_CFLAGS"
       GTEST_LDFLAGS="-L\$(abs_top_builddir)/lib/gmock/lib"
       GTEST_LIBS="$PTHREAD_LIBS"],
       [GTEST_CPPFLAGS="-I\$(abs_top_srcdir)/lib/gmock -DGTEST_HAS_PTHREADS=0"
        GTEST_CXXFLAGS=
        GTEST_LDFLAGS="-L\$(abs_top_builddir)/lib/gmock/lib"
        GTEST_LIBS=])

AC_SUBST([GTEST_VERSION])
AC_SUBST([GTEST_CPPFLAGS])
AC_SUBST([GTEST_CXXFLAGS])
AC_SUBST([GTEST_LDFLAGS])
AC_SUBST([GTEST_LIBS])

#configure with GD
AC_ARG_WITH([gd-dir],[AS_HELP_STRING([--with-gd-dir=path],
            [specify location of single directory installation])],
            [gd_dir=${withval}],[gd_dir=no])

AC_ARG_WITH([gd-inc],[AS_HELP_STRING([--with-gd-inc=path],
            [specify location of header files])],
            [gd_inc=${withval}],[gd_inc=no])

AC_ARG_WITH([gd-lib],[AS_HELP_STRING([--with-gd-lib=path],
            [specify location of libraries])],
            [gd_lib=${withval}],[gd_lib=no])

with_gd=no
AS_IF([test "x${gd_dir}" != xno],
      [with_gd=yes
       gdinc="-I${gd_dir}/include"
       gdlib="-L${gd_dir}/lib -lgd"],
       [test "x${gd_inc}" != xno && test "x${gd_lib}" != xno],
       [with_gd=yes
        gdinc="-I${gd_inc}"
        gdlib="-L${gd_lib} -lgd"])

AS_IF([test $with_gd = yes],
      [AC_DEFINE([HAVE_GD],[1],[build with gd])])

AC_SUBST(GD_CFLAGS,${gdinc})
AC_SUBST(GD_LIBS,${gdlib})

#configure with HDF4
AC_ARG_WITH([hdf4-dir],[AS_HELP_STRING([--with-hdf4-dir=path],
            [specify location of single directory installation])],
            [hdf4_dir=${withval}],[hdf4_dir=no])

AC_ARG_WITH([hdf4-inc],[AS_HELP_STRING([--with-hdf4-inc=path],
            [specify location of header files])],
            [hdf4_inc=${withval}],[hdf4_inc=no])

AC_ARG_WITH([hdf4-lib],[AS_HELP_STRING([--with-hdf4-lib=path],
            [specify location of libraries])],
            [hdf4_lib=${withval}],[hdf4_lib=no])

with_hdf=no
hdf4_libflags="-lmfhdf -ldf -ljpeg -lz"

AS_IF([test "x${hdf4_dir}" != xno],
      [with_hdf=yes
       hdfinc="-I${hdf4_dir}/include"
       hdflib="-L${hdf4_dir}/lib ${hdf4_libflags}"],
       [test "x${hdf4_inc}" != xno && test "x$h{df4_lib}" != xno],
       [with_hdf=yes
        hdfinc="-I${hdf4_inc}"
        hdflib="-L${hdf4_lib} ${hdf4_libflags}"])

AS_IF([test $with_hdf = yes],
      [AC_DEFINE([HAVE_HDF4],[1],[build with hd4])])

AC_SUBST(HDF_CFLAGS,${hdfinc})
AC_SUBST(HDF_LIBS,${hdflib})

#configure with MPI present in user environment
AC_ARG_WITH([mpi],[AS_HELP_STRING([--with-mpi],
            [use when paths to mpi are part of user environment])],
            [with_mpi=yes],[with_mpi=no])


#configure with MPICH
AC_ARG_WITH([mpich-dir],[AS_HELP_STRING([--with-mpich-dir=path],
            [specify location of single directory installation])],
            [mpich_dir=${withval}],[mpich_dir=no])

with_mpich=no
AS_IF([test "x${mpich_dir}" != xno],
      [with_mpich=yes
       mpichinc="-I${mpich_dir}/include"
       mpichlib="-L${mpich_dir}/lib -lmpi -lmpicxx"])

AC_SUBST(MPI_CFLAGS,${mpichinc})
AC_SUBST(MPI_LIBS,${mpichlib})

AS_IF([test $with_mpich = yes],
      [AC_DEFINE([HAVE_MPI],[1],[build with mpich])])

#configure with PETSC (must be built with HYPRE)
AC_ARG_WITH([petsc-dir],[AS_HELP_STRING([--with-petsc-dir=path],
            [specify location of single directory installation])],
            [petsc_dir=${withval}],[petsc_dir=no])

with_petsc=no
petsc_libflags="-lpetsc -lHYPRE -llapack -lblas -ldl -fopenmp -lm -L/usr/X11R6/lib -lX11"

AS_IF([test "x${petsc_dir}" != xno],
      [with_petsc=yes
       petscinc="-I${petsc_dir}/include"
       petsclib="-L${petsc_dir}/lib ${petsc_libflags}"])

AC_SUBST(PETSC_CFLAGS,${petscinc})
AC_SUBST(PETSC_LIBS,${petsclib})

AS_IF([test $with_petsc = yes],
      [AC_DEFINE([HAVE_HYPRE],[1],[build with HYPRE])])

#check if PETSC dependent modules should be built
AM_CONDITIONAL([PETSC_OPT],[test $with_petsc = yes])

AM_COND_IF([PETSC_OPT],
           [AC_CONFIG_FILES([iFluid/Makefile
                            melting/Makefile
                            crystal/Makefile])])

#configure with CGAL
AC_ARG_WITH([cgal-dir],[AS_HELP_STRING([--with-cgal-dir=path],
            [specify location of single directory installation])],
            [cgal_dir=${withval}],[cgal_dir=no])

AC_ARG_WITH([cgal-inc],[AS_HELP_STRING([--with-cgal-inc=path],
            [specify location of header files])],
            [cgal_inc=${withval}],[cgal_inc=no])

AC_ARG_WITH([cgal-lib],[AS_HELP_STRING([--with-cgal-lib=path],
            [specify location of libraries])],
            [cgal_lib=${withval}],[cgal_lib=no])

with_cgal=no
cgal_libflags="-lCGAL_Core -lCGAL_ImageIO -lCGAL -lgmp -lmpfr"

AS_IF([test "x${cgal_dir}" != "xno"],
      [with_cgal=yes
       cgalinc="-I${cgal_dir}/include -frounding-math"
       cgallib="-L${cgal_dir}/lib ${cgal_libflags}"],
       [test "x${cgal_inc}" != xno && test "x${cgal_lib}" != xno],
       [with_cgal=yes
        cgalinc="-I${cgal_inc} -frounding-math"
        cgallib="-L${cgal_lib} ${cgal_libflags}"])

AC_SUBST(CGAL_CFLAGS,${cgalinc})
AC_SUBST(CGAL_LIBS,${cgallib})

#The macro HAVE_CGAL isn't being used right now, but may be useful for debugging
AS_IF([test $with_cgal = yes],
      [AC_DEFINE([HAVE_CGAL],[1],[build with cgal])])

#check if airfoil should be built
AM_CONDITIONAL([CGAL_OPT],
               [test $with_cgal = yes])

AM_COND_IF([CGAL_OPT],
           [AM_COND_IF([PETSC_OPT],
                       [AC_CONFIG_FILES(
                        [Collision/Makefile
                         airfoil/Makefile
                         parachute/Makefile])])])

#enable/disable Collision module ....



#optionally compile the examples (Makefiles always generated)
AC_ARG_WITH([example2d],[AS_HELP_STRING([--with-example2d],
              [build the 2d example directory])],
              [exTwoD=yes],[exTwoD=no])

AM_CONDITIONAL([EX2D_OPT],[test $exTwoD = yes])

AC_ARG_WITH([example3d],[AS_HELP_STRING([--with-example3d],
              [build the 3d example directory])],
              [exThreeD=yes],[exThreeD=no])

AM_CONDITIONAL([EX3D_OPT],[test $exThreeD = yes])

#Output files
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/util/Makefile
                 src/intfc/Makefile
                 src/front/Makefile
                 lib/Makefile
                 lib/gmock/Makefile
                 tests/Makefile
                 cFluid/Makefile
                 example2d/Makefile
                 example3d/Makefile])

#Generate output
AC_OUTPUT

AC_MSG_RESULT([
FronTier++ $VERSION is now configured
Configure Information:
  C Compiler        : $CC
    DEFS            :   $DEFS
    CPPFLAGS        :   $CPPFLAGS
    CFLAGS          :   $CFLAGS

  C++ Compiler      : $CXX
    DEFS            :   $DEFS
    CPPFLAGS        :   $CPPFLAGS
    CXXFLAGS        :   $CXXFLAGS

  Linker            : $LD
    LDFLAGS         :   $LDFLAGS
    LIBS            :   $LIBS

  Google Test 
    GTEST_CPPFLAGS  : $GTEST_CPPFLAGS
    GTEST_CXXFLAGS  : $GTEST_CXXFLAGS
    GTEST_LDFLAGS   : $GTEST_LDFLAGS
    GTEST_LIBS      : $GTEST_LIBS
])
