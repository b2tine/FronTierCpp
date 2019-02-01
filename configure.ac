AC_PREREQ([2.69])
#Initialize autoconf with project name and version number
AC_INIT([FronTier],[1.1])

#Include macro directory
AC_CONFIG_MACRO_DIR([m4])
#Directory for generated build tools
AC_CONFIG_AUX_DIR([build-aux])
#Header file with generated #DEFINE macros.
#   Must be included at the top of include/util/cdecs.h
AC_CONFIG_HEADERS([config.h])

#Initialize automake and specify foreign (non-GNU) standards.
AM_INIT_AUTOMAKE([foreign --warnings=none])
AM_SILENT_RULES([yes])


#Initialize Libtool
AM_PROG_AR
LT_INIT([disable-shared])

#Check for compilers
AC_PROG_CXX([mpic++ mpicxx])


#Check for pthreads
AX_PTHREAD([have_pthread=yes],[have_pthread=no])

GTEST_VERSION="1.8.0"
AS_IF([test "x${have_pthread}" = xyes],
      [GTEST_CPPFLAGS="-I\$(top_srcdir)/lib/gmock -DGTEST_HAS_PTHREADS=1"
       GTEST_CXXFLAGS="$PTHREAD_CFLAGS"
       GTEST_LDFLAGS="-L\$(top_builddir)/lib/gmock/lib"
       GTEST_LIBS="$PTHREAD_LIBS"],
       [GTEST_CPPFLAGS="-I\$(top_srcdir)/lib/gmock -DGTEST_HAS_PTHREADS=0"
        GTEST_CXXFLAGS=
        GTEST_LDFLAGS="-L\$(top_builddir)/lib/gmock/lib"
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
       [test "x${hdf4_inc}" != xno && test "x${hdf4_lib}" != xno],
       [with_hdf=yes
        hdfinc="-I${hdf4_inc}"
        hdflib="-L${hdf4_lib} ${hdf4_libflags}"])

AS_IF([test $with_hdf = yes],
      [AC_DEFINE([HAVE_HDF4],[1],[build with hdf4])])

AC_SUBST(HDF_CFLAGS,${hdfinc})
AC_SUBST(HDF_LIBS,${hdflib})


mpi_libflags="-lmpi"

#configure with MPI present in user environment
AC_ARG_WITH([mpi],[AS_HELP_STRING([--with-mpi],
            [use when mpi provided by user environment])],
            [with_mpi=yes
             mpiinc=
             mpilib=${mpi_libflags}],
            [with_mpi=no])

#configure with OpenMPI
AC_ARG_WITH([openmpi-dir],[AS_HELP_STRING([--with-openmpi-dir=path],
            [specify location of single directory installation])],
            [openmpi_dir=${withval}],[openmpi_dir=no])

AS_IF([test "x${openmpi_dir}" != xno],
      [with_openmpi=yes
       mpiinc="-I${openmpi_dir}/include"
       mpilib="-L${openmpi_dir}/lib ${mpi_libflags}"],
       [with_openmpi=no])

#configure with MPICH
AC_ARG_WITH([mpich-dir],[AS_HELP_STRING([--with-mpich-dir=path],
            [specify location of single directory installation])],
            [mpich_dir=${withval}],[mpich_dir=no])

AS_IF([test "x${mpich_dir}" != xno],
      [with_mpich=yes
       mpiinc="-I${mpich_dir}/include"
       mpilib="-L${mpich_dir}/lib ${mpi_libflags}"],
       [with_mpich=no])

AC_SUBST(MPI_CFLAGS,${mpiinc})
AC_SUBST(MPI_LIBS,${mpilib})

AS_IF([test $with_mpi = yes || test $with_openmpi = yes || test $with_mpich = yes],
      [AC_DEFINE([HAVE_MPI],[1],[build with mpi])])


#TODO: more flexibility in hypre-petsc configuration

#configure with HYPRE (NOT IN USE YET)
AC_ARG_WITH([hypre-dir],[AS_HELP_STRING([--with-hypre-dir=path],
            [specify location of single directory installation])],
            [hypre_dir=${withval}],[hypre_dir=no])

hypre_libflags="-lHYPRE"

with_hypre=no
AS_IF([test "x${hypre_dir}" != xno],
      [with_hypre=yes
       hypreinc="-I${hypre_dir}/include"
       hyprelib="-L${hypre_dir}/lib ${hypre_libflags}"])

#--with-hypre-petsc-dir option

#configure with PETSC
AC_ARG_WITH([petsc-dir],[AS_HELP_STRING([--with-petsc-dir=path],
            [specify location of single directory installation])],
            [petsc_dir=${withval}],[petsc_dir=no])

petsc_libflags="-lpetsc -lHYPRE -llapack -lblas -ldl -fopenmp -lm -L/usr/X11R6/lib -lX11"

with_petsc=no
AS_IF([test "x${petsc_dir}" != xno],
      [with_petsc=yes
       petscinc="-I${petsc_dir}/include"
       petsclib="-L${petsc_dir}/lib ${petsc_libflags}"])

AC_SUBST(PETSC_CFLAGS,${petscinc})
AC_SUBST(PETSC_LIBS,${petsclib})

#hypre should only be defined if it was specified ... 
AS_IF([test $with_petsc = yes],
      [AC_DEFINE([HAVE_HYPRE],[1],[build with HYPRE])])

#check if PETSC dependent modules should be built
AM_CONDITIONAL([PETSC_OPT],[test $with_petsc = yes])

AM_COND_IF([PETSC_OPT],
           [AC_CONFIG_FILES([solver/Makefile
                             solver/objs/Makefile
                             iFluid/Makefile
                             iFluid/objs/Makefile
                             iFluid/tests/Makefile
                             melting/Makefile
                             crystal/Makefile
                             subsurf/Makefile])])

#configure with BOOST
AC_ARG_WITH([boost-dir],[AS_HELP_STRING([--with-boost-dir=path],
            [specify location of single directory installation])],
            [boost_dir=${withval}],[boost_dir=no])

AC_ARG_WITH([boost-inc],[AS_HELP_STRING([--with-boost-inc=path],
            [specify location of header files])],
            [boost_inc=${withval}],[boost_inc=no])

AC_ARG_WITH([boost-lib],[AS_HELP_STRING([--with-boost-lib=path],
            [specify location of libraries])],
            [boost_lib=${withval}],[boost_lib=no])


boost_libflags=

with_boost=no
AS_IF([test "x${boost_dir}" != "xno"],
      [with_boost=yes
       boostinc="-I${boost_dir}/include"
       boostlib="-L${boost_dir}/lib ${boost_libflags}"],
       [test "x${boost_inc}" != xno && test "x${boost_lib}" != xno],
       [with_boost=yes
        boostinc="-I${boost_inc}"
        boostlib="-L${boost_lib} ${boost_libflags}"])

AC_SUBST(BOOST_CFLAGS,${boostinc})
AC_SUBST(BOOST_LIBS,${boostlib})

AS_IF([test $with_boost = yes],
      [AC_DEFINE([HAVE_BOOST],[1],[build with boost])])



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

cgal_libflags="-lCGAL_Core -lCGAL_ImageIO -lCGAL -lgmp -lmpfr"

with_cgal=no
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



#configure with Armadillo
AC_ARG_WITH([arma-dir],[AS_HELP_STRING([--with-arma-dir=path],
            [specify location of single directory installation])],
            [arma_dir=${withval}],[arma_dir=no])

AC_ARG_WITH([arma-inc],[AS_HELP_STRING([--with-arma-inc=path],
            [specify location of header files])],
            [arma_inc=${withval}],[arma_inc=no])

AC_ARG_WITH([arma-lib],[AS_HELP_STRING([--with-arma-lib=path],
            [specify location of libraries])],
            [arma_lib=${withval}],[arma_lib=no])


arma_cflags=-DARMA_DONT_USE_WRAPPER
arma_libflags="-lblas -llapack"
#TODO: add option to specify -lopenblas intead.

with_arma=no
AS_IF([test "x${arma_dir}" != xno],
      [with_arma=yes
       armainc="-I${arma_dir}/include ${arma_cflags}"
       armalib="-L${arma_dir}/lib ${arma_libflags}"],
       [test "x${arma_inc}" != xno && test "x${arma_lib}" != xno],
       [with_arma=yes
        armainc="-I${arma_inc} ${arma_cflags}"
        armalib="-L${arma_lib} ${arma_libflags}"])

#The macro HAVE_ARMA isn't being used right now, but may be useful for debugging
AS_IF([test $with_arma = yes],
      [AC_DEFINE([HAVE_ARMA],[1],[build with armadillo])])

AC_SUBST(ARMA_CFLAGS,${armainc})
AC_SUBST(ARMA_LIBS,${armalib})


#configure with CVODE
AC_ARG_WITH([cvode-dir],[AS_HELP_STRING([--with-cvode-dir=path],
            [specify location of single directory installation])],
            [cvode_dir=${withval}],[cvode_dir=no])

AC_ARG_WITH([cvode-inc],[AS_HELP_STRING([--with-cvode-inc=path],
            [specify location of header files])],
            [cvode_inc=${withval}],[cvode_inc=no])

AC_ARG_WITH([cvode-lib],[AS_HELP_STRING([--with-cvode-lib=path],
            [specify location of libraries])],
            [cvode_lib=${withval}],[cvode_lib=no])

cvode_libflags="-lsundials_cvode -lsundials_nvecserial"

with_cvode=no
AS_IF([test "x${cvode_dir}" != "xno"],
      [with_cvode=yes
       cvodeinc="-I${cvode_dir}/include"
       cvodelib="-L${cvode_dir}/lib ${cvode_libflags}"],
       [test "x${cvode_inc}" != xno && test "x${cvode_lib}" != xno],
       [with_cvode=yes
        cvodeinc="-I${cvode_inc}"
        cvodelib="-L${cvode_lib} ${cvode_libflags}"])

AC_SUBST(CVODE_CFLAGS,${cvodeinc})
AC_SUBST(CVODE_LIBS,${cvodelib})

#The macro HAVE_CVODE isn't being used right now, but may be useful for debugging
AS_IF([test $with_cvode = yes],
      [AC_DEFINE([HAVE_CVODE],[1],[build with cvode])])


#configure with NLOPT
AC_ARG_WITH([nlopt-dir],[AS_HELP_STRING([--with-nlopt-dir=path],
            [specify location of single directory installation])],
            [nlopt_dir=${withval}],[nlopt_dir=no])

AC_ARG_WITH([nlopt-inc],[AS_HELP_STRING([--with-nlopt-inc=path],
            [specify location of header files])],
            [nlopt_inc=${withval}],[nlopt_inc=no])

AC_ARG_WITH([nlopt-lib],[AS_HELP_STRING([--with-nlopt-lib=path],
            [specify location of libraries])],
            [nlopt_lib=${withval}],[nlopt_lib=no])

nlopt_libflags="-lnlopt"

with_nlopt=no
AS_IF([test "x${nlopt_dir}" != "xno"],
      [with_nlopt=yes
       nloptinc="-I${nlopt_dir}/include"
       nloptlib="-L${nlopt_dir}/lib ${nlopt_libflags}"],
       [test "x${nlopt_inc}" != xno && test "x${nlopt_lib}" != xno],
       [with_nlopt=yes
        nloptinc="-I${nlopt_inc}"
        nloptlib="-L${nlopt_lib} ${nlopt_libflags}"])

AC_SUBST(NLOPT_CFLAGS,${nloptinc})
AC_SUBST(NLOPT_LIBS,${nloptlib})

#The macro HAVE_NLOPT isn't being used right now, but may be useful for debugging
AS_IF([test $with_nlopt = yes],
      [AC_DEFINE([HAVE_NLOPT],[1],[build with nlopt])])



#check if airfoil and Collsion should be built
AM_CONDITIONAL([CGAL_ARMA_OPT],
               [test $with_cgal = yes && test $with_arma = yes])


AM_COND_IF([CGAL_ARMA_OPT],
           [AM_COND_IF([PETSC_OPT],
                       [AC_CONFIG_FILES([Collision/Makefile
                                         Collision/objs/Makefile
                                         airfoil/Makefile
                                         airfoil/objs/Makefile
                                         parachute/Makefile])])])

#check if folding should be built
AM_COND_IF([CGAL_ARMA_OPT],
           [AM_CONDITIONAL([FOLDING_OPT],
                           [test $with_cvode = yes && test $with_nlopt = yes])])

AM_COND_IF([FOLDING_OPT],[AC_CONFIG_FILES([folding/Makefile
                                           folding/objs/Makefile])])



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
                 lib/Makefile
                 lib/gmock/Makefile
                 src/Makefile
                 src/util/Makefile
                 src/util/objs/Makefile
                 src/intfc/Makefile
                 src/intfc/objs/Makefile
                 src/front/Makefile
                 src/front/objs/Makefile
                 cFluid/Makefile
                 cFluid/objs/Makefile
                 cFluid/tests/Makefile
                 example2d/Makefile
                 example2d/objs/Makefile
                 example3d/Makefile
                 example3d/objs/Makefile
                 tests/Makefile])

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
