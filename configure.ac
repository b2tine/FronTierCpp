AC_PREREQ([2.69])
#Initialize autoconf with project name and version number
AC_INIT([FronTier],[1.0])

#Include macros from the autoconf-archive
AC_CONFIG_MACRO_DIR([autoconf-archive/m4])
#Directory for generated build tools
AC_CONFIG_AUX_DIR([build-aux])
#Header file with generated #DEFINE statements
AC_CONFIG_HEADERS([build-aux/config.h])

#Initialize automake and specify foreign (non-GNU) standards.
AM_INIT_AUTOMAKE([foreign subdir-objects])

#Initialize Libtool
AM_PROG_AR
LT_INIT([static])

#Check for compilers
#AC_PROG_CC([gcc])
AC_PROG_CXX([g++])
AC_LANG([C++])
#AC_LANG([C],[C++])

#Check for pthreads
ACX_PTHREAD([have_pthread=yes],[have_pthread=no])

GTEST_VERSION="1.8.0"
AS_IF([test "x$have_pthread" = "xyes"],
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

#configure with hdf4
AC_ARG_WITH([hdf4-dir],[AS_HELP_STRING([--with-hdf4-dir=path],
            [specify location of single directory installation)])],
            [hdf4_dir=$withval],[hdf4_dir="no"])

AC_ARG_WITH([hdf4-inc],[AS_HELP_STRING([--with-hdf4-inc=path],
            [specify location of header files])],
            [hdf4_inc=$withval],[hdf4_inc="no"])

AC_ARG_WITH([hdf4-lib],[AS_HELP_STRING([--with-hdf4-lib=path],
            [specify location of libraries])],
            [hdf4_lib=$withval],[hdf4_lib="no"])

with_hdf=no
AS_IF([test "x$hdf4_dir" != "xno"],
      [with_hdf=yes
       hdfinc="-I${hdf4_dir}/include"
       hdflib="-L${hdf4_dir}/lib"],
       [test "x$hdf4_inc" != "xno" && test "x$hdf4_lib" != "xno"],
       [with_hdf=yes
        hdfinc="-I${hdf4_inc}"
        hdflib="-L${hdf_lib}"])

AS_IF([test $with_hdf = yes],
      [AC_DEFINE([USE_HDF4],[1],[building with hd4])
      AC_SUBST(HDF_CFLAGS,${hdfinc})
      AC_SUBST(HDF_LIBS,${hdflib}" -lmfhdf -ldf -ljpeg -lz")])


AC_ARG_WITH([example2d],[AS_HELP_STRING([--with-example2d],
              [build the 2d example directory])],
              [exTwoD="yes"],[exTwoD="no"])

AM_CONDITIONAL([EX2D_OPT],[test "$exTwoD" = "yes"])

AC_ARG_WITH([example3d],[AS_HELP_STRING([--with-example3d],
              [build the 3d example directory])],
              [exThreeD="yes"],[exThreeD="no"])

AM_CONDITIONAL([EX3D_OPT],[test "$exThreeD" = "yes"])

#Output files
AC_CONFIG_FILES([Makefile
                 lib/Makefile
                 lib/gmock/Makefile
                 tests/Makefile
                 example2d/Makefile
                 example3d/Makefile])

#Generate output
AC_OUTPUT

AC_MSG_RESULT([
FronTier_Automake $VERSION is now configured
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
