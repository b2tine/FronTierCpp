#Initialize autoconf with project name and version number
AC_INIT([FronTier],[1.0])

#Directory for generated build tools
AC_CONFIG_AUX_DIR([build-aux])
#Directory for generated macros
#AC_CONFIG_MACRO_DIR([m4])

#Initialize automake and specify foreign (non-GNU) standards.
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

#To create static library
AC_PROG_RANLIB

#Check for C compiler
AC_PROG_CC
#Check for C++ compiler
AC_PROG_CXX
#Check for ar 
AM_PROG_AR

#Header files to produce that include #DEFINE statements
#AC_CONFIG_HEADERS([config.h])

AC_ARG_WITH([hdf],[AS_HELP_STRING([--with-hdf],
            [use hdf (default is no)])],[],[use_hdf=no])

#Makefiles to produce
#AC_CONFIG_FILES([Makefile libs/Makefile example2d/Makefile])
AC_CONFIG_FILES([Makefile libs/Makefile])

#Calls ./configure in example2d recursively
#AC_CONFIG_SUBDIRS([example2d])

#Generate output
AC_OUTPUT

