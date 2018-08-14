#Initialize autoconf with project name and version number
AC_INIT([FronTier],[1.0])

#Store generated build tools
AC_CONFIG_AUX_DIR([build-aux])

#Initialize automake and specify foreign (non-GNU) standards.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

#To create static library
AC_PROG_RANLIB

#Check for C compiler
AC_PROG_CC

#NOTE: will most likely want to use autoheader for this task
#Header files to produce that include #DEFINE statements
#AC_CONFIG_HEADERS([config.h])

#Makefiles to produce
AC_CONFIG_FILES([Makefile src/util/Makefile])

#Generate output
AC_OUTPUT

