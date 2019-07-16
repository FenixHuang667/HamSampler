#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Tsampler], [1.0], [fenixprotoss@gmail.com])
AC_CONFIG_SRCDIR([src/pair_mat.h])
AM_INIT_AUTOMAKE
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow strchr strdup strstr])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
