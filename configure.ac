#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([sarcastic], [5.0], [sarcastic@mac.com])
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/sarclib/image_stats.h])
AC_PREFIX_DEFAULT([.])
AC_CONFIG_MACRO_DIR([./macros.m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_CPP


AC_CHECK_PROG([have_gdal],[gdal-config],[yes],[no])
AM_CONDITIONAL([HAVE_GDAL], [test "x$have_gdal" = xyes])
AM_COND_IF([HAVE_GDAL],[AC_DEFINE([sarclib_HAVE_GDAL])],[AC_MSG_ERROR([required program 'gdal' not found.])])
AC_DEFINE([TRUNCATE])
AC_DEFINE([_FILE_OFFSET_BITS],[64])

AM_PROG_AR

# Checks for libraries.
AC_CHECK_LIB([readline],[history_search_pos])

AX_CHECK_CL
if test "X$no_cl" = "Xyes"; then
    AC_MSG_FAILURE([You need OpenCL])
fi


#AC_CHECK_LIB([CGAL], [refine_Delaunay_mesh_2])
# FIXME: Replace `main' with a function in `-lCGAL_Core':
#AC_CHECK_LIB([CGAL_Core], [main])
# FIXME: Replace `main' with a function in `-lboost_filesystem':
#AC_CHECK_LIB([boost_filesystem], [main])
# FIXME: Replace `main' with a function in `-lboost_system':
#AC_CHECK_LIB([boost_system], [main])
# FIXME: Replace `main' with a function in `-lboost_thread':
#AC_CHECK_LIB([boost_thread], [main])

AC_CHECK_LIB([expat], [XML_ErrorString])

AC_CHECK_LIB([fftw3f], [fftwf_assertion_failed])
AC_CHECK_LIB([fftw3f_threads], [fftwf_cleanup_threads])

# FIXME: Replace `main' with a function in `-lmpfr':
#AC_CHECK_LIB([mpfr], [main])
# FIXME: Replace `main' with a function in `-lreadline':
#AC_CHECK_LIB([readline], [main])
# FIXME: Replace `main' with a function in `-lsarclib':
#AC_CHECK_LIB([sarclib], [main])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stddef.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor gettimeofday memmove memset mkdir pow rint sqrt strcasecmp strchr strdup strncasecmp strrchr strstr])

AC_CONFIG_FILES([Makefile src/sarclib/Makefile src/SARTrace/Makefile src/sarcastic/Makefile])
AC_OUTPUT
