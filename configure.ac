#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(sarcastic, 5.0, sarcastic@mac.com)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR(src/sarclib/image_stats.h)
AC_CONFIG_HEADERS(config.h)

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_CPP

# Checks for libraries.
AC_CHECK_LIB([CGAL], [refine_Delaunay_mesh_2])
# FIXME: Replace `main' with a function in `-lCGAL_Core':
AC_CHECK_LIB([CGAL_Core], [main])
# FIXME: Replace `main' with a function in `-lboost_filesystem':
AC_CHECK_LIB([boost_filesystem], [main])
# FIXME: Replace `main' with a function in `-lboost_system':
AC_CHECK_LIB([boost_system], [main])
# FIXME: Replace `main' with a function in `-lboost_thread':
AC_CHECK_LIB([boost_thread], [main])
# FIXME: Replace `main' with a function in `-lexpat':
AC_CHECK_LIB([expat], [main])
# FIXME: Replace `main' with a function in `-lfftw3f':
AC_CHECK_LIB([fftw3f], [main])
# FIXME: Replace `main' with a function in `-lmpfr':
AC_CHECK_LIB([mpfr], [main])
# FIXME: Replace `main' with a function in `-lreadline':
AC_CHECK_LIB([readline], [main])
# FIXME: Replace `main' with a function in `-lsarclib':
AC_CHECK_LIB([sarclib], [main])

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

AC_CONFIG_FILES([src/Common/tinyxml/Makefile
                 src/SARTrace/Makefile
                 src/bircs/Makefile
                 src/colladaToPlyFile/Makefile
                 src/materialise/Makefile
                 src/sarcastic/Makefile
                 src/sarclib/Makefile])
AC_OUTPUT
