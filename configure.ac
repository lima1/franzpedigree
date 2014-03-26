#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(FRANz, [2.0.0], [markus@bioinf.uni-leipzig.de] )
AC_CONFIG_SRCDIR([src/dataio.h])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h stdbool.h string.h math.h sys/time.h getopt.h])
# libtap
AC_CHECK_HEADERS([pthread.h float.h limits.h])
AC_CHECK_LIB(pthread, main)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.

# uncomment AC_FUNC_MALLOC on MINGW
AC_FUNC_MALLOC
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([memset memmove memcpy pow exp log sqrt strcasecmp strchr strrchr atexit isnan])
AC_CHECK_FUNCS([vasprintf flockfile funlockfile asprintf snprintf strstr floor])

AC_CONFIG_FILES([Makefile
                 libdir/Makefile
                 libdir/dcmt/Makefile
                 libdir/dcmt/lib/Makefile
                 libdir/dcmt/include/Makefile
                 examples/Makefile
                 extras/Makefile
                 extras/input/Makefile
                 extras/installer/Makefile
                 extras/installer/win/Makefile
                 extras/installer/mac/Makefile
                 man/Makefile
                 t/Makefile
                 src/Makefile])

# Enable support for OpenMP
AC_OPENMP([C])
CFLAGS="$CFLAGS $OPENMP_CFLAGS -std=c99 -DADJLIST"
LDFLAGS="$LDFLAGS $OPENMP_LDFLAGS"
if test "$enable_openmp" = "yes"; then
    AC_DEFINE(HAVE_OPENMP,1,[Enable openmp])
fi
AC_ARG_ENABLE([debug],
AC_HELP_STRING([--enable-debug], [Enable debug code (default=no)]),
[enable_debug=$enableval],
[enable_debug=no])
if test "x$enable_debug" = xyes; then
    CFLAGS="$CFLAGS -g  -O0 -pedantic -Wall -fstack-protector -D_FORTIFY_SOURCE=2"
    #CFLAGS="$CFLAGS -DHASH_EMIT_KEYS=3"
    AC_DEFINE(DEBUG,1,[Enable debug code])
else    
    CFLAGS="$CFLAGS -O3 -pedantic -Wall"
#    AC_DEFINE(DEBUG,0,[Disable debug code])
fi

# Mingw
#CFLAGS="$CFLAGS -DWIN32"

AC_ARG_ENABLE([gprof],
AC_HELP_STRING([--enable-gprof], [Enable profiling (default=no)]),
[enable_gprof=$enableval],
[enable_gprof=no])
if test "x$enable_gprof" = xyes; then
    CFLAGS="$CFLAGS -pg "
fi

AC_ARG_ENABLE([assertions],
AC_HELP_STRING([--enable-assertions], [Enable assertions (default=no)]),
[enable_assertions=$enableval],
[enable_assertions=no])
if test "x$enable_assertions" = xno; then
    AC_DEFINE(NDEBUG, [], [Don't execute assertions])
fi

AC_COVERAGE()

AC_PROG_RANLIB
AC_OUTPUT

