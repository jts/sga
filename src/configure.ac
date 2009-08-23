#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(sga, 0.1, js18@sanger.ac.uk)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Resolve/Resolve.cpp])
AC_CONFIG_HEADER([config.h])

# Options to configure.
AC_ARG_WITH(gsl, AS_HELP_STRING([--with-gsl=PATH],
        [specify prefix directory for the gnu scientific library]))

if test "$with_gsl" -a -d "$with_gsl"; then
        gsl_cppflags="-I$with_gsl/include"
        gsl_ldflags="-L$with_gsl/lib"
fi

# Set the compiler flags to include the gsl path
AC_SUBST(CPPFLAGS, "$CPPFLAGS $gsl_cppflags")
AC_SUBST(LDFLAGS, "$LDFLAGS $gsl_ldflags")

# Check for the GNU scientific library
AM_PATH_GSL

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

AC_CONFIG_FILES([Makefile
		Util/Makefile
		Bigraph/Makefile
		Scaffold/Makefile
		Stats/Makefile
		UniEst/Makefile
		BDE/Makefile
		Resolve/Makefile
		IndelFinder/Makefile
		SuffixTools/Makefile
		Tests/Makefile
		SGA/Makefile])

AC_OUTPUT
