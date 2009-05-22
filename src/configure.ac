#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(Resolve, 0.1, js18@sanger.ac.uk)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Resolve/Resolve.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

AC_CONFIG_FILES([Makefile
		Resolve/Makefile
		SeqGraph/Makefile
		 Util/Makefile
		 UniEst/Makefile
		 Scaffold/Makefile])
AC_OUTPUT
