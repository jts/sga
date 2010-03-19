#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(sga, 0.1, js18@sanger.ac.uk)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([SGA/sga.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

AC_CONFIG_FILES([Makefile
		Util/Makefile
		Bigraph/Makefile
		Algorithm/Makefile
		StringGraph/Makefile
		Concurrency/Makefile
		SuffixTools/Makefile
		Tests/Makefile
		SGA/Makefile])

AC_OUTPUT
