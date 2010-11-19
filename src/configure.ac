#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(sga, 0.9.4, js18@sanger.ac.uk)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([SGA/sga.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Ensure the sparse hash is available
AC_CHECK_HEADERS([google/sparse_hash_set],,[AC_MSG_ERROR([google sparse hash library is required: http://code.google.com/p/google-sparsehash/])])

# Check for hash table headers
AC_CHECK_HEADERS([ \
	google/sparse_hash_set google/sparse_hash_map \
	unordered_map tr1/unordered_map ext/hash_map \
	unordered_set tr1/unordered_set ext/hash_set \
])

# Check for libraries
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])

# Check for the hoard memory allocator 
AC_ARG_WITH(hoard, AS_HELP_STRING([--with-hoard=PATH],
	[specify directory containing the hoard library]))
if test "$with_hoard" -a -d "$with_hoard"; then
	hoard_ldflags="-Wl,-rpath,$with_hoard -L$with_hoard"
    enable_hoard=1
fi


# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')
AC_SUBST(LDFLAGS, "$hoard_ldflags $LDFLAGS")

if test "$enable_hoard"; then
    AC_SEARCH_LIBS([hoardmalloc], [hoard])
fi

AC_CONFIG_FILES([Makefile
		Thirdparty/Makefile
		Util/Makefile
		SQG/Makefile
		Bigraph/Makefile
		Algorithm/Makefile
		StringGraph/Makefile
		Concurrency/Makefile
		SuffixTools/Makefile
        Scaffold/Makefile
		SGA/Makefile])

AC_OUTPUT
