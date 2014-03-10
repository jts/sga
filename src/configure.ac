#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(sga, 0.10.13, js18@sanger.ac.uk)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([SGA/sga.cpp])
AC_CONFIG_HEADER([config.h])
m4_include([m4/m4_ax_openmp.m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for libraries
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )

# Check for openmp
AX_OPENMP([openmp_cppflags="-fopenmp" AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])

# Check for the bamtools library path
# Bamtools has two different path formats
# If it is built in place the files will be in {path}/lib/ and {path}/include/
# If it is installed (with make install) it will be in {path}/include/bamtools
# and {path}/lib/bamtools. Both variants are handled here
AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))

if test "$with_bamtools" -a -d "$with_bamtools"; then
    bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib -Wl,-rpath,$with_bamtools/lib/bamtools -L$with_bamtools/lib/bamtools"
    bamtools_include="-I$with_bamtools/include -I$with_bamtools/include/bamtools"
fi

# Support for tcmalloc/jemalloc/hoard

# Check for the jemalloc memory allocator
AC_ARG_WITH(jemalloc, AS_HELP_STRING([--with-jemalloc=PATH],
	[specify directory containing the jemalloc library]))

# Check for the tcmalloc
AC_ARG_WITH(tcmalloc, AS_HELP_STRING([--with-tcmalloc=PATH],
	[specify directory containing the tcmalloc library]))

# Check for the hoard memory allocator
AC_ARG_WITH(hoard, AS_HELP_STRING([--with-hoard=PATH],
	[specify directory containing the hoard memory allocator library]))

# Set library path to user-selected allocator
if test "$with_hoard" -a -d "$with_hoard"; then
	external_malloc_ldflags="-Wl,-rpath,$with_hoard -L$with_hoard"
    enable_hoard=1
fi

if test "$with_tcmalloc" -a -d "$with_tcmalloc"; then
	external_malloc_ldflags="-Wl,-rpath,$with_tcmalloc -L$with_tcmalloc"
    enable_tcmalloc=1
fi

if test "$with_jemalloc" -a -d "$with_jemalloc"; then
	external_malloc_ldflags="-Wl,-rpath,$with_jemalloc -L$with_jemalloc"
    enable_jemalloc=1
fi

# Check for the google sparse hash
AC_ARG_WITH(sparsehash, AS_HELP_STRING([--with-sparsehash=PATH],
	[specify directory containing the google sparsehash headers http://code.google.com/p/google-sparsehash/)]))

if test "$with_sparsehash" -a -d "$with_sparsehash"; then
    sparsehash_include="-I$with_sparsehash/include"
fi

# Warn that multithreading is not available on macosx, since it does not implement unnamed semaphores
AC_MSG_CHECKING(for host type)
host="`uname -a | awk '{print $1}'`";
if test "$host" = Darwin;then
  AC_MSG_RESULT(warning: multi-threading is not available since OSX does not support un-named pthread semaphores.)
else
  AC_MSG_RESULT(you are not using osx so multi-threading should work.);
fi

# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra $fail_on_warning -Wno-unknown-pragmas")
AC_SUBST(CXXFLAGS, "-O3")
AC_SUBST(CFLAGS, "-O3")
AC_SUBST(CPPFLAGS, "$CPPFLAGS $openmp_cppflags $sparsehash_include $bamtools_include")
AC_SUBST(LDFLAGS, "$openmp_cppflags $external_malloc_ldflags $bamtools_ldflags $LDFLAGS")

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS -lbamtools")

#
if test "$enable_hoard"; then
    AC_CHECK_LIB([hoard], [malloc])
fi

if test "$enable_tcmalloc"; then
    AC_CHECK_LIB([tcmalloc], [malloc])
fi

if test "$enable_jemalloc"; then
    AC_CHECK_LIB([jemalloc], [malloc])
fi


# Ensure the sparse hash is available
AC_CHECK_HEADERS([google/sparse_hash_set],,[AC_MSG_ERROR([google sparse hash library is required: http://code.google.com/p/google-sparsehash/])])

# Check for hash table headers
AC_CHECK_HEADERS([ \
	google/sparse_hash_set google/sparse_hash_map \
	unordered_map tr1/unordered_map ext/hash_map \
	unordered_set tr1/unordered_set ext/hash_set \
])

# Make sure the bamtools headers can be found
AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

AC_CONFIG_FILES([Makefile
        bin/Makefile
		Thirdparty/Makefile
		Util/Makefile
		SQG/Makefile
		Bigraph/Makefile
		Algorithm/Makefile
		StringGraph/Makefile
		Concurrency/Makefile
		SuffixTools/Makefile
        GraphDiff/Makefile
        Scaffold/Makefile
		SGA/Makefile])

AC_OUTPUT
