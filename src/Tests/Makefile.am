bin_PROGRAMS = Tests

Tests_CPPFLAGS = \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/Util 

Tests_LDADD = \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/Bigraph/libbigraph.a

Tests_SOURCES = Tests.cpp
