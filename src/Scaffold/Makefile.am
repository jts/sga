bin_PROGRAMS = Scaffold

Scaffold_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/SeqGraph

Scaffold_LDADD = \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/SeqGraph/libseqgraph.a

Scaffold_SOURCES = Scaffold.cpp ScaffoldData.cpp
