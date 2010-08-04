bin_PROGRAMS = sga_scaffold

sga_scaffold_CPPFLAGS = \
	-I$(top_srcdir)/Util \
    -I$(top_srcdir)/Thirdparty \
    -I$(top_srcdir)/Bigraph

sga_scaffold_LDADD = \
	$(top_builddir)/Util/libutil.a \
    $(top_builddir)/Thirdparty/libthirdparty.a

sga_scaffold_SOURCES = sga-scaffold.cpp \
                       ScaffoldVertex.h ScaffoldVertex.cpp \
                       ScaffoldEdge.h ScaffoldEdge.cpp \
                       ScaffoldGraph.h ScaffoldGraph.cpp \
                       ScaffoldVisitors.h ScaffoldVisitors.cpp
