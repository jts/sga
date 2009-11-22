bin_PROGRAMS = sga

sga_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/StringGraph

sga_LDADD = \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/SuffixTools/libsuffixtools.a \
	$(top_builddir)/Bigraph/libbigraph.a \
	$(top_builddir)/StringGraph/libstringgraph.a 

sga_SOURCES = sga.cpp \
              index.cpp index.h \
              overlap.cpp overlap.h \
              assemble.cpp assemble.h \
              oview.cpp oview.h \
              extract.cpp extract.h
