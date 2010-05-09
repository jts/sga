bin_PROGRAMS = sga

sga_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/StringGraph \
	-I$(top_srcdir)/Concurrency \
	-I$(top_srcdir)/Algorithm \
	-I$(top_srcdir)/SQG \
	-I$(top_srcdir)/Thirdparty

sga_LDADD = \
    $(top_builddir)/StringGraph/libstringgraph.a \
    $(top_builddir)/Algorithm/libalgorithm.a \
	$(top_builddir)/SuffixTools/libsuffixtools.a \
	$(top_builddir)/Bigraph/libbigraph.a \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/Concurrency/libconcurrency.a \
	$(top_builddir)/SQG/libsqg.a \
	$(top_builddir)/Thirdparty/libthirdparty.a

sga_LDFLAGS = -pthread

sga_SOURCES = sga.cpp \
              index.cpp index.h \
              overlap.cpp overlap.h \
              assemble.cpp assemble.h \
              oview.cpp oview.h \
			  preprocess.cpp preprocess.h
