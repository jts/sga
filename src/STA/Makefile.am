bin_PROGRAMS = STA

STA_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Stats

STA_LDADD = \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/Stats/libstats.a 


STA_SOURCES = STA.cpp \
			  BWT.cpp BWT.h \
			  STCommon.cpp STCommon.h \
			  SuffixTree.cpp SuffixTree.h \
			  STGlobal.h \
			  Occurance.cpp Occurance.h \
			  SuffixArray.cpp SuffixArray.h \
			  InverseSuffixArray.cpp InverseSuffixArray.h \
			  ReadTable.cpp ReadTable.h
