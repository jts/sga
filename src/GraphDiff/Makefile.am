noinst_LIBRARIES = libgraphdiff.a

libgraphdiff_a_CPPFLAGS = \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/StringGraph \
	-I$(top_srcdir)/SQG \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/Algorithm \
    -I$(top_srcdir)/Concurrency

libgraphdiff_a_SOURCES = \
        DindelUtil.h DindelUtil.cpp \
        DindelRealignWindow.h DindelRealignWindow.cpp \
        DindelHMM.h DindelHMM.cpp \
        OverlapHaplotypeBuilder.h OverlapHaplotypeBuilder.cpp \
        StringHaplotypeBuilder.h StringHaplotypeBuilder.cpp \
        DeBruijnHaplotypeBuilder.h DeBruijnHaplotypeBuilder.cpp \
        PairedDeBruijnHaplotypeBuilder.h PairedDeBruijnHaplotypeBuilder.cpp \
        HapgenUtil.h HapgenUtil.cpp \
        VCFTester.h VCFTester.cpp \
        GraphCompare.h GraphCompare.cpp \
        OverlapExtractorWithCorrection.h OverlapExtractorWithCorrection.cpp \
        DBGPathGuide.h DBGPathGuide.cpp \
        IOverlapExtractor.h
