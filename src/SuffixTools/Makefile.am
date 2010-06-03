noinst_LIBRARIES = libsuffixtools.a

libsuffixtools_a_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Concurrency

libsuffixtools_a_SOURCES = BWT.cpp BWT.h \
                           STCommon.cpp STCommon.h \
                           STGlobals.h \
                           Occurrence.cpp Occurrence.h \
                           SuffixArray.cpp SuffixArray.h \
                           SuffixCompare.cpp SuffixCompare.h \
                           InverseSuffixArray.cpp InverseSuffixArray.h \
                           SACAInducedCopying.h SACAInducedCopying.cpp \
                           BWTAlgorithms.h BWTAlgorithms.cpp \
						   BWTDiskConstruction.h BWTDiskConstruction.cpp \
						   BWTReader.h BWTReader.cpp \
						   BWTWriter.h BWTWriter.cpp \
						   SAReader.h SAReader.cpp \
						   SAWriter.h SAWriter.cpp \
						   GapArray.h GapArray.cpp \
						   RankProcess.h RankProcess.cpp \
                           RLBWT.h RLBWT.cpp
