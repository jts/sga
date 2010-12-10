noinst_LIBRARIES = libconcurrency.a

libconcurrency_a_CPPFLAGS = \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/Algorithm

libconcurrency_a_SOURCES = \
        OverlapProcess.h OverlapProcess.cpp \
        RmdupProcess.h RmdupProcess.cpp \
        FMMergeProcess.h FMMergeProcess.cpp \
        SequenceProcessFramework.h \
        SequenceWorkItem.h \
        ThreadWorker.h \
		MkqsThread.h
