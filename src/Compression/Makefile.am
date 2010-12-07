noinst_LIBRARIES = libcompression.a

libcompression_a_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Thirdparty \
    -I$(top_srcdir)/SuffixTools
libcompression_a_SOURCES = HuffmanTreeCodec.h \
						   HuffmanUtil.h HuffmanUtil.cpp \
						   PackedTableDecoder.h \
                           StreamEncoding.h \
						   HuffmanForest.h HuffmanForest.cpp \
                           BWTCompressor.h BWTCompressor.cpp
