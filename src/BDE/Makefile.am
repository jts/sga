bin_PROGRAMS = BDE

BDE_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Stats

BDE_LDADD = \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/Stats/libstats.a 


BDE_SOURCES = BDE.cpp PairStreamer.cpp
