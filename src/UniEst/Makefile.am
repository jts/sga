bin_PROGRAMS = UniEst

UniEst_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Stats

UniEst_LDADD = \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/Stats/libstats.a \
	$(GSL_LIBS)

UniEst_SOURCES = UniEst.cpp UniData.cpp
