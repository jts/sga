bin_PROGRAMS = sga

sga_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/StringGraph \
	-I$(top_srcdir)/Concurrency \
	-I$(top_srcdir)/Algorithm \
	-I$(top_srcdir)/SQG \
	-I$(top_srcdir)/Scaffold \
	-I$(top_srcdir)/GraphDiff \
	-I$(top_srcdir)/Thirdparty

sga_LDADD = \
	$(top_builddir)/Scaffold/libscaffold.a \
	$(top_builddir)/GraphDiff/libgraphdiff.a \
	$(top_builddir)/StringGraph/libstringgraph.a \
	$(top_builddir)/Concurrency/libconcurrency.a \
	$(top_builddir)/Algorithm/libalgorithm.a \
	$(top_builddir)/SuffixTools/libsuffixtools.a \
	$(top_builddir)/Bigraph/libbigraph.a \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/SQG/libsqg.a \
	$(top_builddir)/Thirdparty/libthirdparty.a

sga_LDFLAGS = -pthread

sga_SOURCES = sga.cpp \
              index.cpp index.h \
              overlap.cpp overlap.h \
              assemble.cpp assemble.h \
              correct.cpp correct.h \
              oview.cpp oview.h \
              preprocess.cpp preprocess.h \
              rmdup.cpp rmdup.h \
              merge.cpp merge.h \
              subgraph.cpp subgraph.h \
              scaffold.cpp scaffold.h \
              scaffold2fasta.cpp scaffold2fasta.h \
              connect.cpp connect.h \
              walk.cpp walk.h \
              filter.cpp filter.h \
              kmer-count.cpp kmer-count.h \
              stats.cpp stats.h \
              fm-merge.cpp fm-merge.h \
              gmap.h gmap.cpp \
              filterBAM.h filterBAM.cpp \
              cluster.h cluster.cpp \
              gen-ssa.h gen-ssa.cpp \
              bwt2fa.h bwt2fa.cpp \
              graph-diff.h graph-diff.cpp \
              graph-concordance.h graph-concordance.cpp \
              gapfill.h gapfill.cpp \
              preqc.h preqc.cpp \
              overlap-long.h overlap-long.cpp \
              variant-detectability.h variant-detectability.cpp \
              rewrite-evidence-bam.h rewrite-evidence-bam.cpp \
              somatic-variant-filters-bam.h somatic-variant-filters.cpp \
              haplotype-filter.h haplotype-filter.cpp \
              OverlapCommon.h OverlapCommon.cpp \
              SGACommon.h 
