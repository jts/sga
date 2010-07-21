#! /usr/bin/perl
use strict;

#
# scaffoldDriver.pl CONTIGS READS
#
# Build scaffolds from the CONTIGS using the pairs in the READS file
#
my $contigsFile = $ARGV[0];
my $readsFile = $ARGV[1];
my $kmer = 100;
my $n = 10;

# Output files
my $splitReads1 = "splitReads.1.fa";
my $splitReads2 = "splitReads.2.fa";
my $sacFile1 = $contigsFile . ".1.sac";
my $sacFile2 = $contigsFile . ".2.sac";
my $bamSEFile = $contigsFile . ".se.bam";
my $bamFile = $contigsFile . ".bam";
my $finalPrefix = $contigsFile . ".final";
my $finalBam = "$finalPrefix.bam";
my $histFile = $contigsFile . ".hist";
my $deFile = $contigsFile . ".de";

# Program paths
my $splitReads = "/nfs/team71/phd/js18/work/devel/sga/tools/splitReads.pl";
my $bwa = "/software/solexa/bin/aligners/bwa/current/bwa";
my $samtools = "/nfs/team71/phd/js18/software/samtools/samtools-dev/samtools"; # this uses Shaun Jackman's patch to calculate ISIZE
my $parseAligns = "/nfs/team71/phd/js18/bin/ParseAligns";
my $distanceEst = "/nfs/team71/phd/js18/bin/DistanceEst";

if(1) {
# Index the contigs
run("$bwa index $contigsFile");

# Split the input reads file into two files, one for each half of the pair
run("$splitReads $readsFile $splitReads1 $splitReads2");

# Align the reads to the contigs and generate a sam file
# Do not allow indels
# Do not use sampe because ParseAligns is used to do the pairing
run("$bwa aln -o 0 $contigsFile $splitReads1 > $sacFile1");
run("$bwa aln -o 0 $contigsFile $splitReads2 > $sacFile2");

run("$bwa sampe -s $contigsFile $sacFile1 $sacFile2 $splitReads1 $splitReads2 | $samtools view -Sb -o $bamSEFile -");
unlink($sacFile1);
unlink($sacFile2);

# Use Shaun Jackman's patched samtools fixmate to fill in the ISIZE field for the pairs
run("$samtools fixmate $bamSEFile $bamFile");
unlink("$bamSEFile");

# Run ParseAligns to generate a new SAM file with only inter-contig pairs
# and the fragment size histogram
run("$parseAligns -k $kmer -h $histFile --sam $bamFile > /dev/null");

#
run("$samtools sort $bamFile $finalPrefix");
run("$samtools index $finalBam");
unlink($bamFile);

# Run DistanceEst to generate estimates between the contigs
run("$distanceEst -n $n -k $kmer $histFile $finalBam > $deFile");

}

sub run
{
    my($cmd) = @_;
    print STDERR $cmd . "\n";
    system($cmd);
}
