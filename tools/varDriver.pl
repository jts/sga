#! /usr/bin/perl
use strict;
use File::Basename;

#
# scaffoldDriver.pl CONTIGS READS
#
# Build scaffolds from the CONTIGS using the pairs in the READS file
#
my $contigsFile = $ARGV[0];
my $readsFile1 = $ARGV[1];
my $readsFile2 = $ARGV[2];

if($contigsFile eq "" || $readsFile1 eq "" || $readsFile2 eq "")
{
    print "Usage: varDriver.pl <contigsFile> <readsFile1> <readsFile2\n";
    exit(1);
}

# Output files
my $sacFile1 = $readsFile1 . ".sac";
my $sacFile2 = $readsFile2 . ".sac";
my $tmpBam1 = $readsFile1 . ".tmp.bam";
my $tmpBam2 = $readsFile2 . ".tmp.bam";

my $bamFile1 = $readsFile1 . ".bam";
my $bamFile2 = $readsFile2 . ".bam";

# Program paths
my $toolsDir = dirname($0);
my $bwa = "bwa"; #"/software/solexa/bin/aligners/bwa/current/bwa";
my $samtools = "samtools"; #"/nfs/team71/phd/js18/software/samtools/samtools-dev/samtools"; # this uses Shaun Jackman's patch to calculate ISIZE

# Index the contigs
run("$bwa index $contigsFile");

# Align the reads to the contigs and generate a sam file
# Do not allow indels
# Do not use sampe because ParseAligns is used to do the pairing
run("$bwa aln $contigsFile $readsFile1 > $sacFile1");
run("$bwa aln $contigsFile $readsFile2 > $sacFile2");

run("$bwa samse $contigsFile $sacFile1 $readsFile1 | $samtools view -Sb -o $tmpBam1 -");
run("$bwa samse $contigsFile $sacFile2 $readsFile2 | $samtools view -Sb -o $tmpBam2 -");
unlink($sacFile1);
unlink($sacFile2);

#
run("$samtools sort $tmpBam1 $readsFile1");
run("$samtools sort $tmpBam2 $readsFile2");
unlink($tmpBam1);
unlink($tmpBam2);

run("$samtools index $bamFile1");
run("$samtools index $bamFile2");

sub run
{
    my($cmd) = @_;
    print STDERR $cmd . "\n";
    system($cmd);
}
