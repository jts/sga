#! /usr/bin/perl

use strict;
use Getopt::Long;

my $n = 5;
my $k = 99;
my $minLength = 0;
my $prefix = "";
my $numThreads = 1;

# Filter the abyss distance est histogram to remove insert sizes
# with fewer than hist_min data points
my $hist_min = 3;

GetOptions("prefix=s" => \$prefix,
           "k=i"      => \$k,
           "n=i"      => \$n,
           "m=i"      => \$minLength,
           "t=i"      => \$numThreads);


my $bFail = 0;
if($prefix eq "")
{
    print("Error: A prefix must be specified\n");
    $bFail = 1;
}
my $bamFile = $ARGV[0];
if($bamFile eq "")
{
    print "Error: A bam file must be provided\n";
    $bFail = 1;
}

if($bFail)
{
    usage();
    exit(1);
}

my $cmd;

# fix-mate info
$cmd = "abyss-fixmate -h $prefix.tmp.hist $bamFile | samtools view -Sb - > $prefix.diffcontigs.bam";
runCmd($cmd);

# DistanceEst is very slow when the learned insert size distribution is very wide
# To work around this, remove singleton data points from the histogram
runCmd("awk \'\$2 >= $hist_min\' $prefix.tmp.hist > $prefix.hist");

# sort 
$cmd = "samtools sort $prefix.diffcontigs.bam $prefix.diffcontigs.sorted";
runCmd($cmd);

# distance est
$cmd = "DistanceEst -s $minLength -n $n -k $k -j $numThreads -s $minLength -o $prefix.de $prefix.hist $prefix.diffcontigs.sorted.bam";
runCmd($cmd);

sub usage
{
    print "sga-bam2de.pl - Make a distance estimate file from a bam file of reads aligned to contigs\n";
    print "sga-bam2de.pl -n N --prefix lib300 lib300.bam\n";
    print "Options:\n";
    print "                -n N             Minimum number of pairs required to consider two contigs linked\n";
    print "                -m LEN           Only find links between contigs with length at least LEN bp\n";
    print "                -t NUM           Use NUM threads for computing the distance estimates\n";
    print "                --prefix NAME    Use NAME as the prefix for the outfiles\n";
}

sub runCmd
{
    my($c) = @_;
    print "$c\n";
    system($c);
}
