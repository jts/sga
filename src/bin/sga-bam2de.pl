#! /usr/bin/perl

use strict;
use Getopt::Long;

my $n = 5;
my $minLength = 0;
my $prefix = "";
my $numThreads = 1;

GetOptions("prefix=s" => \$prefix,
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
$cmd = "abyss-fixmate -h $prefix.hist $bamFile | samtools view -Sb - > $prefix.diffcontigs.bam";
runCmd($cmd);

# sort 
$cmd = "samtools sort $prefix.diffcontigs.bam $prefix.diffcontigs.sorted";
runCmd($cmd);

# distance est
$cmd = "DistanceEst -n $n -k 99 -j $numThreads -s $minLength -o $prefix.de $prefix.hist $prefix.diffcontigs.sorted.bam";
runCmd($cmd);

sub usage
{
    print "sga-bam2de.pl - Make a distance estimate file from a bam file of reads aligned to contigs\n";
    print "sga-bam2de.pl -n N --prefix lib300 lib300.bam\n";
    print "Options:\n";
    print "                -n N             Minimum number of pairs required to consider two contigs linked\n";
    print "                -m LEN           Only find links between contigs with length at least LEN bp\n";
    print "                --prefix NAME    Use NAME as the prefix for the outfiles\n";
}

sub runCmd
{
    my($c) = @_;
    print "$c\n";
    system($c);
}
