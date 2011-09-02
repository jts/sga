#! /usr/bin/perl
#
# Rename the sga index files for the given set of reads
#
use strict;
use File::Basename;

my $reads = $ARGV[0];
my $outprefix = $ARGV[1];
die("error: prefix cannot be null") if $outprefix eq "";
my ($file, $path, $suffix) = fileparse($reads, (".fa", ".fastq"));

my $inbase = $path . $file;
my $outbase = $path . $outprefix;

# Make sure the prefix is not already used
die("error: prefix $outprefix already in use") if(-f "$outbase.bwt");

foreach my $ext (".sai", ".rsai", ".bwt", ".rbwt", $suffix)
{
    run("mv $inbase$ext $outbase$ext");
}

sub run
{
    my ($cmd) = @_;
    print $cmd . "\n";
    system($cmd);
}
