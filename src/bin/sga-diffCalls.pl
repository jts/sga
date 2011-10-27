#! /usr/bin/perl

use strict;
use Getopt::Long;

my $baseFile;
my $variantFile;

GetOptions("variant=s" => \$variantFile,
           "base=s"    => \$baseFile);

die("A --base file must be provided") if($baseFile eq "");
die("A --variant file must be provided") if($variantFile eq "");

open(B, $baseFile) || die("Cannot open $baseFile");
open(V, $variantFile) || die("Cannot open $variantFile");

# Constants
my $PASS_COL = 6;

# Stats
my %baseCounts;
my %variantCounts;

my $numKept = 0;

my $lb;
my $lv;
while( ($lb = <B>) && ($lv = <V>) )
{
    chomp $lb;
    chomp $lv;

    if($lv =~ /^#/)
    {
        die unless($lb =~ /^#/);
        print $lv . "\n";
        next;
    }

    my @fb = split(' ', $lb);
    my @fv = split(' ', $lv);

    die("Error: Mismatching positions " . $fb[1] . " != " . $fv[1]) unless ($fb[1] == $fv[1]);
    $baseCounts{$fb[$PASS_COL]}++;
    $variantCounts{$fv[$PASS_COL]}++;

    # Output variants that PASS in the variant vcf and are no call in the base
    if($fb[6] eq "NoCall" && $fv[6] eq "PASS")
    {
        # Add the quality of the normal as an additional tag
        $fv[7] .= ";BQ=". $fb[5];
        print join("\t", @fv) . "\n";
    }
}

#
print STDERR "Filter stats for base file\n";
foreach my $k (keys %baseCounts)
{
    print STDERR "Num $k in base file: " . $baseCounts{$k} . "\n";
}

print STDERR "Filter stats for variant file\n";
foreach my $k (keys %variantCounts)
{
    print STDERR "Num $k in variant file: " . $variantCounts{$k} . "\n";
}
