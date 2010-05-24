#! /usr/bin/perl

# Merge driver - make commands for merging BWTs together

use strict;

my @files = @ARGV;
my $n = scalar(@files);
my $finalName = "final";

# In optimize memory mode, we load the smaller of the two bwts into memory
# In optimize time mode, we load the larger of the two into memory
my $MODE_OPT_MEMORY = 0;
my $MODE_OPT_TIME = 1;
my $mode = $MODE_OPT_MEMORY;

my $finalParam = "";
if($n == 2)
{
    $finalParam = "-p $finalName";
}

# Sort the input files by size, smallest to largest
my @sorted = sort { getFilesize($a) <=> getFilesize($b) } @files;
#print join(" ", @sorted) . "\n";
# Merge the largest file with the smallest, etc
my $half = $n / 2;
my $i = 0;
my $j = $n - 1;
while($i < $j)
{
    print makeMergeLine($sorted[$i], $sorted[$j]);
    ++$i;
    --$j;
}

sub makeMergeLine
{
    my($f1, $f2) = @_;

    my $larger;
    my $smaller;
    if(-s $f1 > -s $f2)
    {
        $larger = $f1;
        $smaller = $f2;
    }
    else
    {
        $larger = $f2;
        $smaller = $f1;
    }

    if($mode == $MODE_OPT_MEMORY)
    {
        return qw($SGA) . " merge -r $finalParam $larger $smaller\n";
    }
    else
    {
        return qw($SGA) . " merge -r $finalParam $smaller $larger\n";
    }
}

sub getFilesize
{
    my($a) = @_;
    return -s $a;
}
