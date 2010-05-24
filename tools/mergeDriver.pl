#! /usr/bin/perl

# Merge driver - make commands for merging BWTs together

use strict;

my @files = @ARGV;
my $n = scalar(@files);
my $finalName = "final";

my $finalParam = "";
if($n == 2)
{
    $finalParam = "-p $finalName";
}

for(my $i = 0; $i < $n; $i += 2)
{
    next if($i == $n - 1);
    print qw($SGA) . " merge -r $finalParam $files[$i] $files[$i + 1]\n";
}
