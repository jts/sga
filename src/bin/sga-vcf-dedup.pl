#! /usr/bin/perl
use strict;
# Remove duplicated entries from a VCF file

my %seen_hash;

while(<>) {
    # Print header lines
    if(/^#/) {
        print;
        next;
    }

    chomp;
    my @f = split('\t');

    my $key = sprintf("%s.%d.%s.%s", $f[0], $f[1], $f[3], $f[4]);
    if(!defined($seen_hash{$key})) {
        print $_ . "\n";
        $seen_hash{$key} = 1;
    }
}
