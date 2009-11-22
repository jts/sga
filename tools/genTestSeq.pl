#! /usr/bin/perl

use strict;
use Getopt::Long;

my $l = 1000;

GetOptions("length=i" => \$l);
my @bases = ('A', 'C', 'G', 'T');

print ">generated_seq $l\n";
for(my $i = 0; $i < $l; ++$i)
{
	print randomBase();
}
print "\n";

sub randomBase
{
	my $i = int(rand(4));
	return $bases[$i];
}
