#! /usr/bin/perl

use strict;

while(<>)
{
	chomp;
	next if /Drawing/;
	my @fields = split;
	next unless scalar(@fields) == 5;
	print ">$fields[4]\n";
	print "$fields[0]\n";

}
