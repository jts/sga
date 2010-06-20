#! /usr/bin/perl

use strict;

while(<>)
{
	if(/VT/)
	{
		my @record = split;
        print ">" . $record[1] . "\n";
        print $record[2] . "\n";
	}
}
