#! /usr/bin/perl

use strict;

print "graph G {\n";

while(<>)
{
	if(/VT/)
	{
		my @record = split;
		print quoteStr($record[1]) . " " . makeLabel($record[1]) . ";\n" 	
	}

	if(/ED/)
	{
		my @record = split;
		my $ol = $record[4] - $record[3] + 1;
		print quoteStr($record[1]) . " -- " . quoteStr($record[2]) . " " . makeLabel($ol) . ";\n";
	}
}

print "}\n";

sub quoteStr
{
	my ($s) = @_;
	return qq(") . $s . qq(");
}

sub makeLabel
{
	my($l) = @_;
	my $str = qq([label=") . $l . qq("]);
	return $str;
}

