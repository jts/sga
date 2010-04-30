#! /usr/bin/perl

use strict;

print "digraph G {";

while(<>)
{
	if(/VT/)
	{
		my @record = split;
		print $record[1] . " " . makeLabel($record[1]) . "\n;" 	
	}

	if(/ED/)
	{
		my @record = split;
		my $ol = $record[4] - $record[3] + 1;
		print qq(") . $record[1] . qq(") . " -> " . qq(") . $record[2] . qq(") . " " . makeLabel($ol) . "\n";
	}
}

print "}\n";

sub makeLabel
{
	my($l) = @_;
	my $str = qq([label=") . $l . qq("]);
	return $str;
}

