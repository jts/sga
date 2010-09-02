#! /usr/bin/perl

use strict;

print "digraph G {\n";

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
        my $s1 = $record[3];
        my $s2 = $record[6];
		print quoteStr($record[1]) . " -> " . quoteStr($record[2]) . " " . makeLabel($ol, getColor($s1)) . ";\n";
        print quoteStr($record[2]) . " -> " . quoteStr($record[1]) . " " . makeLabel($ol, getColor($s2)) . ";\n";
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
	my($l, $color) = @_;
    my $c_str;
    if($color ne "")
    {
        $c_str = "color=\"$color\"";
    }
	my $str = qq([label=") . $l . qq(" ) . $c_str. qq(]);
	return $str;
}

sub getColor
{
    my($s) = @_;
    if($s == 0)
    {
        return "red";
    }
    else
    {
        return "black";
    }
}

