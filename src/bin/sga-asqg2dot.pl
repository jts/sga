#! /usr/bin/perl

use strict;
use Getopt::Long;

my $attrFile = "";

# Params
my $bDirected = 1;
my $bUseLabels = 0;
my $bUseShapes = 0;
my $defaultShape = "circle";
my $noVertexNames = 0;

GetOptions("suppress-labels" => \$noVertexNames);

my %vertexColors;
my %vertexShapes;
my %vertexSizes;

my $graphAttr;
if($bDirected)
{
    $graphAttr = "digraph";
}
else
{
    $graphAttr = "graph";
}

print "$graphAttr G {\n";

while(<>)
{
    if(/VT/)
    {
        my @record = split;
        my $id = $record[1];
        my $len = length($record[2]);
        my @attributes;
        my $label = $noVertexNames ? quoteStr("") : quoteStr($id);
        push @attributes, makeAttribute("label", $label);
        my $attrStr = attributes2string(@attributes);
        print quoteStr($id) . " " . $attrStr . ";\n"    
    }

    if(/ED/)
    {
        my @record = split;
        my $ol = $record[4] - $record[3] + 1;
        my $s1 = $record[3];
        my $s2 = $record[6];
        my @attributes; 

        my $edgeStr = quoteStr($record[1]) . " " . getArrow() . " " . quoteStr($record[2]);
        
        push @attributes, makeAttribute("color", getEdgeColor($s1));
        my $attrStr = attributes2string(@attributes);

        print $edgeStr . " " . $attrStr . ";\n";
        if($bDirected)
        {
            print quoteStr($record[2]) . " " . getArrow() . " " . quoteStr($record[1]) . " " . makeLabel($ol, getEdgeColor($s2)) . ";\n";
        }
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

sub makeAttribute
{
    my($name, $value) = @_;
    return $name . "=" . $value;
}

sub attributes2string
{
    return "[" . join(" ", @_) . "]";
}

sub getArrow
{
    if($bDirected)
    {
        return "->";
    }
    else
    {
        return "--"
    }
}

sub getEdgeColor
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

sub getVertexColor
{
    my($id) = @_;
    return "white" if !defined($vertexColors{$id});
    return $vertexColors{$id};
}

sub getVertexShape
{
    my($id) = @_;
    return "circle" if !$bUseShapes || !defined($vertexShapes{$id});
    return $vertexShapes{$id};
}

sub getVertexSize
{
    my($id) = @_;
    return 1 if !defined($vertexSizes{$id});
    return $vertexSizes{$id};
}
