#! /usr/bin/perl
#
# Convert an abyss DistanceEst file to a bambus XML evidence file
# Usage: de2evidence.pl <contigs fasta> <de file>
#
use strict;
use FindBin;
use lib $FindBin::Bin;
use ParseDE;

my $NUM_OPEN_TAGS = 0;

my $contigFile = $ARGV[0];
my $deFile = $ARGV[1];

print qq(<?xml version="1.0" ?>\n);
startTag("EVIDENCE", qq(ID="1" DATE="" PROJECT="" PARAMETERS=""));

# Build the contig records
open(CTGS, $contigFile);
my $currID = "";
my $currLen = 0;
while(my $line = <CTGS>)
{
    chomp $line;
    if($line =~ /^>/)
    {
        if($currID ne "")
        {
            writeContigRecord($currID, $currLen);
        }
        ($currID) = ($line =~ />(\S+).*/);
        $currLen = 0;
    }
    else
    {
        $currLen += length($line);
    }
}
writeContigRecord($currID, $currLen);
close(CTGS);

# Build the link records
open(DE, $deFile);

# Skip the first line, it is a header
my $skip = <DE>;
while(my $line = <DE>)
{
    my ($baseCtg, $linkArray) = ParseDE::parseDERecord($line);
    foreach my $linkRec (@{$linkArray})
    {
        #print "L $base $id $dir $dist " . ($base < $id) . "\n";
        # We only want to output one link per pair but the DE file
        # contains a record for both contigs. Only output the record
        # for the lexographically lower of the two
        if($baseCtg lt $linkRec->{ctg})
        {
            writeLinkRecord($baseCtg, $linkRec->{ctg}, $linkRec->{strand}, $linkRec->{dist});
        }
    }
}
endTag("EVIDENCE");

sub writeLinkRecord
{
    my ($id1, $id2, $orientation, $dist) = @_;
    my $lid = $id1 . "_" . $id2;
    my $type = "DistanceEst";
    my $info = qq(ID=") . $lid . qq(" SIZE=") . $dist . qq(" TYPE=") . $type . qq(");
    startTag("LINK", $info);
    my $cline1 = qq(ID=") . $id1 . qq(" ORI=") . "BE" . qq(");
    my $cline2 = qq(ID=") . $id2 . qq(" ORI=") . (($orientation eq "+") ? "BE" : "EB") . qq(");
    simpleTag("CONTIG", $cline1);
    simpleTag("CONTIG", $cline2);
    endTag("LINK");
}

sub writeContigRecord
{
    my($id, $len) = @_;
    my $info = qq(ID=") . $id . qq(" NAME=") . $id . qq(" LEN=") . $len . qq(");
    startTag("CONTIG", $info);
    endTag("CONTIG");
}

sub simpleTag
{
    my($tag, $meta) = @_;
    print "\t" x $NUM_OPEN_TAGS . "<$tag $meta/>\n";

}

sub startTag
{
    my($tag, $meta) = @_;
    print "\t" x $NUM_OPEN_TAGS . "<$tag $meta>\n";
    $NUM_OPEN_TAGS++;
}

sub endTag
{
    my($tag) = @_;
    $NUM_OPEN_TAGS--;
    print "\t" x $NUM_OPEN_TAGS . "</$tag>\n";
}
