#! /usr/bin/perl
#
# Convert an abyss DistanceEst file to a graphviz dot file
# Usage: de2dot.pl <contig file> <de file>
#
use strict;
use Bio::Perl;
use POSIX qw(ceil);
use FindBin;
use lib $FindBin::Bin;
use ParseDE;

my $contigFile = $ARGV[0];
my $deFile = $ARGV[1];

my $contigFile  = Bio::SeqIO->new(-file => $contigFile, '-format' => 'Fasta');
# read contig lengths
my %contigLength;
while(my $seq = $contigFile->next_seq())
{
    print STDERR $seq->id . " " . $seq->length . "\n";
    $contigLength{$seq->id} = $seq->length;
}


open my $deFH, $deFile;

print "graph {\n";
# Skip the first line, it is a header
my $skip = <$deFH>;
while(my $line = <$deFH>)
{
    my ($baseCtg, $links) = ParseDE::parseDERecord($line);
    if(scalar(@{$links}) > 0)
    {
        my $len = $contigLength{$baseCtg};
        #print STDERR "$baseCtg $len\n";
        my $size = 0.5;
        my $label = makeLabel($baseCtg);
        if($len < 500)
        {
            $size = 0.1;
            $label = makeLabel("<500");
        }
        elsif($len < 1000)
        {
            $size = 0.5;
            $label = makeLabel("<1kb");
        }
        elsif($len < 5000)
        {
            $size = 1.0;
        }
        else
        {
            $size = 2;
        }
        my $nodeText = sprintf("%s [fixedsize=1 label=%s width=%f]\n", makeLabel($baseCtg), $label, $size);
        print $nodeText;
    }

    foreach my $l (@{$links})
    {
        if($baseCtg le $l->{ctg})
        {
            my $colText = "color=" . makeLabel(($l->{dir} == 0 ? "red" : "black"));
            my $d = $l->{dist};
            my $t = 3;
            my $sd = $l->{sd};
            my $labText = makeLabel("(" . ($d - ceil($t*$sd)) . ":" . ($d + ceil($t*$sd)) . ")");
            print makeLabel($baseCtg) . " -- " . makeLabel($l->{ctg}) . qq( [label=) . $labText . " $colText];\n"
        }
    }
}
print "}\n";

sub makeLabel
{
    my($str) = @_;
    return qq(") . $str . qq(");
}
