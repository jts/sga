#! /usr/bin/perl

use strict;
use Bio::Perl;

if(@ARGV < 1)
{
    print "usage: breakScaffolds.pl <scaffolds.fa>\n";
    print "Break scaffolds into a set of contigs and generate a layout .scaf file\n";
}

my $scaffoldFile = $ARGV[0];

my $scaffolds  = Bio::SeqIO->new(-file => $scaffoldFile, '-format' => 'Fasta');

my $outfile = "scaffoldContigs.fa";
my $seqout = Bio::SeqIO->new(-file   => ">$outfile",
                             -format => "fasta");
open(SCAF, ">scaffoldContigs.scaf");

my $numCtgs = 0;
my $numGaps = 0;
my $numScaffolds = 0;
while(my $scaffold = $scaffolds->next_seq())
{
    ++$numScaffolds;
    my $currGap = 0;
    my $currIdx = 0;
    
    my $seq = $scaffold->seq;
    my @contigs = split('(N+)', $seq);
    foreach my $c (@contigs)
    {
        if($c =~ /N/)
        {
            $currGap = length($c);
            ++$numGaps;
        }
        else
        {
            my $id = $scaffold->display_id . ".$currIdx";
            my $outobj = Bio::PrimarySeq->new ( -seq => $c,
                                                -id  => $id );
            
            if($currIdx == 0)
            {
                print SCAF "$id";
            }
            else
            {
                print SCAF "\t" . join(",", $id, $currGap, 0, 1, 0, "R");
            }
            $seqout->write_seq($outobj); 

            ++$currIdx;
            ++$numCtgs;
        }
    }
    print SCAF "\n";
}

print "Broke $numScaffolds scaffolds into $numCtgs contigs and $numGaps gaps\n";
