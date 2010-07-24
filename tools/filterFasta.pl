#! /usr/bin/perl
#
# Remove sequences from a fasta file that are less than a threshold length
# Usage: filterFasta.pl -n <N> <file>
#
use strict;
use Bio::Perl;
use Getopt::Long;

my $n = 0;
GetOptions("n=i" => \$n);

my $contigFile  = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'Fasta');
my $seq_out = Bio::SeqIO->new('-fh' => \*STDOUT,
                              '-format' => 'Fasta');
while(my $seq = $contigFile->next_seq())
{
    if($seq->length >= $n)
    {
        $seq_out->write_seq($seq);
    }
}
