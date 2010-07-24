#! /usr/bin/perl

use strict;

my $contig_file = $ARGV[0];
my $ref_prefix = $ARGV[1]; #"~/work/devel/sga/data/ecoli_k12.fa";
my $bwa = "/software/solexa/bin/aligners/bwa/current/bwa";
my $samtools = "/software/solexa/bin/aligners/samtools/current/samtools";
my $sam_file = $contig_file . ".sam";
my $sam_md_file = $contig_file. ".md.sam";
system("$bwa bwasw $ref_prefix $contig_file > $sam_file");
system("$samtools calmd -S $sam_file $ref_prefix > $sam_md_file");

# Parse the results
open(SAM, $sam_md_file) || die("Cannot open $sam_md_file");

my $correct_threshold = 0.95;
my $num_wrong = 0;
my $num_right = 0;
my $num_base_mismatch = 0;
my $num_base_softclip = 0;
my $sum_misassembled = 0;
my $len_filter = 500;
my $missing_threshold = 100; # if more than this value is missing from the alignment the contig is incorrect
my %best_match;

open(OUT, ">divergent.fa");
while(<SAM>)
{
    # skip header
    next if(/^@/);

    my @fields = split(/\t/);

    my $name = $fields[0];
    my $cigar = $fields[5];
    my $seq = $fields[9];
    my $seq_length = length($seq);
    next if $seq_length < $len_filter;
    my $matched = parseCigar($cigar);
    my $align_frac = ($matched / $seq_length);
    my $nm = findField(\@fields, "NM");

    #print "$name $matched $seq_length " . ($matched / $seq_length) . "\n";
    if(!defined($best_match{$name}) || $best_match{$name}->{af} < $align_frac)
    {
        $best_match{$name}->{af} = $align_frac;
        $best_match{$name}->{seq} = $seq;
        $best_match{$name}->{a_str} = "$fields[2] $fields[3]";
        $best_match{$name}->{cigar} = $cigar;
        $best_match{$name}->{nm} = $nm;
        $best_match{$name}->{matched} = $matched;
    }
}

foreach my $k (keys %best_match)
{
    my $af = $best_match{$k}->{af};
    my $len =  length($best_match{$k}->{seq});
    my $num_sc = countSoftclip($best_match{$k}->{cigar});
    my $matched = $best_match{$k}->{matched};
    if($len >= $len_filter)
    {
        print "af: $af softclipped: $num_sc matched: $matched\n";
        if(($len - $matched) >= $missing_threshold)
        {
            print "$k (len: " . $len . ") is potentially misassembled, best match $af\n";
            print OUT ">" . join(" ", ($k, length($best_match{$k}->{seq}), $best_match{$k}->{a_str}, $af)) . "\n" . $best_match{$k}->{seq} . "\n";
            ++$num_wrong;
            $sum_misassembled += $len;
        }
        else
        {
            ++$num_right;
            $num_base_mismatch += countMismatch($best_match{$k}->{nm});
            $num_base_softclip += $num_sc;
        }
    }
}

close(OUT);
print "Num wrong: $num_wrong\n";
print "Num right: $num_right\n";
print "Sum misassembled: $sum_misassembled\n";
print "Num mismatches: $num_base_mismatch\n";
print "Num softclip: $num_base_softclip\n";
print "Total mismatches: " . ($num_base_mismatch + $num_base_softclip) . "\n";

sub parseCigar
{
    my($cs) = @_;
    my $done = 0;
    my $amount_matched = 0;
    while(!$done)
    {
        my ($len, $code, $rest) = ($cs =~ /(\d+)(\w)(.*)/);
        if($code eq "M")
        {
            $amount_matched += $len;
        }
        $done = 1 if $rest eq "";
        $cs = $rest;
    }
    print "matched: $amount_matched\n";
    return $amount_matched;
}

# for contigs we calculate the number of mismatches as the value in the NM field plus
# the number of bases that were softclipped off the alignment
sub countMismatch
{
    my($nm) = @_;

    # Count the number of soft-clipped bases
    my ($v) = ($nm =~ /NM:i:(\d+)/);
    return $v;
}

# Count the number of bases that were softclipped off the end of the contig
sub countSoftclip
{
    my($cigar) = @_;

    # Count the number of soft-clipped bases
    my $soft_count = countCigarCode($cigar, 'S');
    return $soft_count;
}

sub countCigarCode
{
    my($cs, $expected_code) = @_;
    my $done = 0;
    my $count = 0;
    my $amount_matched = 0;
    while(!$done)
    {
        my ($len, $code, $rest) = ($cs =~ /(\d+)(\w)(.*)/);
        if($code eq $expected_code)
        {
            $count += $len;
        }
        $done = 1 if $rest eq "";
        $cs = $rest;
    }
    return $count;
}

# return the field with the given tag
sub findField
{
    my($fields, $tag) = @_;

    foreach my $f (@{$fields})
    {
        if($f =~ /$tag/)
        {
            return $f;
        }
    }
    return "";
}
