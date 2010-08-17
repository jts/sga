#! /usr/bin/perl

use strict;
use Set::IntSpan;
use Bio::Perl;

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
my %matched_coords;
my %contig_lengths;

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
    my $left_clip = cigarLeftClipped($cigar);
    my $right_clip = cigarRightClipped($cigar);
    print "$cigar lc: $left_clip rc: $right_clip\n";
    my $start = $left_clip;
    my $end = $seq_length - $right_clip;
    my $coord = "$start-$end";
    $contig_lengths{$name} = $seq_length;

    if(!defined($matched_coords{$name}))
    {
        $matched_coords{$name} = new Set::IntSpan($coord);
    }
    else
    {
        print "s before: $matched_coords{$name}\n";
        $matched_coords{$name}->U($coord);
    }
    print "k: $name\n";
    print "p: $fields[2] $fields[3]\n";
    print "l: $seq_length\n";
    print "c: $coord\n";
    print "s after: $matched_coords{$name}\n";
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

# If any contigs did not have an alignment at all, add them here
my $file  = Bio::SeqIO->new(-file => $contig_file, '-format' => 'Fasta');
while(my $elem = $file->next_seq())
{
    my $name = $elem->id;
    my $seq = $elem->seq;
    if(!defined($best_match{$name}))
    {
        $best_match{$name}->{af} = 0;
        $best_match{$name}->{seq} = $seq;
        $best_match{$name}->{a_str} = "unaligned";
        $best_match{$name}->{cigar} = "";
        $best_match{$name}->{nm} = "N/A";
        $best_match{$name}->{matched} = 0;
        $contig_lengths{$name} = length($seq);
        $matched_coords{$name} = new Set::IntSpan;
        print $name . " unaligned\n";
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
        print "$k af: $af softclipped: $num_sc matched: $matched\n";
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

open(FA, ">$contig_file" . ".unmatched.fa");
open(SPLIT, ">$contig_file" . ".split.fa");
foreach my $k (keys %matched_coords)
{
    my $match = $matched_coords{$k};
    my $len = $contig_lengths{$k};
    my $full = "0-$len";
    my $full_set = new Set::IntSpan $full;
    my $unmatched = $full_set - $match;
    my $len_unmatched = $unmatched->cardinality();
    my $seq = $best_match{$k}->{seq};
    print "$k len: $len matched: $match unmatched: $unmatched ($len_unmatched)\n";
    print "$k LEN: $len_unmatched\n";

    if($len_unmatched > $len_filter)
    {
        print FA ">$k unmatched: $unmatched\n";
        print FA $best_match{$k}->{seq} . "\n";
    }

    my @spans = spans $unmatched;
    foreach my $s (@spans)
    {
        my $l = $s->[1] - $s->[0];
        if($l >= $len_filter)
        {
            my $sub = substr($seq, $s->[0], $l);
            print SPLIT ">$k:$s->[0]-$s->[1]\n";
            print SPLIT "$sub\n";
        }
    }
}
close(FA);
close(SPLIT);

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

sub cigarLeftClipped
{
    my($cs) = @_;

    # Parse the first code
    my ($len, $code, $rest) = ($cs =~ /(\d+)(\w)(.*)/);
    if($code eq 'S')
    {
        return $len;
    }
    return 0;
}

sub cigarRightClipped
{
    my($cs) = @_;
    my ($len, $code) = ($cs =~ /(\d+)(\w)$/);
    if($code eq 'S')
    {
        return $len;
    }
    return 0;
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
