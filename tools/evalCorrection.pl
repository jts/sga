#! /usr/bin/perl

use strict;

# The reference file contains the reads without errors, the truth
my $referenceReadsFile = $ARGV[0];
my $uncorrectReadsFile = $ARGV[1];
my $correctedReadsFile = $ARGV[2];

open(REF, $referenceReadsFile);
open(UNC, $uncorrectReadsFile);
open(COR, $correctedReadsFile);

my $refSeq;
my $uncSeq;
my $corSeq;

my $total_corrected = 0; # num bases successfully corrected
my $total_corrupted = 0; # num bases incorrectly changed
my $total_wrong_before = 0; # num bases that were wrong pre-correct
my $total_wrong_after = 0; # num bases that are wrong after-correction
my $total_seq = 0;
my %distBefore;
my %distAfter;

while(($refSeq = <REF>) && ($uncSeq= <UNC>) && ($corSeq = <COR>))
{
    chomp $refSeq;
    chomp $uncSeq;
    chomp $corSeq;
    if($refSeq =~ />/)
    {
        my $id1 = getID($refSeq);
        my $id2 = getID($uncSeq);
        my $id3 = getID($corSeq);
        die("Mismatched id lines") if $id1 ne $id2 || $id1 ne $id3;
    }
    else
    {
        my ($a, $b, $c, $d) = threeWayCompare($refSeq, $uncSeq, $corSeq);
        $total_corrected += $a;
        $total_corrupted += $b;
        $total_wrong_before += $c;
        $total_wrong_after += $d;
        $total_seq += length($refSeq);
        $distBefore{$c}++;
        $distAfter{$d}++;
    }
}

print "Num corrected: $total_corrected\n";
print "Num corrupted: $total_corrupted\n";
printf "Total errors before correction: %d (%1.3f%)\n", ($total_wrong_before, $total_wrong_before / $total_seq);
printf "Total errors after correction: %d (%1.3f%)\n", ($total_wrong_after, $total_wrong_after / $total_seq);

for(my $i = 0; $i < 10; ++$i)
{
    my $vb = (not defined $distBefore{$i}) ? 0 : $distBefore{$i};
    my $va = (not defined $distAfter{$i}) ? 0 : $distAfter{$i};

    print join("\t", ($i, $vb, $va)) . "\n";
}

sub threeWayCompare
{
    my($r, $u, $c) = @_;
    my $l = length($r);
    my $num_corrected = 0; # num bases successfully corrected
    my $num_corrupted = 0; # num bases incorrectly changed
    my $num_original = 0; # num bases that were wrong pre-correct
    my $num_after = 0; # num bases that are wrong after-correction
    my $nc = 0;

    for(my $i = 0; $i < $l; ++$i)
    {
        my $unc_correct = substr($r, $i, 1) eq substr($u, $i, 1);
        my $cor_correct = substr($r, $i, 1) eq substr($c, $i, 1);
       
        ++$num_corrected if(!$unc_correct && $cor_correct);
        ++$num_corrupted if($unc_correct && !$cor_correct);
        ++$num_original if(!$unc_correct);
        ++$num_after if(!$cor_correct);

    }
    return ($num_corrected, $num_corrupted, $num_original, $num_after);

}

sub classifyDifferences
{
    my($a, $b) = @_;
    my $la = length($a);
    my $nd = 0;
    for(my $i = 0; $i < $la; ++$i)
    {
        ++$nd if(substr($a, $i, 1) ne substr($b, $i, 1));
    }
    return $nd;
}

sub getID
{
    my($l) = @_;
    my ($id) =~ (/>(.*) /);
    return $id;
}
