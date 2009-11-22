#! /usr/bin/perl

use strict;
use Getopt::Long;

my $error_rate = 0.01;
my @bases = ('A', 'C', 'G', 'T');
my $errors_created = 0;
GetOptions("error-rate=f" => \$error_rate);

my $header ="";
while(<>)
{
	if(/^>/)
	{
		$header = $_;
		chomp $header;
	}
	else
	{
		my $seq = $_;
		chomp $seq;
		my ($p, $n) = permuteSeq($seq);
		print $header . " " . "$n\n";
		print $p . "\n";
	}
}

print STDERR "Created $errors_created errors\n";
sub permuteSeq
{
	my($s) = @_;
	my @arr = split('', $s);
	my $out = "";
	my $n = 0;
	for(my $i = 0; $i < @arr; $i++)
	{
		if(rand() < $error_rate)
		{
			$out .= randBaseNot($arr[$i]);
			++$n;
		}
		else
		{
			$out .= $arr[$i];
		}
	}
	$errors_created += $n;
	return ($out,$n);
}

sub randBaseNot
{
	my($b) = @_;
	my $o;
	while(1)
	{
		$o = randBase();
		return $o if $o ne $b;
	}
}

sub randBase
{
	my $i = int(rand(4));
	return $bases[$i];
}
