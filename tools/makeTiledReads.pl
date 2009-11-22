#! /usr/bin/perl
#

use strict;
use Getopt::Long;

my $file = "";
my $rl = 36;
my $overlap = -1;
my $bHelp = 0;
my $bAltRC = 0;
GetOptions("len=i" => \$rl,
			"overlap=i" => \$overlap,
			"alternate_strand" => \$bAltRC,
			"help" => \$bHelp);

if($bHelp)
{
	print STDERR "makeTiledReads.pl --len n --overlap n <FILE>\n";
	exit(1);
}

my $step = 1;
if($overlap >= 0)
{
	$step = $rl - $overlap;
}

my $total = 0;
my @buffer;
# slurp in the file
while(<>)
{
	chomp;
	if(/^>/)
	{
		if(scalar(@buffer) > 0)
		{
			outputReads(\@buffer);
			@buffer = ();
		}
	}
	else
	{
		push @buffer, split('');
	}
}

outputReads(\@buffer);


sub outputReads
{
	my($buffer) = @_;
	my $gl = scalar(@$buffer);
	for(my $i = 0; $i < ($gl - $rl + 1); $i += $step)
	{
		my $end = $i + $rl - 1;
		my $seq = uc(join("",@$buffer[$i..$end]));
		$seq = rc($seq) if($bAltRC && $total % 2 == 1);
		print ">$total:$i\n" . $seq . "\n";
		++$total;
	}
}

sub rc
{
	my($seq) = @_;
	$seq =~ tr/ACGT/TGCA/;
	return reverse $seq;
}
