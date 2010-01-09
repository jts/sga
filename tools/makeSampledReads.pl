#! /usr/bin/perl
#

use strict;
use Getopt::Long;

my $file = "";
my $rl = 36;
my $pe_mean = 200;
my $pe_sd = 0.1;
my $coverage = 30;
my $bHelp = 0;
my $bVerbose = 0;
my $bSameStrand = 0;
my $bTrackPos = 0;
my $bSingleEnd = 0;

GetOptions("length=i" => \$rl,
			"pe_mean=i" => \$pe_mean,
			"pe_sd=f" => \$pe_sd,
			"coverage=i" => \$coverage,
			"same-strand" => \$bSameStrand,
			"single-end" => \$bSingleEnd,
			"track" => \$bTrackPos,
			"help" => \$bHelp,
			"verbose" => \$bVerbose);

if($bHelp)
{
	print STDERR "makeSampledReads.pl [OPTIONS] <FILE>\n";
	print STDERR "Sample reads from FILE\n";
	print STDERR "Options:\n";
	print STDERR "--coverage=i        The depth of coverage per base in the genome (for example: 30X would be --cov 30)\n";
	print STDERR "--length=i          The length of the reads\n";
	print STDERR "--same_strand       Force both reads to be from the same strand\n";
	print STDERR "--single-end        Do not make paired reads\n";
	print STDERR "--track             Generate unpaired reads and use the position the read was sampled from as the ID\n";
	print STDERR "--pe_mean=i         The mean fragment size for paired-reads\n";
	print STDERR "--pe_sd=f           The standard deviation of the fragment distribution\n";
	print STDERR "                       if this is <= 1.0 it is interpreted to be a fraction\n";
	print STDERR "                       of the mean, otherwise it is the actual SD in base pairs\n";
	exit(1);
}

if($pe_sd <= 1.0)
{
	$pe_sd = $pe_sd * $pe_mean;
}

if($bVerbose)
{
	print STDERR "Coverage: $coverage\n";
	print STDERR "Fragment size mean: $pe_mean\n";
	print STDERR "Fragment size sd: $pe_sd\n";
	print STDERR "Read length: $rl\n";
	print STDERR "Tracking: $bTrackPos\n";
}

my $total = 0;
my $buffer = "";

# slurp in the file
while(<>)
{
	chomp;
	if(/^>/)
	{
		if(length($buffer) > 0)
		{
			if($bTrackPos || $bSingleEnd)
			{
				outputSEReads(\$buffer);
			}
			else
			{
				outputPEReads(\$buffer);
			}
			$buffer = "";
		}
	}
	else
	{
		$buffer .= $_;
	}
}

if($bTrackPos || $bSingleEnd)
{
	outputSEReads(\$buffer);
}
else
{
	outputPEReads(\$buffer);
}

sub outputPEReads
{

	die("PE mode is deprecated - Normal Distribution has been hacked out");

	my($buffer) = @_;
	my $gl = length($$buffer);
	my $num_reads = $gl * $coverage /  (2 * $rl);

	for(my $i = 0; $i < $num_reads; ++$i)
	{
		my $start_1 = int(rand($gl));
		my $start_2 = 200 + $start_1 - $rl;
		my $end_1 = $start_1 + $rl - 1;
		my $end_2 = $start_2 + $rl - 1;
		next if $end_2 >= $gl;

		my $seq_1 = getRead($buffer, $start_1, $end_1); 
		my $seq_2 = getRead($buffer, $start_2, $end_2);
		$seq_2 = rc($seq_2) unless $bSameStrand;
		next if($seq_1 =~ /N/ || $seq_2 =~ /N/);
		print join("\n", (">$total/A", $seq_1, ">$total/B", $seq_2)) . "\n";
		++$total;
	}
}

sub outputSEReads
{
	my($buffer) = @_;
	my $gl = length($$buffer);
	my $num_reads = $gl * $coverage /  $rl;
	
	for(my $i = 0; $i < $num_reads; ++$i)
	{
		my $start = int(rand($gl));
		my $end = $start + $rl - 1;
		next if $end >= $gl;

		my $seq = getRead($buffer, $start, $end); 
		my $name = $bTrackPos ? "$total:$start" : $total;
		$seq = rc($seq) unless ($bSameStrand || rand() < 0.5);
		next if($seq =~ /N/ || length($seq) < $rl);
		print join("\n", (">$name", $seq)) . "\n";
		++$total;
	}

}

sub getRead
{
	my($buffer, $start, $end) = @_;
	#return uc(join("",@$buffer[$start..$end]));
	return substr($$buffer, $start, $end - $start + 1);
}

sub rc
{
	my($seq) = @_;
	$seq =~ tr/ACGT/TGCA/;
	return reverse $seq;
}

