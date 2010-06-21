#! /usr/bin/perl

use strict;
use Getopt::Long;

my $error_rate = 0.01;
my @bases = ('A', 'C', 'G', 'T');
my $errors_created = 0;
my $profile_file = "";

GetOptions("error-rate=f" => \$error_rate,
           "profile-file=s" => \$profile_file);

my @error_profile;
my $profile_len = 0;
if($profile_file ne "")
{
    open(FILE, $profile_file);
    while(<FILE>)
    {
        my @fields = split;
        push @error_profile, $fields[3];
    }
    $profile_len = scalar(@error_profile);

    # Calculate the error rate implied by the profile
    my $sum_error = calculateTotalErrorRate();
    print STDERR "Input profile error rate: $sum_error\n";

    my $scale = $error_rate / $sum_error;

    # Uniformly rescale the profile probabilities so the total error rate is $error_rate
    for(my $i = 0; $i < $profile_len; ++$i)
    {
        $error_profile[$i] *= $scale
    }

    $sum_error = calculateTotalErrorRate();
    print STDERR "Scaled profile error rate: $sum_error\n";
}

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
        my $error_p = $error_rate;
        if($profile_len > 0)
        {
            die("profile length is less than read length") if $i > $profile_len;
            $error_p = $error_profile[$i];
        }

		if(rand() < $error_p)
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

sub calculateTotalErrorRate
{
    # Calculate the error rate implied by the profile
    my $sum_error = 0;
    for(my $i = 0; $i < $profile_len; ++$i)
    {
        $sum_error += (1.0 / $profile_len) * $error_profile[$i];
    }
    return $sum_error;
}
