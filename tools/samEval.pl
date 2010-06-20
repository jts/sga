#! /usr/bin/perl

while(<>)
{
	chomp;
	next if /^@/;
	my @fields = split(/\t/);
	my $cigar = $fields[5];
	my $m_len = substr($cigar, 0, length($cigar) - 1);
	my $s_len = length($fields[9]);
	my $correct = $m_len eq $s_len;
	
	if(!$correct)
	{
		#print ">$fields[0]\n";
		#print "$fields[9]\n";
		print "$m_len $s_len\n";
	}
}
