#! /usr/bin/perl
# Cat a set of fasta/fastq files together and build an index of read index -> original file

use strict;
use Getopt::Long;
use Bio::Perl;
use File::Basename;

my $prefix = "final";
my $sequence_filename = "$prefix.fa";
my $index_filename = "$prefix.popidx";

# Open filehandles for the output sequences and the output index
open(OUT, ">$sequence_filename") || die("Cannot open $sequence_filename");
open(IDX, ">$index_filename") || die("Cannot open $index_filename");

# Track position data
my $current_index = 0;
my $current_start_index = 0;
my $current_name = "";

# Iterate over every file
foreach my $in_filename (@ARGV) {

    # Write out the previous file's data to the index, if any
    if($current_name ne "") {
        print IDX join("\t", ($current_start_index, $current_index - 1, $current_name)) . "\n";
    }
    
    # Reset position data
    $current_name = basename($in_filename, (".fa", ".fastq"));
    $current_start_index = $current_index;

    print "Processing $in_filename with $current_name\n";

    # Iterate over all reads in the file
    open(IN, $in_filename) || die("Cannot open $in_filename");
    while(my $line = <IN>) {
        chomp $line;
        my ($header) = split(' ', $line);

        my $record = "";
        if($header =~ /^>/) {
            # parse fasta, assume 1 line per sequence
            $record = $header . "\n" . <IN>;
        } elsif($header =~ /^@/) {
            # parse fastq
            $record = $header . "\n";
            $record .= <IN>;
            $record .= <IN>;
            $record .= <IN>;
        } else {
            die("Unexpected format\n");
        }

        $current_index += 1;
        print OUT $record;
    }
}

# Write the last element of the index
print IDX join("\t", ($current_start_index, $current_index - 1, $current_name)) . "\n";

close(IDX);
close(OUT);
