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
my $seq_out = Bio::SeqIO->new('-file' => ">$sequence_filename", '-format' => 'Fasta');
open(IDX, ">$index_filename") || die("Cannot open $index_filename");

# Track position data
my $current_index = 0;
my $current_start_index = 0;
my $current_name = "";

# Iterate over every file
foreach my $in_filename (@ARGV) {

    # Write out the previous file's data to the index, if any
    if($current_name ne "") {
        print IDX join("\t", ($current_start_index, $current_index - 1, $current_name));
    }
    
    # Reset position data
    $current_name = basename($in_filename, (".fa", ".fastq"));
    $current_start_index = $current_index;

    print "Processing $in_filename with $current_name\n";

    # Iterate over all reads in the file
    my $in  = Bio::SeqIO->new(-file => $in_filename);
    while(my $seq = $in->next_seq()) {
        $seq_out->write_seq($seq);
        $current_index += 1;
    }
}

# Write the last element of the index
print IDX join("\t", ($current_start_index, $current_index - 1, $current_name));

close(IDX);
$seq_out->close();
