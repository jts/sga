#! /bin/bash

#
# Example assembly of 100bp C. elegans data set
#

# We assume the data is downloaded from the SRA and converted to fastq files
# Set IN1 and IN2 to be the paths to the data on your filesystem
IN1=SRR065390_1.fastq
IN2=SRR065390_2.fastq

# Parameters
SGA_BIN=sga

# The number of threads to use
CPU=8

# To save memory, we index $D reads at a time then merge the indices together
D=4000000

# First, preprocess the data to remove ambiguous basecalls
$SGA_BIN preprocess --pe-mode 1 -o SRR065390.fastq $IN1 $IN2

# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
$SGA_BIN index -d $D -t $CPU --no-reverse SRR065390.fastq

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
$SGA_BIN correct -k 41 --learn -t $CPU -o reads.ec.k41.fastq SRR065390.fastq

# Index the corrected data.
$SGA_BIN index -d $D -t $CPU reads.ec.k41.fastq

# Remove exact-match duplicates and reads with low-frequency k-mers
$SGA_BIN filter -x 2 -t $CPU reads.ec.k41.fastq

# Merge simple, unbranched chains of vertices
$SGA_BIN fm-merge -m 45 -t $CPU -o merged.m45.k41.fa reads.ec.k41.filter.pass.fa

# Build an index of the merged sequences
$SGA_BIN index -d 1000000 -t $CPU merged.m45.k41.fa

$SGA_BIN rmdup -t $CPU merged.m45.k41.fa

# Compute the structure of the string graph
$SGA_BIN overlap -m 45 -t $CPU merged.m45.k41.rmdup.fa

# Perform the contig assembly
$SGA_BIN assemble -m 75 -o assemble.m75 merged.m45.k41.rmdup.asqg.gz
