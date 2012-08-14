#! /bin/bash -x

#
# Example assembly of 100bp C. elegans data set. The only argument
# this script takes is the overlap length used for the final contig assembly.
#

# We assume the data is downloaded from the SRA and converted to fastq files
# Set IN1 and IN2 to be the paths to the data on your filesystem
IN1=SRR065390_1.fastq
IN2=SRR065390_2.fastq

# Parameters
SGA_BIN=sga-0.9.31

# Overlap parameter used for the final assembly. This is the only argument
# to the script
OL=75

# The number of threads to use
CPU=8

# To save memory, we index $D reads at a time then merge the indices together
D=4000000

# Correction k-mer value
CK=41

# The minimum k-mer coverage for the filter step. Each 27-mer
# in the reads must be seen at least this many times
COV_FILTER=2

# Overlap parameter used for FM-merge. This value must be no greater than the minimum
# overlap value you wish to try for the assembly step.
MOL=55

# Parameter for the small repeat resolution algorithm
R=10

# The number of pairs required to link two contigs into a scaffold
MIN_PAIRS=5

# The minimum length of contigs to include in a scaffold
MIN_LENGTH=200

# Distance estimate tolerance when resolving scaffold sequences
SCAFFOLD_TOLERANCE=1

# Turn off collapsing bubbles around indels
MAX_GAP_DIFF=0

# First, preprocess the data to remove ambiguous basecalls
$SGA_BIN preprocess --pe-mode 1 -o SRR065390.fastq $IN1 $IN2

#
# Error correction
#
# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
$SGA_BIN index -a ropebwt -t $CPU --no-reverse SRR065390.fastq

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
$SGA_BIN correct -k $CK --discard --learn -t $CPU -o reads.ec.k$CK.fastq SRR065390.fastq

#
# Contig assembly
#

# Index the corrected data.
$SGA_BIN index -a ropebwt -t $CPU reads.ec.k$CK.fastq

# Remove exact-match duplicates and reads with low-frequency k-mers
$SGA_BIN filter -x $COV_FILTER -t $CPU --homopolymer-check --low-complexity-check reads.ec.k$CK.fastq

# Merge simple, unbranched chains of vertices
$SGA_BIN fm-merge -m $MOL -t $CPU -o merged.k$CK.fa reads.ec.k$CK.filter.pass.fa

# Build an index of the merged sequences
$SGA_BIN index -d 1000000 -t $CPU merged.k$CK.fa

# Remove any substrings that were generated from the merge process
$SGA_BIN rmdup -t $CPU merged.k$CK.fa

# Compute the structure of the string graph
$SGA_BIN overlap -m $MOL -t $CPU merged.k$CK.rmdup.fa

# Perform the contig assembly without bubble popping
$SGA_BIN assemble -m $OL -g $MAX_GAP_DIFF -r $R -o assemble.m$OL merged.k$CK.rmdup.asqg.gz

#
# Scaffolding/Paired end resolution
# 
CTGS=assemble.m$OL-contigs.fa
GRAPH=assemble.m$OL-graph.asqg.gz

# Realign reads to the contigs
~/work/devel/sga/src/bin/sga-align --name celegans.pe $CTGS $IN1 $IN2

# Make contig-contig distance estimates
~/work/devel/sga/src/bin/sga-bam2de.pl -n $MIN_PAIRS --prefix libPE celegans.pe.bam

# Make contig copy number estimates
~/work/devel/sga/src/bin/sga-astat.py -m $MIN_LENGTH celegans.pe.refsort.bam > libPE.astat

$SGA_BIN scaffold -m $MIN_LENGTH --pe libPE.de -a libPE.astat -o scaffolds.n$MIN_PAIRS.scaf $CTGS
$SGA_BIN scaffold2fasta -m $MIN_LENGTH -a $GRAPH -o scaffolds.n$MIN_PAIRS.fa -d $SCAFFOLD_TOLERANCE --use-overlap --write-unplaced scaffolds.n$MIN_PAIRS.scaf
