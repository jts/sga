#! /bin/sh -x

IN=$1

if [ $# > 1 ]; then
    REF_PREFIX=$2
else
    REF_PREFIX=/nfs/team71/phd/js18/work/devel/sga/data/ecoli_k12.fa
fi

# These environment variables must be defined, probably better
# to require BWA/SAMTOOLS to be on the path
#BWA=/software/solexa/bin/aligners/bwa/current/bwa
#SAMTOOLS=/software/solexa/bin/aligners/samtools/current/samtools

# Run bwa	
$BWA aln $REF_PREFIX $IN > $IN.sai

# Convert to SAM
$BWA samse $REF_PREFIX $IN.sai $IN > $IN.sam

# Convert to BAM
$SAMTOOLS view -Sb $IN.sam > $IN.tmp.bam

# Sort BAM
$SAMTOOLS sort $IN.tmp.bam $IN.sorted

mv $IN.sorted.bam $IN.bam

# Index BAM
$SAMTOOLS index $IN.bam

# Cleanup
rm $IN.tmp.bam
rm $IN.sam
rm $IN.sai
