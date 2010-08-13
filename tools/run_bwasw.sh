#! /bin/sh

IN=$1
REF_PREFIX=$2 #/nfs/team71/phd/js18/work/devel/sga/data/ecoli_k12.fa
BWA=/software/solexa/bin/aligners/bwa/current/bwa
SAMTOOLS=/software/solexa/bin/aligners/samtools/current/samtools

# Run bwa	
$BWA bwasw $REF_PREFIX $IN > $IN.sam

