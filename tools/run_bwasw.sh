#! /bin/sh -x

IN=$1
REF_PREFIX=$2 #/nfs/team71/phd/js18/work/devel/sga/data/ecoli_k12.fa
BWA=/software/solexa/bin/aligners/bwa/current/bwa
SAMTOOLS=/software/solexa/bin/aligners/samtools/current/samtools

# Run bwa	
$BWA bwasw $REF_PREFIX $IN > $IN.sam

# Convert to SAM
#$BWA samse $REF_PREFIX $IN.sai $IN > $IN.sam

# Convert to BAM
#$SAMTOOLS view -Sb $IN.sam > $IN.tmp.bam

# Sort BAM
#$SAMTOOLS sort $IN.tmp.bam $IN.sorted

#mv $IN.sorted.bam $IN.bam

# Index BAM
#$SAMTOOLS index $IN.bam

# Cleanup
#rm $IN.tmp.bam
#rm $IN.sam
#rm $IN.sai
