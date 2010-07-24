#! /bin/sh

SEQ=$1
REF_PREFIX=/nfs/team71/phd/js18/work/devel/sga/data/ecoli_k12.fa
TMPFILE=bwasw-input.fa
echo ">testseq" > $TMPFILE
echo $SEQ >> $TMPFILE

# Run bwa	
$BWA /software/solexa/bin/aligners/bwa/current/bwa bwasw $REF_PREFIX $TMPFILE

