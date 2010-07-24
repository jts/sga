#! /bin/sh -x

IN=$1

if [ $# > 1 ]; then
    REF_PREFIX=$2
else
    REF_PREFIX=/nfs/team71/phd/js18/work/devel/sga/data/ecoli_k12.fa
fi

~/work/devel/sga/tools/run_bwa.sh $IN $REF_PREFIX
~/work/devel/sga/tools/samQC.py -s -p -n 100000 -o $IN.analyze $IN.bam $REF_PREFIX
