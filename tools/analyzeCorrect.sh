#! /bin/sh -x

IN=$1

if [ $# > 1 ]; then
    REF_PREFIX=$2
else
    REF_PREFIX=/nfs/team71/phd/js18/work/devel/sga/data/ecoli_k12.fa
fi
TOOLS_DIR=`dirname $0`
echo $TOOLS_DIR
$TOOLS_DIR/run_bwa.sh $IN $REF_PREFIX
$TOOLS_DIR/samQC.py -s -p -n 100000 -o $IN.analyze $IN.bam $REF_PREFIX
