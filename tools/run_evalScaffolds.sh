#! /bin/sh

IN=$1
REF=$2
TOOLS_DIR=`dirname $0`
$TOOLS_DIR/breakScaffolds.pl $IN
$TOOLS_DIR/run_bwasw.sh scaffoldContigs.fa $REF
python $TOOLS_DIR/evalScaffolds.py scaffoldContigs.scaf scaffoldContigs.fa.sam
