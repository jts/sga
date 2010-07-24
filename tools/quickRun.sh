#! /bin/sh -x

IN=$1
PREFIX=$2

ML=50 #min length
COL=25 # overlap length to use for correction
AOL=40 # overlap length to use for assembly

$SGA preprocess -q 15 -m $ML $IN > $PREFIX.fastq
$SGA index $PREFIX.fastq
$SGA correct -m $COL -e 0.06 -l 16 -t 4 -r 4 -o $PREFIX.ec.fa $PREFIX.fastq
$SGA index $PREFIX.ec.fa
$SGA rmdup -e 0.04 -t 4 $PREFIX.ec.fa
$SGA index $PREFIX.ec.rmdup.fa
$SGA overlap -m $AOL -e 0 -t 4 -i $PREFIX.ec.rmdup.fa
$SGA assemble -r -t 10 -b 2 -o $PREFIX.contigs.fa $PREFIX.ec.rmdup.asqg.gz
