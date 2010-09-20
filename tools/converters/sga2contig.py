#! /nfs/team71/phd/js18/software/Python-2.6.4/python

import pysam
import sys
import getopt

# Functions
def usage():
    print 'usage: sga2contig.py <contigs fasta> <alignments bam>'
    print 'Create an AMOS-style contig file using a fasta file of contigs and a bam file'
    print 'with reads aligned to the contigs.'
    print 'Options:'
    print '    -h       Print this message and exit'

# Print a sequence, break up long lines
def printSequence(seq):
    step = 80
    l = len(seq)
    for x in range(0, l, step):
        print seq[x:x+step-1]

#
def countAlignedReads(contigID, samFile):
    iter = samFile.fetch(contigID)
    n = 0
    for currRead in iter:
        n += 1
    return n

# 
def processRecord(contigID, contigSeq, samFile):

    # Count the number of reads aligned to this contig
    n = countAlignedReads(contigID, samFile)

    # Output the header line
    l = len(contigSeq)
    print '##{0} {1} {2} bases, {3} checksum'.format(contigID, n, l, '00000000')
    printSequence(contigSeq)

    # Process all the reads aligned to this contig
    iter = samFile.fetch(contigID)
    n = 0
    for currRead in iter:
        rcStr = "";
        clearCoord = [1, currRead.rlen]
        if currRead.is_reverse:
            rcStr = "RC"
            clearCoord.reverse()

        alignCoord = [currRead.pos, currRead.pos + currRead.rlen - 1]
        outName = currRead.qname
        if currRead.is_read1:
            outName += "/1"
        elif currRead.is_read2:
            outName += "/2"

        print('#{0}({1}) [{2}] {3} bases, {4} checksum. {{{5} {6}}} <{7} {8}>'.format(
               outName, currRead.pos, rcStr, currRead.rlen, '00000000', 
               clearCoord[0], clearCoord[1], alignCoord[0], alignCoord[1]))
        printSequence(currRead.seq)

#
# Main
#
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'h:')
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
	
for (oflag, oarg) in opts:
        if oflag == '-h':
            usage()
            sys.exit(0)

if len(args) == 2:
    contigsFilename = args[0]
    samFilename = args[1]
else:
    usage()
    sys.exit(2)

samFile = pysam.Samfile(samFilename, "rb")

# Process every contig and write out a record
contigFile = open(contigsFilename)
currID = ""
currSeq = ""

for line in contigFile:
    if line[0] == '>':
        if currID != "":
            processRecord(currID, currSeq, samFile)
        currID = line[1:].split()[0]
        currSeq = ""
    else:
        currSeq += line.rstrip()

processRecord(currID, currSeq, samFile)
