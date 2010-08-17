#! /usr/bin/python

import pysam
import sys
import getopt
import math

def usage():
    print 'usage: evalScaffolds.py SAM'
    print 'evaluate the correctness of contigs aligned to a reference using the alignments in SAM'
    print 'Options:'
    print '               --help      print help and exit'
    print ''

def op2char(op):
    
    if op == 0:
        return 'M'
    elif op == 1:
        return 'I'
    elif op == 2:
        return 'D'
    elif op == 3:
        return 'N'
    elif op == 4:
        return 'S'
    elif op == 5:
        return 'H'
    elif op == 6:
        return 'P'

def cigar2string(cigar):
    out = list()
    for t in cigar:
        op = op2char(t[0])
        l = t[1]
        s = str(l) + str(op)
        out.append(s)
    return ''.join(out)
        
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'h', ["help"])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    
for (oflag, oarg) in opts:
        if oflag == '-h' or oflag == '--help':
            usage()
            sys.exit(0)

if len(args) == 1:
    samFilename = args[0]
else:
    usage()
    sys.exit(2)

# Parse the contig alignments
samFile = pysam.Samfile(samFilename, "r")

alignments = dict()

for record in samFile:

    # Save the longest alignment seen for this contig
    if record.qname not in alignments:# or alignments[record.qname].alen < record.alen:
        alignments[record.qname] = record

for qname, record in alignments.items():
    if record.alen != record.rlen:
        cs = cigar2string(record.cigar)
        print '%s is possibly misassembled. alen: %d rlen: %d cigar: %s' % (record.qname, record.alen, record.rlen, cs) 
        
