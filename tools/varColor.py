#! /usr/bin/python
#
# varColor.py - Compute relative contributions of different individuals to a contig graph constructed from the population
#
import pysam
import sys
import getopt
import math

class ContigData:
    def __init__(self, name, len):
        self.name = name
        self.len = len
        self.nlen = len
        self.n = 0
        self.reads1 = 0
        self.reads2 = 0

# Params

def usage():
    print 'usage: varColor.py lib1.bam lib2.bam'
    print 'Compute relative contributions of different individuals to a contig graph constructed from the population'
#    print 'Options:'

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '', ['help'])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

minLength = 0

for (oflag, oarg) in opts:
        #if oflag == '-n':
        #    numIterations = int(oarg)
        if oflag == '--help':
            usage()
            sys.exit(1)

if len(args) < 2:
    print 'Two BAM files must be provided'
    usage()
    sys.exit(2)

filename1 = args[0]
filename2 = args[1]

# Read the contig names and lengths from the first bam
# It should be identical in both
bamFile1 = pysam.Samfile(filename1, "rb")
bamFile2 = pysam.Samfile(filename2, "rb")
contigData = list()

for (name, len) in zip(bamFile1.references, bamFile1.lengths):
    t = name, len, 0
    contigData.append(ContigData(name, len))
    #print 'Name: ' + name
    #print 'Length: ' + str(len)
print 'ID\tlength\treads1\treads2\n'

# Read the alignment data from each bam
# It may be faster to iterate over the alignments and lookup
# the contigs in the list but pysam throws a too many open files error on doing so
for cd in contigData:
    iter1 = bamFile1.fetch(cd.name)
    for alignment in iter1:
        cd.reads1 += 1

    iter2 = bamFile2.fetch(cd.name)
    for alignment in iter2:
        cd.reads2 += 1

    total = cd.reads1 + cd.reads2

    if total > 0:
        frac1 = float(cd.reads1) / total
        frac2 = float(cd.reads2) / total
    else:
        frac1 = 0.5
        frac2 = 0.5
    print '%s\t%d\t%d\t%d\t%d\t%f\t%f' % (cd.name, cd.len, cd.reads1, cd.reads2, total, frac1, frac2)
