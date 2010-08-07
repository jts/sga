#! /usr/bin/python
#
# a-stat.py - Compute Myers' a-statistic for a set of contigs using read alignments
# in a bam file
#
import pysam
import sys
import getopt
import math

class ContigData:
    def __init__(self, name, len):
        self.name = name
        self.len = len
        self.n = 0
        self.astat = 0.0
        self.bUnique = True

def computeAStat(arrivalRate, len, n):
    return arrivalRate * len - (n * math.log(2))

# Params
numContigsForInitialEstimate = 20
numIterations = 3
singleCopyThreshold = 17 # Myers' default

def usage():
    print 'usage: a-stat.py in.bam'
    print 'Compute Myers\' a-statistic for a set of contigs using the read alignments in in.bam'
    print 'Options:'
    print '    -m=INT   only compute a-stat for contigs at least INT bases in length'
    print '    -b=INT   use the longest INT contigs to perform the initial estimate'
    print '             of the arrival rate (default: ' + str(numContigsForInitialEstimate) + ')' 
    print '    -n=INT   perform INT bootstrap iterations of the estimate'

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'm:b:n:', ['help'])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

minLength = 0

for (oflag, oarg) in opts:
        if oflag == '-m':
            minLength = int(oarg)
        if oflag == '-b':
            numContigsForInitialEstimate = int(oarg)
        if oflag == '-n':
            numIterations = int(oarg)
        if oflag == '--help':
            usage()
            sys.exit(1)

if len(args) == 0:
    print 'A BAM file must be provided'
    usage()
    sys.exit(2)

bamFilename = args[0]

sys.stderr.write('Reading alignments from ' + bamFilename + '\n')
contigData = list()

# Read the contig names and lengths from the bam
bamFile = pysam.Samfile(bamFilename, "rb")

for (name, len) in zip(bamFile.references, bamFile.lengths):
    t = name, len, 0
    contigData.append(ContigData(name, len))
    #print 'Name: ' + name
    #print 'Length: ' + str(len)

# Read the alignments and populate the counts
totalReads = 0
iter = bamFile.fetch()
for alignment in iter:
    ref_idx = alignment.rname
    ref_name = bamFile.getrname(ref_idx)

    # Update the count
    cd = contigData[ref_idx]
    cd.n += 1
    totalReads += 1
    assert(cd.name == ref_name)

contigData.sort(key=lambda c : c.len, reverse=True)

# Estimate the initial arrival rate using the longest contigs
bootstrapLen = 0;
bootstrapReads = 0;

for i in range(0, numContigsForInitialEstimate):
    cd = contigData[i]
    bootstrapLen += cd.len
    bootstrapReads += cd.n

arrivalRate = float(bootstrapReads) / float(bootstrapLen)
genomeSize = int(totalReads / arrivalRate)

sys.stderr.write('Initial arrival rate: ' + str(arrivalRate) + '\n');
sys.stderr.write('Initial genome size estimate: ' + str(genomeSize) + '\n')

for i in range(0, numIterations):

    bootstrapLen = 0;
    bootstrapReads = 0;
    for cd in contigData:
        cd.astat = computeAStat(arrivalRate, cd.len, cd.n)
        cd.bUnique = cd.astat > singleCopyThreshold

        if cd.len >= minLength and cd.bUnique:
            bootstrapLen += cd.len
            bootstrapReads += cd.n

    # Estimate arrival rate based on unique contigs
    arrivalRate = float(bootstrapReads) / float(bootstrapLen)
    genomeSize = int(totalReads / arrivalRate)
    sys.stderr.write('Iteration ' + str(i) + ' arrival rate: ' + str(arrivalRate) + '\n');
    sys.stderr.write('Iteration ' + str(i) + ' genome size estimate: ' + str(genomeSize) + '\n')

for cd in contigData:
    if cd.len >= minLength:
        print '%s\t%d\t%d\t%f\t%f' % (cd.name, cd.len, cd.n, cd.n / (cd.len * arrivalRate), cd.astat)
