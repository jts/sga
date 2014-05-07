#! /usr/bin/env python
#
# sga-astat.py - Compute Myers' a-statistic for a set of contigs using read alignments
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
        self.nlen = len
        self.n = 0
        self.astat = 0.0
        self.bUnique = True

def computeAStat(arrivalRate, len, n):
    return arrivalRate * len - (n * math.log(2))

# Params
numContigsForInitialEstimate = 20
numIterations = 3
singleCopyThreshold = 30
bKeepDuplicates = True
minLength = 0
genomeSize = 0
arrivalRate = 0

def usage():
    print 'usage: sga-astat.py in.bam'
    print 'Compute Myers\' a-statistic for a set of contigs using the read alignments in in.bam'
    print 'Options:'
    print '    -m=INT          only compute a-stat for contigs at least INT bases in length'
    print '    -b=INT          use the longest INT contigs to perform the initial estimate'
    print '                    of the arrival rate (default: ' + str(numContigsForInitialEstimate) + ')' 
    print '    -n=INT          perform INT bootstrap iterations of the estimate'
    print '    -g=INT          use INT as the genome size instead of estimating it'
    print '    --no-duplicates do not use duplicate reads to calculate statistics'

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'm:b:n:g:', ['help', 'no-duplicates'])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

for (oflag, oarg) in opts:
        if oflag == '-m':
            minLength = int(oarg)
        if oflag == '-b':
            numContigsForInitialEstimate = int(oarg)
        if oflag == '-n':
            numIterations = int(oarg)
        if oflag == '-g':
            genomeSize = int(oarg)
        if oflag == '--no-duplicates':
            bKeepDuplicates = False
        if oflag == '--help':
            usage()
            sys.exit(1)

if len(args) == 0:
    print 'Error: a BAM file must be provided\n'
    usage()
    sys.exit(2)

bamFilename = args[0]

sys.stderr.write('Reading alignments from ' + bamFilename + '\n')
contigData = list()

# Read the contig names and lengths from the bam
bamFile = pysam.Samfile(bamFilename, "rb")

for (name, length) in zip(bamFile.references, bamFile.lengths):
    t = name, length, 0
    contigData.append(ContigData(name, length))
    #print 'Name: ' + name
    #print 'Length: ' + str(len)

# Read the alignments and populate the counts
totalReads = 0
sumReadLength = 0

last_ref_idx = -1
last_pos = -1

for alignment in bamFile:

    if alignment.is_unmapped:
        continue

    ref_idx = alignment.rname
    ref_name = bamFile.getrname(ref_idx)

    pos = alignment.pos

    if ref_idx == last_ref_idx and pos == last_pos and not bKeepDuplicates:
        continue

    # Update the count
    cd = contigData[ref_idx]
    cd.n += 1
    totalReads += 1
    sumReadLength += alignment.rlen
    assert(cd.name == ref_name)
    last_ref_idx = ref_idx
    last_pos = pos

avgReadLen = sumReadLength / totalReads
contigData.sort(key=lambda c : c.len, reverse=True)

# Compute the length of the contigs in number of positions
# that can generate a read of length avgReadLen. Using
# when calculating the expected number of reads is a better
# approximation for small contigs
for cd in contigData:
    cd.nlen = cd.len - avgReadLen + 1 if (cd.len > avgReadLen) else 0

# Estimate the initial arrival rate using the longest contigs
# if the genome size was not provided
bootstrapLen = 0;
bootstrapReads = 0;

if genomeSize == 0:
    for i in range(0, min(numContigsForInitialEstimate, len(contigData))):
        cd = contigData[i]
        bootstrapLen += cd.nlen
        bootstrapReads += cd.n

    arrivalRate = float(bootstrapReads) / float(bootstrapLen)
    genomeSize = int(totalReads / arrivalRate)

    sys.stderr.write('Initial arrival rate: ' + str(arrivalRate) + '\n');
    sys.stderr.write('Initial genome size estimate: ' + str(genomeSize) + '\n')

    for i in range(0, numIterations):

        bootstrapLen = 0;
        bootstrapReads = 0;
        for cd in contigData:
            cd.astat = computeAStat(arrivalRate, cd.nlen, cd.n)
            cd.bUnique = cd.astat > singleCopyThreshold

            if cd.len >= minLength and cd.bUnique:
                bootstrapLen += cd.nlen
                bootstrapReads += cd.n

        # Estimate arrival rate based on unique contigs
        arrivalRate = float(bootstrapReads) / float(bootstrapLen)
        genomeSize = int(totalReads / arrivalRate)
        sys.stderr.write('Iteration ' + str(i) + ' arrival rate: ' + str(arrivalRate) + '\n');
        sys.stderr.write('Iteration ' + str(i) + ' genome size estimate: ' + str(genomeSize) + '\n')

arrivalRate = float(totalReads) / genomeSize
sys.stderr.write('Using genome size: ' + str(genomeSize) + '\n')
sys.stderr.write('Using arrival rate: ' + str(arrivalRate) + '\n')

# Compute final astat values
for cd in contigData:
    cd.astat = computeAStat(arrivalRate, cd.nlen, cd.n)
    cd.bUnique = cd.astat > singleCopyThreshold

    if cd.len >= minLength and cd.bUnique:
        bootstrapLen += cd.nlen
        bootstrapReads += cd.n

sumUnique = 0
sumRepeat = 0
for cd in contigData:
    if cd.len >= minLength and cd.nlen > 0:
        print '%s\t%d\t%d\t%d\t%f\t%f' % (cd.name, cd.len, cd.nlen, cd.n, cd.n / (cd.nlen * arrivalRate), cd.astat)
        
        if cd.bUnique:
            sumUnique += cd.len
        else:
            sumRepeat += cd.len

sys.stderr.write('Sum unique bases in contigs >= %d bp in length: %d\n' % (minLength, sumUnique))
sys.stderr.write('Sum repeat bases in contigs >= %d bp in length: %d\n' % (minLength, sumRepeat))

