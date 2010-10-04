#! /usr/bin/python
# Convert a file of contigs and a bam file of read alignments to an AMOS afg file.
# Inspired by abyss2afg
import pysam
import sys
import getopt
from string import maketrans

#
# Classes
#
class Contig:
    def __init__(self, eid, seq):
        l = len(seq)
        qual = 'K' * l

        global currIID
        self.iid = currIID
        currIID += 1
        self.eid = eid
        self.seq = seq
        self.qual = qual
        self.tle = list()

    def addTLE(self, clr, off, readIID, readEID):
        data = (clr, off, readIID, readEID)
        self.tle.append(data)

    def printTiles(self):
        for tile in self.tle:
            print '{TLE'
            print 'clr:%d,%d' % (tile[0][0], tile[0][1])
            print 'off:%d' % (tile[1])
            print 'src:%d' % (tile[2])
            print 'eid:%s' % (tile[3])
            print '}'
    def printAFG(self):
        print '{CTG\niid:%d\neid:%s' % (self.iid, self.eid)
        print 'seq:'
        printSequence(self.seq)
        print '.'
        print 'qlt:'
        printSequence(self.qual)
        print '.'
        self.printTiles()
        print '}'

# 
class Fragment:
    def __init__(self, eid):
        global currIID
        self.iid = currIID
        currIID += 1
        self.eid = eid
        self.reads = list()

    def addRead(self, read):
        if self.getNumReads() == 2:
            print 'Error too many reads for fragment ' + self.eid
            sys.exit(0)
        self.reads.append(read)
    def getNumReads(self):
        return len(self.reads)
    
    def printAFG(self):
        self.reads.sort(key=lambda read: read.eid)
        read0 = self.reads[0]
        read1 = self.reads[1]
        read0.printAFG(self.iid)
        read1.printAFG(self.iid)

        print '{FRG'
        print 'rds:%d,%d' % (read0.iid, read1.iid)
        print 'lib:1'
        print 'eid:%s' % self.eid
        print 'iid:%d' % self.iid
        print 'typ:I\n}'

    def __repr__(self):
        return 'Fragment iid: ' + str(self.iid)

#
class Read:
    def __init__(self, eid, seq, qual):
        global currIID
        self.iid = currIID
        currIID += 1
        self.eid = eid
        self.seq = seq
        self.qual = qual

    def printAFG(self, fragIID):
        print '{RED'
        print 'clr:0,%d' % (len(self.seq))
        print 'eid:%s' % self.eid
        print 'iid:%d' % self.iid
    
        if fragIID is not None:
            print 'frg:%d' % (fragIID)

        print 'seq:'
        printSequence(self.seq)
        print '.\nqlt:'
        printSequence(self.qual)
        print '.\n}'
        
    def __repr__(self):
        return 'Read iid: ' + str(self.iid)

#
# Functions
#
def usage():
    print 'usage: sga2afg.py --pe-mean=N --pe-sd=SD <contigs fasta> <alignments bam>'
    print 'Create an AMOS message file using a fasta file of contigs and a bam file'
    print 'Options:'
    print '    -h       Print this message and exit'

def reverseComplement(seq):
    intab = "ACGT"
    outtab = "TGCA"
    transtab = maketrans(intab, outtab)
    t = seq.translate(transtab)
    t = t[::-1]
    return t

# Add a read to the fragment database
def addReadToFragment(fragmentName, read):
    global fragmentDB
    global currIID
    frag = None
    if fragmentName not in fragmentDB:
        frag = Fragment(fragmentName)
        fragmentDB[fragmentName] = frag
        currIID += 1
    else:
        frag = fragmentDB[fragmentName]
    frag.addRead(read)
    return frag

# If both reads of a fragment have been added, write the reads
# and the fragment as an AFG record and delete it from the DB
def processFragment(fragmentName):
    global fragmentDB
    frag = fragmentDB[fragmentName]
    if frag.getNumReads() == 2:
        frag.printAFG()
        del(fragmentDB[fragmentName])

# Print a sequence, break up long lines
def printSequence(seq):
    step = 80
    l = len(seq)
    for x in range(0, l, step):
        print seq[x:x+step]

# Read the alignments to the given contig
# from the bam file and build the AFG record
def processContig(contigEID, contigSeq, bamFile):
    global currIID

    contig = Contig(contigEID, contigSeq)
    contigLen = len(contigSeq)
    # Process all the reads aligned to this contig
    iter = bamFile.fetch(contigEID)
    n = 0

    for currRead in iter:
        
        if currRead.is_unmapped:
            continue

        fragmentName = currRead.qname
        readName = fragmentName
        if currRead.is_read1:
            readName += "/1"
        elif currRead.is_read2:
            readName += "/2"
        
        seqStr = currRead.seq
        if currRead.is_reverse:
            seqStr = reverseComplement(seqStr)

        qualStr = currRead.qual
        if qualStr is None:
            qualStr = 'K' * len(seqStr)

        readObj = Read(readName, seqStr, qualStr)

        # If the read is a pair and the mate is mapped, don't print it
        # until its mate has been added to the fragment DB
        if currRead.is_paired and not currRead.mate_is_unmapped:
            # Get the fragment IID
            addReadToFragment(fragmentName, readObj)
            processFragment(fragmentName)
        else:
            readObj.printAFG(None)

        # Add the tile to the contig
        clearCoord = [0, currRead.rlen]
        if currRead.is_reverse:
            start = clearCoord[0]
            end = clearCoord[1]
            clearCoord = [currRead.rlen - end, currRead.rlen - start]
            clearCoord.reverse();
        offset = currRead.pos
        contig.addTLE(clearCoord, offset, readObj.iid, readObj.eid)

    # Write the contig record
    contig.printAFG()

#
def printMeta():
    print '{UNV'
    print 'eid:afg'
    print 'com:'
    print 'generated by sga2afg'
    print '.'
    print '}\n'


#
def printLibrary(pe_mean, pe_sd):
    print '{LIB'
    print 'eid:1'
    print 'iid:1'
    print '{DST'
    print 'mea:', pe_mean
    print 'std:', pe_sd
    print '}'
    print '}\n'

#
# Main
#

pe_mean = 2000
pe_sd = 200

# Parse args
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'h', ["help", "pe-mean=", "pe-sd="])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
	
for (oflag, oarg) in opts:
        if oflag == '-h' or oflag == '--help':
            usage()
            sys.exit(0)
        elif oflag == '--pe-mean':
            pe_mean = int(oarg)
        elif oflag == '--pe-sd':
            pe_sd = int(oarg)
        else:
            print 'Unhandled argument: ' + str(oflag)
            sys.exit(2)

if len(args) == 2:
    contigsFilename = args[0]
    bamFilename = args[1]
else:
    usage()
    sys.exit(2)

# 
printMeta()
printLibrary(pe_mean, pe_sd)

bamFile = pysam.Samfile(bamFilename, "rb")

# Process every contig and write an AFG record
contigFile = open(contigsFilename)
currCID = ""
currSeq = ""
currIID = 2
fragmentDB = dict()

for line in contigFile:
    if line[0] == '>':
        if currCID != "":
            processContig(currCID, currSeq, bamFile)
        currCID = line[1:].split()[0]
        currSeq = ""
    else:
        currSeq += line.rstrip()

processContig(currCID, currSeq, bamFile)
