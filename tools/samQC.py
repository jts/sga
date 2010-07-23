#! /nfs/team71/phd/js18/software/Python-2.6.4/python

import pysamhack as pysam
import sys
import getopt
import math
import numpy
import pylab
from string import maketrans

class ReportType:
    POS=1
    COUNT=2
    QUAL=3
    LENGTH=4

def reverseComplement(seq):
    intab = "ACGT"
    outtab = "TGCA"
    transtab = maketrans(intab, outtab)
    t = seq.translate(transtab)
    t = t[::-1]
    return t

def readFasta(filename):
    sequence = ""
    name = ""

    infile = open(filename)
    out = dict()
    while 1:
        line = infile.readline()
        if line.startswith(">"):
            if len(sequence) > 0:
                out[name] = sequence
                sequence = ""
            name = line[1:].split()[0]
        else:
            if line == "":
                break
            sequence += line.rstrip()

    if len(sequence) > 0:
        out[name] = sequence
    return out

def getReadCoverage(data, s1):
    for i, x in enumerate(s1):
        if i not in data:
            data[i] = 1
        else:
            data[i] += 1

def getErrorPositions(data, s1, s2):
    count = 0
    for i, x in enumerate(zip(s1, s2)):
        if x[0] != x[1]:
            if i not in data:
                data[i] = 1
            else:
                data[i] += 1
            count += 1
    return count

def getErrorCounts(data, s1, s2):
    count = 0;
    for i, x in enumerate(zip(s1, s2)):
        if x[0] != x[1]:
            count += 1

    if count not in data:
        data[count] = 1
    else:
        data[count] += 1
    return count

def getErrorQual(data, s1, s2, q):
    count = 0
    for x in zip(s1, s2, q):
        t = [0, 0]
        if x[2] in data:
            t = data[x[2]]
        t[0] += 1
        if x[0] != x[1]:
            t[1] += 1
            count += 1
        data[x[2]] = t
    return count

def illuminaQual2Score(c):
    s = ord(c) - 64
    q =  10 * math.log(1 + 10 ** s / 10.0) / math.log(10)
    return q

def sangerQual2Score(c):
    q = ord(c) - 33
    return q

def score2prob(q):
    p = 10 ** (-q / 10.0)
    return p

def hasIndel(a):
    for i in a.cigar:
        if i[0] != 0: # opcode for match
            return 1
    return 0

def usage():
    print 'usage: samQC.py [-p|-c|-s|-l] [-t val] [-n num] [-d filename] [-o filename] [--range chr:start:stop] in.bam|in.sam ref.fasta'
    print 'Options:'
    print '    -p       Output summary of number of base calling errors per read-position'
    print '    -c       Output summary of number of base calling errors per read'
    print '    -l       Output summary of read lengths'
    print '    -q       Output the number of mismatches at each quality symbol'
    print '    -s       Output general summary (total reads, total error rate, etc)'
    print '    -t=INT   Only use the first INT bases of each read'
    print '    -n=INT   Only use INT reads to calculate statistics'
    print '    -o=FILE  Write per-read statitics to FILE'
    print '    -d=FILE  Plot the result and save the figure in FILE'
    print '             If FILE is "-" the figure will be displayed immediately'
    print '    --range  The region of the reference to process stats for, in format chr:start:end'

# Set defaults and parse options
reportType = ReportType.POS
bOutputSummary = False

chr = None
start = None
stop = None
numMaxReads = 10000
trim = None
drawFile = None
total_errors = 0
outFilename = None

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'pcshqtl:d:n:o:', ["range="])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    
for (oflag, oarg) in opts:
        if oflag == '-p':
            reportType = ReportType.POS
        if oflag == '-c':
            reportType = ReportType.COUNT
        if oflag == '-q':
            reportType = ReportType.QUAL
        if oflag == '-l':
            reportType = ReportType.LENGTH
        if oflag == '-s':
            bOutputSummary = True
        if oflag == '-t':
            trim = int(oarg)
        if oflag == '-d':
            drawFile = oarg
        if oflag == '-o':
            outFilename = oarg
        if oflag == '-n':
            numMaxReads = int(oarg)
        if oflag == '-h':
            usage()
            sys.exit(0)
        if oflag == '--range':
            chr,start,stop = oarg.split(":")

if len(args) == 2:
    samFilename = args[0]
    refFilename = args[1]
else:
    usage()
    sys.exit(2)

# Read in the reference sequence
reference = readFasta(refFilename)

for i,s in reference.items():
    print i, len(s)

# Load the samfile
samfile = pysam.Samfile(samFilename, "rb")

# read stats output file
outFile = None
if outFilename is not None:
    outFile = open(outFilename, 'w')

# Intialize the count
data = dict()
readCoverage = dict()
readLength = dict()

nb = 0
nr = 0
num_indel = 0

iter = samfile.fetch(chr, start, stop)
for alignment in iter:
    if nr >= numMaxReads:
        break

    # Get the sequence in the reference that this read aligns to
    if alignment.is_unmapped or alignment.is_secondary:
        continue
    ref_name = samfile.getrname(alignment.rname)
    ref_slice = reference[ref_name][alignment.pos:alignment.pos+alignment.rlen]
    seq_slice = alignment.seq
    qual_slice = alignment.qual

    if trim != None:
        ref_slice = ref_slice[0:trim]
        seq_slice = seq_slice[0:trim]
        qual_slice = qual_slice[0:trim]

    if alignment.is_reverse:
        seq_slice = reverseComplement(seq_slice)
        ref_slice = reverseComplement(ref_slice)

    if hasIndel(alignment):
        num_indel += 1
    else:
        if alignment.seq.find("N") == -1:
            rl = len(seq_slice)
            nb += rl
            nr += 1

            # Parse the stats
            num_errors = 0
            if reportType == ReportType.POS:
                num_errors = getErrorPositions(data, seq_slice, ref_slice)
            elif reportType == ReportType.COUNT:
                num_errors = getErrorCounts(data, seq_slice, ref_slice)
            elif reportType == ReportType.QUAL:
                num_errors = getErrorQual(data, seq_slice, ref_slice, qual_slice)

            if outFile is not None:
                outFile.write(">" + alignment.qname + ' ' + str(num_errors) + '\n')
                outFile.write(alignment.seq + '\n')

            total_errors += num_errors

            # Accumulate read coverage
            getReadCoverage(readCoverage, seq_slice)

            # Accumulate read length
            if rl not in readLength:
                readLength[rl] = 1
            else:
                readLength[rl] += 1

# Set up the output and plotting
plot_x = []
plot_y = []

if reportType == ReportType.POS:
    so = 0
    eo = max(data.keys())
    sum = 0
    
    plot_x = [0] * (eo + 1)
    plot_y = [0] * (eo + 1)

    for i in range(so, eo+1):
        v1 = data[i] if i in data else 0
        v2 = readCoverage[i] if i in readCoverage else 1
        er = float(v1) / float(v2)
        sum +=  v1
        print i,v1,v2,er
        plot_x[i] = i
        plot_y[i] = er

if reportType == ReportType.COUNT:
    so = 0
    eo = max(data.keys())
    sum = 0
    
    plot_x = [0] * (eo + 1)
    plot_y = [0] * (eo + 1)

    for i in range(so, eo+1):
        v = data[i] if i in data else 0
        er = float(v) / float(nr)
        sum += v
        print i,v,er
        plot_x[i] = i
        plot_y[i] = er

if reportType == ReportType.LENGTH:
    so = 0
    eo = max(readLength.keys())
    sum = 0
    
    plot_x = [0] * (eo + 1)
    plot_y = [0] * (eo + 1)

    for i in range(so, eo+1):
        v = readLength[i] if i in readLength else 0
        print i,v
        plot_x[i] = i
        plot_y[i] = v       

if reportType == ReportType.QUAL:
    # Sort the dicttion
    s_keys = data.keys()
    s_keys.sort()
    plot_x = [0] * len(s_keys)
    plot_y = [0] * len(s_keys)

    for idx,i in enumerate(s_keys):
        t = data[i]
        if t[1] > 0:
            score = sangerQual2Score(i)
            prop = float(t[1]) / float(t[0])
            prob = score2prob(score)
            lprop = math.log(prop)
            lprob = math.log(prob)
            print i, score, prob, prop, lprob, lprop
            plot_x[idx] = lprop
            plot_y[idx] = lprob

if bOutputSummary:
    ter = float(total_errors) / float(nb)
    print "Total reads used:", nr
    print "Total bases:", nb
    print "Total errors:", total_errors
    print "Total error rate:", ter
    print "Num reads with indels:", num_indel

if drawFile != None:
    pylab.scatter(plot_x, plot_y)
    if drawFile != '-':
        pylab.savefig(drawFile)
    else:
        pylab.show()

