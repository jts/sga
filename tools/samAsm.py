#! /nfs/team71/phd/js18/software/Python-2.6.4/python

import pysamhack as pysam
import sys
import getopt
import math
import numpy
import pylab
from string import maketrans

class ReportType:
	OVR=1
	COV=2

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

def usage():
	print 'usage: samAsm.py [-n INT] [--range chr:start:stop] in.bam|in.sam ref.fasta'
	print 'Report read overlaps steps based on alignments to a reference. The input bam/sam must be sorted'
	print 'Options:'
	print '    --range  The region of the reference to process stats for, in format chr:start:end'
	print '    -n       Only use the first INT reads'
	print '    -h       Print this message and exit'


# Set defaults and parse options
reportType = ReportType.COV
bOutputSummary = False

chr = None
start = None
stop = None
numMaxReads = -1

try:
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'hn:', ["range="])
except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)
	
for (oflag, oarg) in opts:
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

# Set up variables
numUnmapped = 0
numMapped = 0
prevRead = None
coverage = dict()

# Load the samfile
samfile = pysam.Samfile(samFilename, "rb")
iter = samfile.fetch(chr, start, stop)

for currRead in iter:
	if numMaxReads > 0 and nr >= numMaxReads:
		break

	# Get the sequence in the reference that this read aligns to
	if currRead.is_unmapped:
		numUnmapped += 1
		continue
	else:
		numMapped += 1

	start = currRead.pos;
	stop = start + currRead.rlen
	for i in range(start, stop+1):
		if i not in coverage:
			coverage[i] = 1
		else:
			coverage[i] += 1
	
	if reportType == ReportType.OVR and prevRead != None:
		prev_start = prevRead.pos
		prev_stop = prev_start + prevRead.rlen
		overlap = min(prev_stop, stop) - max(prev_start, start)
		print currRead.qname, overlap
	
	prevRead = currRead

s_keys = coverage.keys()
s_keys.sort()
start = 0
stop = s_keys[-1]

for i in xrange(start, stop):
	if i in coverage:
		print i, coverage[i]
	else:
		print i, 0
#print 'Num mapped:', numMapped
#print 'Num unmapped:', numUnmapped
