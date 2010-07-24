#! /nfs/team71/phd/js18/software/Python-2.6.4/python

import pysamhack as pysam
import sys
import getopt

def usage():
	print 'usage: extractRegion.py [--range chr:start:stop] in.bam'
	print 'Outputs all the reads in the specified region of in.bam as a FASTQ file'
	print 'Options:'
	print '    --range  The region of the reference to process stats for, in format chr:start:end'
	print '    -h       Print this message and exit'

# defaults
chr = None
start = None
stop = None

try:
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'h:', ["range="])
except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)
	
for (oflag, oarg) in opts:
		if oflag == '-h':
			usage()
			sys.exit(0)
		if oflag == '--range':
			chr,start,stop = oarg.split(":")

if len(args) == 1:
	samFilename = args[0]
else:
	usage()
	sys.exit(2)

print chr,start,stop
# Load the samfile
samfile = pysam.Samfile(samFilename, "rb")
iter = samfile.fetch(chr, int(start), int(stop))
n = 0
for currRead in iter:

	print '@%s' % n
	print currRead.seq
	print '+'
	print currRead.qual
	n += 1
