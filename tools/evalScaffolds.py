#! /usr/bin/python

import pysam
import sys
import getopt
import math

class ContigPosition:
    def __init__(self, start, end):
        self.start = start
        self.end = end # not inclusive

    def __repr__(self):
        return str(self.start) + '-' + str(self.end)

    def len(self):
        return self.end - self.start

class ScaffoldNode:
    def __init__(self, id, dist, sd, comp, type):
        self.id = id
        self.comp = comp
        self.dist = dist 
        self.sd = sd
        self.type = type

def calcSum(data):
    sum = 0
    for l in data:
        sum += l
    return sum

def calcMax(data):
    data.sort()
    return data[-1]

def calcMean(data):
    return calcSum(data) / len(data)

def calcN50(data):
    sum = calcSum(data)
    data.sort()
    half_sum = sum / 2
    running = 0
    for i in data:
        running += i
        if running >= half_sum:
            return i
    return 0

def parseStartNode(text):
    return ScaffoldNode(text, 0, 0.0, 0, 'N')

def parseInteriorNode(text):
    fields = text.split(",")
    return ScaffoldNode(fields[0], int(fields[1]), float(fields[2]), int(fields[4]), fields[5])

# Calculate the end-to-start distance between the two alignments
def calcInnerDistance(a1, a2):
    if a1.start < a2.start:
        return a2.start - a1.end
    else:
        return a1.start - a2.end

# Calculate the start-to-end distance between the two alignments
def calcOuterDistance(a1, a2):
    if a1.start < a2.start:
        return a2.end - a1.start
    else:
        return a1.end - a2.start

# Find the closest alignment of id to prevAlign
def findClosestAlign(id, alignments, prevAlign):
    alignList = alignments[id]
    bestDist = 2 ** 32
    bestAlign = None
    for l in alignList:
        d = calcInnerDistance(l, prevAlign)
        if d < bestDist:
            bestDist = d
            bestAlign = l
    return bestAlign

def calculateDistanceVector(positions):
    distances = list()
    distances.append(0)

    prev = positions[0]
    for i in positions[1:]:
        d = calcInnerDistance(i, prev)
        distances.append(d)
        prev = i
    return distances

def evaluateScaffoldOrdering(idx, scaffold, alignments):
    print '\n*** Evaluating scaffold %d***' %(idx)
    # Get the alignment of the first component
    startNode = scaffold[0]

    # Get the primary alignment of the starting node
    prevAlign = alignments[startNode.id][0]

    positions = list()
    positions.append(prevAlign);
    for sn in scaffold[1:]:
        # Get the closest alignment of this node to the previous
        bestAlign = findClosestAlign(sn.id, alignments, prevAlign)
        positions.append(bestAlign)

    prev = positions[0]
    all_increasing = True
    all_decreasing = True
    for s in positions[1:]:
        if s is None:
            continue

        if s.start > prev.start:
            all_decreasing = False
        if s.start < prev.start:
            all_increasing = False
        prev = s
    
    ordered = all_increasing or all_decreasing
    correct_str = 'correct'
    if not ordered:
        correct_str = 'not correct'

    print 'Scaffold %d ordering is %s ' %(idx, correct_str)

    sum_bases = 0
    for p in positions:
        sum_bases += p.len()

    span = calcOuterDistance(positions[0], positions[-1])

    if ordered and len(scaffold) > 1:
        ref_distances = calculateDistanceVector(positions)
        
        est_dist = list()
        est_sd = list()
        differences = list()

        max_diff = 0
        sum_diff = 0
        n = 0
        
        for (s, d) in zip(scaffold, ref_distances):
            est_dist.append(s.dist)
            est_sd.append(s.sd)
            diff = abs(s.dist - d)
            differences.append(diff)

            n += 1
            sum_diff += diff
            if diff > max_diff:
                max_diff = diff

        print 'Positions: ' + str(positions)
        print 'Ref distances: ' + str(ref_distances)
        print 'Est distances: ' + str(est_dist)
        #print 'Est sd:        ' + str(est_sd)
        print 'Difference:    ' + str(differences)

        print '\nEstimate accuracy: Maximum difference: %d Mean difference: %f' %(max_diff, sum_diff / n)
        print 'Scaffold has %d components, contains %d bases and spans %d reference bases' %(n, sum_bases, span)
    return (ordered, sum_bases, span)

def usage():
    print 'usage: evalScaffolds.py SCAFFOLDS SAM'
    print 'evaluate the correctness of the SCAFFOLDS using the alignments of the components to the reference in SAM'
    print 'Options:'
    print '            --no-singleton Skip evaluating singleton scaffolds'
    print ''

bEvalSingletons = True
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'h', ["help", "no-singleton"])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    
for (oflag, oarg) in opts:
        if oflag == '-h' or oflag == '--help':
            usage()
            sys.exit(0)
        if oflag == '--no-singleton':
            bEvalSingletons = False

if len(args) == 2:
    scaffoldFilename = args[0]
    samFilename = args[1]
else:
    usage()
    sys.exit(2)

# Parse the contig alignments
samFile = pysam.Samfile(samFilename, "r")

alignments = dict()

iter = samFile.fetch()
for record in iter:
    print '%s pos: %d' %(record.qname, record.pos)
    if record.qname not in alignments:
        alignments[record.qname] = list()
    alignments[record.qname].append(ContigPosition(record.pos, record.pos + record.rlen))

# Evalulate the scaffolds, one by one
idx = 0
total = 0
num_correct = 0
num_wrong = 0

correct_bases = list()
correct_span = list()

for line in open(scaffoldFilename,'r'):
    line = line.rstrip('\n');
    fields = line.split("\t")
    
    scaffold = list()
    scaffold.append(parseStartNode(fields[0]))
    
    for f in fields[1:]:
        scaffold.append(parseInteriorNode(f))

    # Evaluate the scaffold
    if len(scaffold) > 1 or bEvalSingletons:
        [ordered, bases, ref_span] = evaluateScaffoldOrdering(idx, scaffold, alignments)

        if ordered:
            num_correct += 1
            correct_bases.append(bases)
            correct_span.append(ref_span)
        else:
            num_wrong += 1
        total += 1
    idx += 1

print '\n===\nTotal scaffolds: %d Correctly ordered: %d Incorrectly ordered: %d' %(total, num_correct, num_wrong)

print '\nStatistics for correct scaffolds:'
print 'Sum bases: %d' %(calcSum(correct_bases))
print 'Mean bases: %d' %(calcMean(correct_bases))
print 'N50 bases: %d' %(calcN50(correct_bases))
print 'Max bases: %d' %(calcMax(correct_bases))
print ''
print 'Sum reference span: %d' %(calcSum(correct_span))
print 'Mean reference span: %d' %(calcMean(correct_span))
print 'N50 span: %d' %(calcN50(correct_span))
print 'Max span: %d' %(calcMax(correct_span))
