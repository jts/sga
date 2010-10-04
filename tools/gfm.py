#! /usr/bin/python

# Prototype aligner for polymorphic reference sequences
import getopt
import sys
import re
import copy

# Classes

class Vertex:
    def __init__(self, idx, seq, preds):
        self.idx = idx;
        self.pos = 0
        self.seq = seq
        self.preds = preds
        self.b = 'X'
    def __repr__(self):
        out = str()
        out = "%d %s %s" % (self.idx, self.seq, self.preds)
        return out

# FMIndex class. Creates and stores the BWT, Occurrence and Predecessor arrays
class FMIndex:
    def __init__(self, filename):
        self.generateIndex(filename)
    
    def __repr__(self):
        out = str()
        #out += "Original table:\n"
        #for i in xrange(0, len(self.origTable)):
        #    out += "%d\t%s\n" % (i, self.origTable[i])
        
        out += "\nSuffix array table:\n"
        for i in xrange(0, len(self.saTable)):
            out += "%d\t%c\t%s\n" % (i, self.saTable[i].b, self.saTable[i])

        out += self.bwt
        return out

    # Index the reference
    def generateIndex(self, filename):
        # Convert the reference to a list of component strings
        # These are the vertices in the variation graph
        self.vertices = self.generateVertices(filename)
        print self.vertices
        # Compute the BWT by performing a simple block-sort of the vert strings
        self.blockSort(self.vertices)

    def blockSort(self,vertices):
        rotateList = list()

        # Generate all rotations of all the input strings
        # To support the branches, if the first character of the rotated string
        # is a split, we duplicate the entry with the multiplicty of the branch
        # This allows us to uniquely map from any entry describing a split into the predecessor
        # entry
        for v in vertices:
            rotations = self.generateRotations(v)
            for r in rotations:
                if ("ACGT$".find(r.seq[0]) == -1) and len(r.preds) > 1 and r.pos == 0:
                    print "Duplicating ", r
                    for i in xrange(0, len(r.preds)):
                        dup = Vertex(r.idx, r.seq, [r.preds[i]])
                        rotateList.append(dup)
                else:
                    rotateList.append(r)
        
        # Sort the rotation list into lexo order
        rotateList.sort(key=lambda c : c.seq)
        self.saTable = rotateList
        
        # Set up the mapping of split symbols in the BWT string to their corresponding entry
        # in the split section

        # Get the start index of the S section and initialize the BWT
        s_idx = -1
        for i in xrange(0, len(rotateList)):
            if rotateList[i].seq[0] == "S" and s_idx == -1:
                s_idx = i
            
        
        # Write the symbol of the predecessor nodes for each S in the BWT
        #for r in rotateList:
        #    if r.seq[0] == "S":
        #        predSeq = self.vertices[r.preds[0]].seq
        #        b = predSeq[-2]
        #        rotateList[s_idx].b = b
        #        s_idx += 1

        # Get the BWT string
        self.bwt = str()
        for r in rotateList:
            self.bwt += r.b
        print self

    def generateRotations(self, vertex):
        n = len(vertex.seq)
        rotated = copy.copy(vertex)
        rotations = list();
        for i in xrange(0, n):
            rotated = copy.copy(vertex)
            rotated.seq = rotated.seq[i:n] + "$"
            rotated.pos = i
            if i == 0:
                rotated.b = "$"
            else:
                rotated.b = vertex.seq[i - 1]
            rotations.append(copy.deepcopy(rotated))
        return rotations

    def generateVertices(self, filename):
        file = open(filename)
        seq = str()
        for line in file:
            if line[0] != '>':
                seq += line.rstrip()

        # Scan the string left to right. Variants are encoded regex-style within square brackets '[C|G]'
        # Each linear segment is indexed seperately and bounded by ambiguity codes, indicating the first
        # base of the branch.
        # Example: A[C|G]T would yield the components $AS, SCS, SGS, ST$.
        components = re.split('[\[\]]', seq)
        n = len(components)
        vertices = list()
        predecessors = list()
        prevCode = "$"
        nextCode = "$"
        idx = 0
        for i in xrange(0, n):
            currComp = components[i]
            isVariant = currComp.find("|") != -1
            isFirst = i == 0
            isLast = i == n - 1

            # Determine the left and right code for this segment
            # If the segment is a variant, the code is the ambiguity code
            # of the variant. Otherwise, the left code is the ambi code of the previous
            # variant, and the right code is the ambi code of the next variant
            if isVariant and (isFirst or isLast):
                print 'Error, reference cannot start or end with a variant in the prototype'
                sys.exit(1)
            
            leftCode = str()
            rightCode = str()

            if isVariant:
                ambiCode = getAmbiCode(currComp)
                leftCode = ambiCode
                rightCode = ambiCode
            else:
                leftCode = prevCode
                if not isLast:
                    rightCode = getAmbiCode(components[i+1])
                else:
                    rightCode = "$"

            # Fill in the internal sequence
            next = list()
            sub = currComp.split("|")
            for v in sub:
                out = Vertex(idx, leftCode + v + rightCode, predecessors)
                vertices.append(out)
                next.append(idx)
                idx += 1
            predecessors = next
            prevCode = rightCode
        return vertices

# Functions
def getAmbiCode(ambiStr):
    ambiDict = dict({'A':'A', 'C':'C', 'G':'G', 'T':'T', 'AC':'M', 'AG':'R', 'AT':'W', 'CG':'S', 'CT':'Y', 'GT':'K', 'ACG':'V', 'ACT':'H', 'AGT':'D', 'CGT':'B', 'ACGT':'N'})
    if ambiStr.find("|") == -1:
        print "Error: No ambiguity in segment", ambiStr
        sys.exit(1)
    print ambiStr
    sortedChars = ambiStr.split("|")
    sortedChars.sort()
    ambiLookup = "".join(sortedChars)
    return ambiDict[ambiLookup]

def usage():
    print 'usage: gfm.py target.fa query.fa'
    print 'Align the sequences in query.fa against the polymorphic reference in target.fa'
    print 'Polymorphic bases should be encoded in target.fa like A[C|G]A which can be ACA or AGA'
    print 'Options:'
    print '     --verbose     Print details of the algorithm'
#    print '    -m=INT   only compute a-stat for contigs at least INT bases in length'
#    print '    -b=INT   use the longest INT contigs to perform the initial estimate'
#    print '             of the arrival rate (default: ' + str(numContigsForInitialEstimate) + ')' 
#    print '    -n=INT   perform INT bootstrap iterations of the estimate'

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '', ['verbose', 'help'])
except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

#
bVerbose = False
for (oflag, oarg) in opts:
        if oflag == '--verbose':
            bVerbose = True
        if oflag == '--help':
            usage()
            sys.exit(1)

targetFilename = args[0]
queryFilename = args[1]

if bVerbose:
    print 'Target file: ' + targetFilename
    print 'Query file: ' + queryFilename

fmIndex = FMIndex(targetFilename)
