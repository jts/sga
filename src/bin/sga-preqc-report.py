import sys, os.path
import pylab as pl
import numpy as np
from collections import namedtuple
from operator import attrgetter
from matplotlib.backends.backend_pdf import PdfPages

TAG_UNIPATH_LENGTH    = 'UPL'
TAG_KMER_DISTRIBUTION = 'KMD'
TAG_GRAPH_COMPLEXITY = 'LGC'
TAG_RANDOM_WALK = 'RWL'
TAG_FIRST_ERROR_POSITION = 'PFE'
KMER_DISTRIBUTION_MAX = 80

def test():
    X = np.linspace(-np.pi, np.pi, 256, endpoint=True);
    C, S = np.cos(X), np.sin(X)
    pl.plot(X, C)
    pl.plot(X, S)

    pp = PdfPages("test.pdf")
    pl.savefig(pp, format='pdf')
    pp.close()

# Return the N50 of the list of numbers
def n50(values):
    values.sort(reverse=True)
    target = sum(values) / 2
    total = 0
    for v in values:
        total += v
        if total >= target:
            return v
    return 0

#
def plot_mean_unipath_lengths(pp, data):
    names = data.keys()
    
    for name in names:
        kmers = data[name].keys()
        kmers.sort()
        means = list()
        for k in kmers:
            means.append(np.mean(data[name][k]))

        pl.plot(kmers, means, 'o-')

    pl.xlabel("k")
    pl.ylabel("Mean unipath length")
    pl.title("Mean length of unambiguous segments of the k-de Bruijn graph")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

#
def plot_random_walk(pp, data):
    names = data.keys()
    
    for name in names:
        kmers = data[name].keys()
        kmers.sort()
        means = list()
        for k in kmers:
            values = map(attrgetter('walk_length'), data[name][k])
            means.append(np.mean(values))
        pl.plot(kmers, means, 'o-')

    pl.xlabel("k")
    pl.ylabel("Mean Random Walk Length")
    pl.title("Mean length of a random walk through the k-de Bruijn graph")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_kmer_distribution(pp, data):
    names = data.keys()
    for name in names:
        k = data[name].keys()[0]
        
        x = range(1, KMER_DISTRIBUTION_MAX);

        y = [ data[name][k][v] if v in data[name][k] else 0 for v in x ]
        s = sum(y)

        # Normalize y
        ny = [ float(v) / s for v in y ]

        pl.plot(x, ny)

    pl.xlabel(str(k) + "-mer count")
    pl.ylabel("Proportion")
    pl.title(str(k) + "-mer count distribution")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_first_error_position(pp, data):
    names = data.keys()
    for name in names:
        positions = sorted(data[name].keys())
        proportion_error = [ float(data[name][p][0].errors) / data[name][p][0].samples for p in positions ]
        pl.plot(positions, proportion_error)

    pl.xlabel("k-mer Position")
    pl.ylabel("Proportion")
    pl.title("Position of first error in the read")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_graph_complexity(pp, data):
    names = data.keys()
    out = list()
    for name in names:

        k = data[name].keys()
        x = sorted(data[name].keys())
        y = list()

        for k in x:
            y.append(float(data[name][k][0].num_branches) / data[name][k][0].num_kmers)
        pl.plot(x,y, 'o-')
    pl.yscale('log')
    pl.xlabel("k")
    pl.ylabel("High-coverage branch frequency")
    pl.title("Frequency of high-coverage branches in the graph")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

#
def parse_unipath_length(data, fields):

    # Parse a unipath length record
    k,l = int(fields[1]), int(fields[2])

    # Initialize dictionaries
    if k not in data:
        data[k] = list()

    data[k].append(l)

#
def parse_kmer_distribution(data, fields):

    # Parse a unipath length record
    k,occ,count = int(fields[1]), int(fields[2]), int(fields[3])

    if k not in data:
        data[k] = dict()
    data[k][occ] = count

def parse_graph_complexity(data, fields):
    output = map(int, fields[1:])
    k = output[0]
    if k not in data:
        data[k] = list()
    d = GraphComplexityTuple(output[0], output[1], output[2])
    data[k].append(d)

def parse_random_walk(data, fields):
    output = map(int, fields[1:])
    k = output[0]
    if k not in data:
        data[k] = list()
    d = RandomWalkTuple(output[0], output[1])
    data[k].append(d)

def parse_first_error_position(data, fields):
    output = map(int, fields[1:])
    position = output[0]
    if position not in data:
        data[position] = list()
    d = FirstErrorPositionTuple(output[0], output[1], output[2])
    data[position].append(d)

def load_data(data, name, filename):
    file = open(filename, 'r')

    for line in file:
        line = line.rstrip().split('\t')
        
        # Parse the type of data
        tag = line[0]

        if len(line[0]) != 3:
            continue

        if tag not in data:
            data[tag] = dict()
        if name not in data[tag]:
            data[tag][name] = dict()

        if tag == TAG_UNIPATH_LENGTH:
            parse_unipath_length(data[tag][name], line)
        elif tag == TAG_KMER_DISTRIBUTION:
            parse_kmer_distribution(data[tag][name], line)
        elif tag == TAG_GRAPH_COMPLEXITY:
            parse_graph_complexity(data[tag][name], line)
        elif tag == TAG_RANDOM_WALK:
            parse_random_walk(data[tag][name], line)
        elif tag == TAG_FIRST_ERROR_POSITION:
            parse_first_error_position(data[tag][name], line)

#
# Start of program
#

# Describe tuples for type of data
GraphComplexityTuple = namedtuple('GraphComplexityTuple', 'k num_kmers num_branches')
RandomWalkTuple = namedtuple('RandomWalkTuple', 'k walk_length')
FirstErrorPositionTuple = namedtuple('FirstErrorPositionTuple', 'position samples errors')

# Load the data files
data = {}
for f in sys.argv[1:]:
    name = os.path.splitext(os.path.basename(f))[0]
    load_data(data, name, f)

# Make the plots
pp = PdfPages("test_report.pdf")
if TAG_KMER_DISTRIBUTION in data:
    plot_kmer_distribution(pp, data[TAG_KMER_DISTRIBUTION])

if TAG_FIRST_ERROR_POSITION in data:
    plot_first_error_position(pp, data[TAG_FIRST_ERROR_POSITION])

if TAG_RANDOM_WALK in data:
    plot_random_walk(pp, data[TAG_RANDOM_WALK])

if TAG_GRAPH_COMPLEXITY in data:
    plot_graph_complexity(pp, data[TAG_GRAPH_COMPLEXITY])

if TAG_UNIPATH_LENGTH in data:
    plot_mean_unipath_lengths(pp, data[TAG_UNIPATH_LENGTH])
pp.close()
