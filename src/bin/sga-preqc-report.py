import sys, os.path
import pylab as pl
import numpy as np
import json
from collections import namedtuple
from operator import attrgetter
from matplotlib.backends.backend_pdf import PdfPages

# String constants
KMER_DISTRIBUTION_NAME = "KmerDistribution"
FIRST_ERROR_NAME = "FirstErrorPosition"
PCR_DUPLICATE_NAME = "PCRDuplicates"
ERRORS_PER_BASE_NAME = "ErrorsPerBase"
UNIPATH_LENGTH_NAME = "UnipathLength"
GRAPH_COMPLEXITY_NAME =  "LocalGraphComplexity"
RANDOM_WALK_NAME = "RandomWalkLength"
KMER_DISTRIBUTION_MAX = 80

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

# Returns true if any of the data sets has a field
# with the given key
def any_set_has_key(data, key):
    for n in data:
        if key in data[n]:
            return 1
    return 0

#
def plot_mean_unipath_lengths(pp, data):
    names = data.keys()
    
    for name in names:
        kmers = list()
        values = list()
        for t in data[name][UNIPATH_LENGTH_NAME]:
            kmers.append(t['k'])
            values.append(np.mean(t['walk_lengths']))
        pl.plot(kmers, values, 'o-')

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

        kmers = list()
        values = list()

        for t in data[name]['RandomWalkLength']:
            kmers.append(t['k'])
            values.append(np.mean(t['walk_lengths']))
        pl.plot(kmers, values, 'o-')

    pl.xlabel("k")
    pl.ylabel("Mean Random Walk Length")
    pl.title("Mean length of a random walk through the k-de Bruijn graph")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_kmer_distribution(pp, data):
    names = data.keys()
    for name in names:
        k = data[name][KMER_DISTRIBUTION_NAME]['k']
        x = list()
        y = list()
        for t in data[name][KMER_DISTRIBUTION_NAME]['distribution']:
            x.append(t['kmer-depth'])
            y.append(t['count'])

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
        
        indices = data[name][FIRST_ERROR_NAME]['indices']
        base_count = data[name][FIRST_ERROR_NAME]['base_count']
        error_count = data[name][FIRST_ERROR_NAME]['error_count']
        proportion_error = [ float(e) / b for b,e in zip(base_count, error_count) ]
        pl.plot(indices, proportion_error)

    pl.xlabel("k-mer Position")
    pl.ylabel("Proportion")
    pl.title("k-mer position of first error")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_errors_per_base(pp, data):
    names = data.keys()
    for name in names:

        base_count = data[name][ERRORS_PER_BASE_NAME]['base_count']
        error_count = data[name][ERRORS_PER_BASE_NAME]['error_count']
        positions = np.arange(len(base_count))
        proportion_error = [ float(e) / b for b,e in zip(base_count, error_count) ]
        pl.plot(positions, proportion_error)

    pl.xlabel("base position")
    pl.ylabel("Error rate")
    pl.title("Per-position error rate")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()    

def plot_graph_complexity(pp, data):
    names = data.keys()
    out = list()
    for name in names:

        kmers = list()
        branch_rate = list()

        for t in data[name][GRAPH_COMPLEXITY_NAME]:
            kmers.append(t['k'])
            branch_rate.append(float(t['num_branches']) / t['num_kmers'])

        pl.plot(kmers,branch_rate, 'o-')
    pl.yscale('log')
    pl.xlabel("k")
    pl.ylabel("High-coverage branch frequency")
    pl.title("Frequency of high-coverage branches in the graph")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_pcr_duplicates(pp, data):
    names = data.keys()
    out = list()
    for name in names:
        n_dups = data[name][PCR_DUPLICATE_NAME]['num_duplicates']
        n_pairs = data[name][PCR_DUPLICATE_NAME]['num_pairs']
        out.append(float(n_dups) / n_pairs)
    width = 0.35
    ind = np.arange(len(names))
    pl.bar(ind, out, width)
    pl.xticks(ind + width/2, names)
    pl.savefig(pp, format='pdf')
    pl.close()    

#
# Start of program
#

data = {}
for f in sys.argv[1:]:
    if os.path.getsize(f) > 0:
        name = os.path.splitext(os.path.basename(f))[0]
        data[name] = json.load(open(f, 'r'))

pp = PdfPages("test_report.pdf")

plot_kmer_distribution(pp, data) if any_set_has_key(data, KMER_DISTRIBUTION_NAME) else 0
plot_first_error_position(pp, data) if any_set_has_key(data, FIRST_ERROR_NAME) else 0
plot_pcr_duplicates(pp, data) if any_set_has_key(data, PCR_DUPLICATE_NAME) else 0
plot_errors_per_base(pp, data) if any_set_has_key(data, ERRORS_PER_BASE_NAME) else 0
plot_mean_unipath_lengths(pp, data) if any_set_has_key(data, UNIPATH_LENGTH_NAME) else 0
plot_random_walk(pp, data) if any_set_has_key(data, RANDOM_WALK_NAME) else 0
plot_graph_complexity(pp, data) if any_set_has_key(data, GRAPH_COMPLEXITY_NAME) else 0

pp.close()
