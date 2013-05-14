import sys, os.path
import matplotlib
matplotlib.use('Agg')
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
BRANCH_CLASSIFICATION_NAME =  "BranchClassification"
RANDOM_WALK_NAME = "RandomWalkLength"
FRAGMENT_SIZE_NAME = "FragmentSize"
QUALITY_SCORE_NAME = "QualityScores"
GC_DISTRIBUTION_NAME = "GCDistribution"
GENOME_SIZE_NAME = "GenomeSize"

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
    pl.ylim([0, 10000])
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
    pl.ylim([0, 40000])
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_kmer_distribution(pp, data):
    CUTOFF = 0.99
    names = data.keys()
    for name in names:
        k = data[name][KMER_DISTRIBUTION_NAME]['k']
        x = list()
        y = list()
        for t in data[name][KMER_DISTRIBUTION_NAME]['distribution']:
            x.append(t['kmer-depth'])
            y.append(t['count'])

        s = sum(y)
        # Normalize y and apply a cutoff
        nx = list()
        ny = list()
        cumulative_sum = 0

        for a,b in zip(x,y):
            fb = float(b) / s
            cumulative_sum += fb
            if a < 80 and cumulative_sum < CUTOFF:
                nx.append(a)
                ny.append(fb)
        pl.plot(nx, ny)
    
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
    pl.legend(names, loc="upper left")
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

def plot_branch_classification(pp, data):
    names = data.keys()
    out = list()

    data_type = ('variant', 'repeat', 'error')

    # Make a plot for each branch classification
    for dt in data_type:

        type_key = "num_%s_branches" % dt
        for name in names:

            kmers = list()
            branch_rate = list()
            for t in data[name][BRANCH_CLASSIFICATION_NAME]:

                if t['k'] >= 21:
                    kmers.append(t['k'])
                    branch_rate.append(float(t[type_key]) / t['num_kmers'])

            pl.plot(kmers, branch_rate, 'o-')

        pl.yscale('log')
        pl.xlabel("k")

        ylab = "Frequency of %s branches" % dt
        title = "Frequency of %s branches in the k-de Bruijn graph" % dt 
        pl.ylabel(ylab)
        pl.title(title)
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
    pl.ylabel("Duplicate Proportion")
    pl.xlabel("Sample")
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_genome_size(pp, data):
    names = data.keys()
    out = list()
    for name in names:
        s = data[name][GENOME_SIZE_NAME]['size']
        out.append(s / 1000000000.)
    width = 0.35
    ind = np.arange(len(names))
    pl.bar(ind, out, width)
    pl.xticks(ind + width/2, names)
    pl.title("Genome Size")
    pl.ylabel("Size (Gbp)")
    pl.xlabel("Sample")
    pl.savefig(pp, format='pdf')
    pl.close()

def plot_gc_distribution(pp, data):
    names = data.keys()
    gc_bin_size = 0.02

    for name in names:

        bins = data[name][GC_DISTRIBUTION_NAME]['gc_bins']
        read_data = data[name][GC_DISTRIBUTION_NAME]['read_gc_prop']
        base_line, = pl.plot(bins, read_data)

        if 'ref_gc_prop' in data[name][GC_DISTRIBUTION_NAME]:
            ref_data = data[name][GC_DISTRIBUTION_NAME]['ref_gc_prop']
            pl.plot(bins, ref_data, '--', color="red")
    pl.xlabel("GC content")
    pl.ylabel("Proportion")
    pl.title("GC histogram")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()

    # Plot the 2D histogram of coverage vs gc
    for name in names:
        x = [ i * 100 for i in data[name][GC_DISTRIBUTION_NAME]['gc_samples'] ]
        y = data[name][GC_DISTRIBUTION_NAME]['cov_samples']

        hist,xedges,yedges = np.histogram2d(x,y, bins=[20, 50], range=[ [0, 100.0], [0, 100] ])
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
        pl.imshow(hist.T,extent=extent,interpolation='nearest',origin='lower')
        pl.colorbar()
        pl.title(name + ' Read coverage vs GC content')
        pl.xlabel("GC %")
        pl.ylabel("Read coverage")
        pl.savefig(pp, format='pdf')
        pl.close()

def plot_fragment_sizes(pp, data):

    # Trim outliers from the histograms
    DENSITY_CUTOFF = 0.98
    names = data.keys()
    for name in names:
        h = data[name][FRAGMENT_SIZE_NAME]['sizes']
        sizes = {}
        for i in h:
            if i not in sizes:
                sizes[i] = 1
            else:
                sizes[i] += 1
        n = len(h)
        x = list()
        y = list()
        sum  = 0
        for i,j in sorted(sizes.items()):
            if sum < DENSITY_CUTOFF:
                f = float(j) / n
                x.append(i)
                y.append(f)
                sum += f
        pl.plot(x, y)

    pl.xlabel("Fragment Size (bp)")
    pl.ylabel("Proportion")
    pl.title("Estimated Fragment Size Histogram")
    pl.legend(names)
    pl.savefig(pp, format='pdf')
    pl.close()    

def plot_quality_scores(pp, data):
    names = data.keys()

    # Plot mean quality
    for name in names:

        mean_quality = data[name][QUALITY_SCORE_NAME]['mean_quality']
        indices = range(0, len(mean_quality))
        pl.plot(indices, mean_quality, linewidth=2)

    pl.ylim([0, 40])
    pl.xlabel("Base position")
    pl.ylabel("Mean Phred Score")
    pl.title("Mean quality score by position")
    pl.legend(names, loc="lower left")
    pl.savefig(pp, format='pdf')
    pl.close()

    # Plot >q30 fraction
    for name in names:

        q30_fraction = data[name][QUALITY_SCORE_NAME]['fraction_q30']
        indices = range(0, len(q30_fraction))
        pl.plot(indices, q30_fraction)

    pl.xlabel("Base position")
    pl.ylabel("Fraction at least Q30")
    pl.title("Fraction of bases at least Q30")
    pl.legend(names, loc="lower left")
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

# Configure the plot
matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

pp = PdfPages("test_report.pdf")


# Genome Characteristics
plot_genome_size(pp, data) if any_set_has_key(data, GENOME_SIZE_NAME) else 0
plot_branch_classification(pp, data) if any_set_has_key(data, BRANCH_CLASSIFICATION_NAME) else 0

# Quality/Error rate plots
plot_quality_scores(pp, data) if any_set_has_key(data, QUALITY_SCORE_NAME) else 0
plot_first_error_position(pp, data) if any_set_has_key(data, FIRST_ERROR_NAME) else 0
#plot_errors_per_base(pp, data) if any_set_has_key(data, ERRORS_PER_BASE_NAME) else 0
plot_pcr_duplicates(pp, data) if any_set_has_key(data, PCR_DUPLICATE_NAME) else 0
plot_fragment_sizes(pp, data) if any_set_has_key(data, FRAGMENT_SIZE_NAME) else 0

# Coverage plots
plot_kmer_distribution(pp, data) if any_set_has_key(data, KMER_DISTRIBUTION_NAME) else 0
plot_random_walk(pp, data) if any_set_has_key(data, RANDOM_WALK_NAME) else 0
plot_gc_distribution(pp, data) if any_set_has_key(data, GC_DISTRIBUTION_NAME) else 0

# Graph topology plots
plot_graph_complexity(pp, data) if any_set_has_key(data, GRAPH_COMPLEXITY_NAME) else 0
plot_mean_unipath_lengths(pp, data) if any_set_has_key(data, UNIPATH_LENGTH_NAME) else 0

pp.close()
