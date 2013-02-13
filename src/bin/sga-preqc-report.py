import sys, os.path
import pylab as pl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

TAG_UNIPATH_LENGTH    = 'UPL'
TAG_KMER_DISTRIBUTION = 'KMD'
KMER_DISTRIBUTION_MAX = 80

def test():
    X = np.linspace(-np.pi, np.pi, 256, endpoint=True);
    C, S = np.cos(X), np.sin(X)
    pl.plot(X, C)
    pl.plot(X, S)

    pp = PdfPages("test.pdf")
    pl.savefig(pp, format='pdf')
    pp.close()

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


#
# Start of program
#
data = {}
for f in sys.argv[1:]:
    name = os.path.splitext(os.path.basename(f))[0]
    load_data(data, name, f)

# Make the plots
pp = PdfPages("test_report.pdf")
if TAG_KMER_DISTRIBUTION in data:
    plot_kmer_distribution(pp, data[TAG_KMER_DISTRIBUTION])

if TAG_UNIPATH_LENGTH in data:
    plot_mean_unipath_lengths(pp, data[TAG_UNIPATH_LENGTH])
pp.close()
