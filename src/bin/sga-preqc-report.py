#!/usr/bin/env python
"""Generate a readable report from preqc output.
"""

from __future__ import print_function, division

import sys, os.path
import matplotlib as MPL
import argparse
MPL.use('Agg') # this needs to be changed if using --show
import pylab as pl
import numpy as np
import json
import inspect
from collections import namedtuple
from operator import attrgetter
import matplotlib.transforms as transforms
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable

# String constants
KMER_DISTRIBUTION_NAME = "KmerDistribution"
FIRST_ERROR_NAME = "FirstErrorPosition"
PCR_DUPLICATE_NAME = "PCRDuplicates"
ERRORS_PER_BASE_NAME = "ErrorsPerBase"
UNIPATH_LENGTH_NAME = "UnipathLength"
DE_BRUIJN_SIMULATION_NAME = "SimulateAssembly"
GRAPH_COMPLEXITY_NAME =  "LocalGraphComplexity"
BRANCH_CLASSIFICATION_NAME =  "BranchClassification"
RANDOM_WALK_NAME = "RandomWalkLength"
FRAGMENT_SIZE_NAME = "FragmentSize"
QUALITY_SCORE_NAME = "QualityScores"
GC_DISTRIBUTION_NAME = "GCDistribution"
GENOME_SIZE_NAME = "GenomeSize"
RANDOM_WALK_LENGTH_NAME = "RandomWalkLength"
# other constants
PLOT_MARKERS = ['o','s','^','v','<','>','h','d','*'] # don't specify 12
# module scope variables
_PLOT_FUNCS = [] # filled in with functions for making individual plots


# top-level main
def main(argv=None):
    """Top-level main.  Called when files is used as a script.
    Primarially just handles the arguments using argpasre.
    :param argv: pass these to parse_args as instead of using sys.argv
        (default: None).
    """
    # options and argument definitions and defaults
    parser = argparse.ArgumentParser(description='Generate pdf report for sga-preqc.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    # input
    parser.add_argument('preqc_file', metavar='PREQC_FILE', nargs='+',
                    help='PreQC file(s) used to generate pdf report')
    # output
    parser.add_argument('-o', '--output', metavar='OUTPUT_PFX', default='preqc_report',
                        help='Report output prefix. (Default: preqc_report)')
    parser.add_argument('--png', help="Output a png file per figure",
                        action="store_true")
    # @TCC '--show' option causes difficulties with non-graphical terminals
    parser.add_argument('--show', help="Probably won't work; Show (interactive) figures using pylab",
                        action="store_true")
    # modes
    parser.add_argument('-p', '--page_per_plot', help="One plot per page/figure",
                        action="store_true")
    # selecting plots (overrides modes)
    parser.add_argument('-P', '--plot',
                        help="Select plot name/type\nimplies page-per-plot\ncan be given multiple times\nuse --list_plots to list valid values",
                        action="append")
    # misc
    parser.add_argument('-L', '--list_plots', help="List available plot names and exit",
                        action="store_true")
    # parse arguments
    args = parser.parse_args(argv)
    # misc: listing available plots
    if( args.list_plots ):
        plot_funcs, plotsample_funcs = _get_available_plot_functions()
        print("Available plots:")
        print('\n'.join(['\t'+k for k in sorted(plot_funcs)]))
        print("Available per-sample plots:")
        print('\n'.join(['\t'+k for k in sorted(plotsample_funcs)]))
        sys.exit(0)
    # main mode
    output_pfx = args.output
    preqc_files = args.preqc_file
    print('output_pfx: ', output_pfx)
    print('preqc_files:')
    print('\t' + '\n\t'.join(preqc_files))
    # load data
    data = load_preqc_datafiles(preqc_files)
    # make the report
    if( args.plot ):
        make_report_plot_per_page(args.plot, output_pfx, data,
                                  save_png=args.png, pylab_show=args.show)
    else:
        if( args.page_per_plot ):
            make_report_plot_per_page('all', output_pfx, data,
                                      save_png=args.png, pylab_show=args.show)
        else:
            make_report_with_subplots(output_pfx, data,
                                      save_png=args.png, pylab_show=args.show)


#######################################
## utility functions ##

def n50(values):
    """Return the N50 of the list of numbers"""
    values.sort(reverse=True)
    target = sum(values) / 2.
    total = 0
    for v in values:
        total += v
        if total >= target:
            return v
    return 0

def get_distinct_colors(nr):
    """get nr distinct fairly colorblind-safe colors
        :params nr: should be 1 to 12, above 12 returns 12
        :returns: list of colours in HTML hex
        adapted form work by Paul Tol Pieter van der Meer at
            SRON - Netherlands Institute for Space Research"""
    # color table in HTML hex format
    hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
               '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466',
               '#4477AA']
    # indicies in hexcols (combinations) for different numbers requested
    xarr = [[12],
            [12, 6],
            [12, 6, 5],
            [12, 6, 5, 3],
            [0, 1, 3, 5, 6],
            [0, 1, 3, 5, 6, 8],
            [0, 1, 2, 3, 5, 6, 8],
            [0, 1, 2, 3, 4, 5, 6, 8],
            [0, 1, 2, 3, 4, 5, 6, 7, 8],
            [0, 1, 2, 3, 4, 5, 9, 6, 7, 8],
            [0, 10, 1, 2, 3, 4, 5, 9, 6, 7, 8],
            [0, 10, 1, 2, 3, 4, 5, 9, 6, 11, 7, 8]]
    # check number requested
    if nr < 1 :
        nr = 1
    elif nr > 12:
        nr = 12
    return [ hexcols[xarr[nr-1][i]] for i in range(nr) ]


def path_splitall(path):
    """split a path into a list with entry for each part"""
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


def unique_names_from_filenames(filenames, splitext=False):
    """build unique names based on filenames adding path parts as needed"""
    if( splitext ):
        trim_func = lambda x: os.path.splitext(x)[0]
    else:
        trim_func = str
    # filenames need to be abspath to avoid ambiguity
    filenames = [os.path.abspath(f) for f in filenames]
    # start with the basenames
    names = [trim_func(os.path.basename(f)) for f in filenames]
    pidx = [1]*len(names) # how much of each filename to use for name
    # loop over all names resolving collisions
    for i in range(len(names)):
        while True:
            dups = [j for j,x in enumerate(names) if x == names[i] and i != j ]
            if( not dups ):
                break
            else:
                dups.append(i)
                for j in dups:
                    pidx[j] += 1
                    tmp = path_splitall(filenames[j])[-pidx[j]:]
                    names[j] = trim_func(os.path.join(*tmp))
    return names


## private helpers and code-fragment reuse functions ##

def _get_available_plot_functions_with_prefix(prefix):
    plot_funcs = {}
    for tmp in dir(sys.modules[__name__]):
        if( tmp.startswith(prefix) ):
            name = tmp[len(prefix):]
            tmp2 = getattr(sys.modules[__name__],tmp)
            if( hasattr(tmp2,'__call__') ): # callable, presumably a function
                assert( name not in plot_funcs )
                plot_funcs[name] = tmp2
    return plot_funcs

def _get_available_plot_functions():
    # make dicts of all the "plot_" and "plotsample_" functions associate with names
    # @TCC this will need to be changed if moved from module to class
    # basic plots (one plot for all samples)
    plot_funcs = _get_available_plot_functions_with_prefix('plot_')
    # per-sample plots (one plot for each samples)
    plotsample_funcs = _get_available_plot_functions_with_prefix('plotsample_')
    return plot_funcs, plotsample_funcs

def _create_fig_and_subplots(figname, subplot_rows, subplot_cols,
                             suptitle_base=None,
                             figs=None, subplots=None, extend_subplots=False):
    """creates a new figure with grid of subplots
        add fig to figs list and subplots to subplots list"""
    new_subplots = []
    # create figure
    fig = pl.figure()
    # add the figname to the figure object
    assert ( not hasattr(fig, 'figname') )
    fig.figname = figname
    # create the subplot axes
    for i in range(subplot_rows):
        for j in range(subplot_cols):
            new_subplots.append(fig.add_subplot(subplot_rows, subplot_cols, len(new_subplots)+1))
    # add the title (optionally)
    if( suptitle_base is not None ):
        fig.suptitle(suptitle_base+figname, size=16) # figure level title
    if( figs is not None ):
        figs.append(fig)
    if( subplots is not None ):
        if( extend_subplots ):
            subplots.extend(new_subplots)
        else:
            subplots.append(new_subplots)
    return fig, new_subplots

def _final_fig_layout(fig):
    fig.tight_layout()
    fig.subplots_adjust(top=0.93) # to make room for suptitle

def _finish_plot(ax, names, legend_loc=None, no_info_message="No Information"):
    """show a message in the axes if there is no data (names is empty)
       optionally add a legend
       return Fase if names is empty, True otherwise"""
    if( not names ):
        ax.text(0.5,0.5, no_info_message,
                fontweight='bold', va='center', ha='center',
                transform=ax.transAxes)
        return False
    if( legend_loc is not None ):
        ax.legend(names, loc=legend_loc)
    return True

def _fancy_barh(ax, values, data, val_fmt='', is_legend=False):
    """fancy-ish horizontal bar plot
       values must be same len as data, use np.nan if no values for sample
       :param is_legend: if True, will include legend like line segments
       :param val_fmt: is format string for value;
                       None don't show vlaue, '' just str()
    """
    assert( len(values) == len(data) )
    names = []
    for d in data:
        names.append(d['name'])
    width = 0.90
    bar_pos = np.arange(len(names))
    rects = ax.barh(bar_pos, values, width)
    ax.set_yticks([])
    ax.set_ylim([min(bar_pos)-(1-width)/2., max(bar_pos)+width+(1-width)/2.])
    ax.invert_yaxis()
    # set bar colors and annotate with sample name : size
    xmax = ax.get_xlim()[1]
    for ii,rect in enumerate(rects):
        rect.set_facecolor(data[ii]['plot_color'])
        x = 0.01
        y = rect.get_y()+rect.get_height()/2.
        if( val_fmt is None ):
            label = str(names[ii])
        else:
            if( np.isnan(values[ii]) ):
                label = '%s : No Info'%(names[ii])
            elif( val_fmt == '' ):
                label = '%s : %s'%(names[ii], str(values[ii]))
            else:
                label = ('%s : '+val_fmt)%(names[ii], values[ii])
        ax.text(x, y, label,
                va='center', ha='left',
                transform=transforms.blended_transform_factory(
                ax.transAxes, ax.transData),
                bbox=dict(boxstyle="round,pad=0.2", alpha=0.65, fc='w', lw=0) )
        if( is_legend ):
            # legend like line segment
            linex = [-0.12, -0.04]
            box_h= 0.5
            box_w_pad = 0.025
            ax.add_patch(MPL.patches.FancyBboxPatch((linex[0]-box_w_pad,y-box_h/2.),
                        linex[1]-linex[0]+2*box_w_pad, box_h,
                        ec='w',
                        fc='w',
                        boxstyle="square,pad=0",
                        transform=transforms.blended_transform_factory(
                                ax.transAxes, ax.transData),
                        clip_on=False) )
            line, = ax.plot([-0.12,-0.04], [y]*2,
                    '-', color=data[ii]['plot_color'],
                    marker=data[ii]['plot_marker'],
                    transform=transforms.blended_transform_factory(
                    ax.transAxes, ax.transData),
                    clip_on=False)
    ax.set_xlim((0,xmax))
    ax.set_ylabel(' \n \n ') # @TCC hack - fake ylabel so tight_layout adds spacing
    return True


## Plots of the data ##
# plot_<whatever> functions return False if no real data/info is plotted

def plot_legend(ax, data):
    """make and overall legend which takes up an axes by itself"""
    # create proxy artits
    proxy_arts = []
    names = []
    for d in data:
        #proxy_arts.append(MPL.patches.Rectangle((0,0),1,1, fc=d['plot_color']))
        proxy_arts.append(MPL.lines.Line2D((0,0),(1,1),
                linestyle='-', marker=d['plot_marker'],
                color=d['plot_color']))
        names.append(d['name'])
    ax.legend(proxy_arts, names, loc=2, bbox_to_anchor=(0,1), borderaxespad=0.)
    ax.set_frame_on(False)
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    return True


def plot_genome_size(ax, data, is_legend=False):
    values = []
    for d in data:
        s = np.nan
        if( GENOME_SIZE_NAME in d and
            'size' in d[GENOME_SIZE_NAME] ):
            s = d[GENOME_SIZE_NAME]['size'] / 1.0e6
        values.append(s)
    _fancy_barh(ax, values, data, val_fmt='%0.1f', is_legend=is_legend)
    ax.set_title("Est. Genome Size")
    ax.set_xlabel("Size (Mbp)")
    if( not is_legend ):
        ax.set_ylabel("Sample")
    return True


def plot_branch_classification_variant(ax, data, legend_loc=0):
    return _plot_branch_classification('variant', ax, data, legend_loc)

def plot_branch_classification_repeat(ax, data, legend_loc=0):
    return _plot_branch_classification('repeat', ax, data, legend_loc)

def plot_branch_classification_error(ax, data, legend_loc=0):
    return _plot_branch_classification('error', ax, data, legend_loc)

def _plot_branch_classification(branch_type, ax, data, legend_loc):
    type_key = "num_%s_branches" % branch_type
    names = []
    for d in data:
        if( BRANCH_CLASSIFICATION_NAME in d ):
            kmers = list()
            branch_rate = list()
            for t in d[BRANCH_CLASSIFICATION_NAME]:
                # We require at least 2 branches to be classified
                # in this type to plot this k-mer. This is to prevent
                # super small rates (like when classifying corrected data)
                # making the scale huge
                if t['k'] >= 21 and t[type_key] > 2:
                    kmers.append(t['k'])
                    branch_rate.append(float(t[type_key]) / t['num_kmers'])
            ax.plot(kmers, branch_rate, '-',
                    marker=d['plot_marker'], color=d['plot_color'])
            names.append(d['name'])
            have_something_to_plot = True
    # scale, lables, ect
    ax.set_yscale('log')
    ax.set_xlabel("k")
    ax.set_ylabel("Frequency of %s branches"%(branch_type))
    ax.set_title("%s branches in k-de Bruijn graph"%(branch_type))
    return _finish_plot(ax, names, legend_loc,
                        "No %s braches information"%(branch_type))


def plot_graph_complexity(ax, data, legend_loc=0):
    names = []
    for d in data:
        if( GRAPH_COMPLEXITY_NAME in d):
            kmers = list()
            branch_rate = list()
            for t in d[GRAPH_COMPLEXITY_NAME]:
                kmers.append(t['k'])
                branch_rate.append(float(t['num_branches']) / t['num_kmers'])
            ax.plot(kmers,branch_rate, '-',
                    marker=d['plot_marker'], color=d['plot_color'])
            names.append(d['name'])
    ax.set_yscale('log')
    ax.set_xlabel("k")
    ax.set_ylabel("High-coverage branch frequency")
    ax.set_title("High-coverage branches in graph")
    return _finish_plot(ax, names, legend_loc,
                        "No graph complexity information")


def plot_mean_unipath_lengths(ax, data, legend_loc=0):
    names = []
    for d in data:
        if UNIPATH_LENGTH_NAME in d:
            kmers = list()
            values = list()
            for t in d[UNIPATH_LENGTH_NAME]:
                if( t['walk_lengths'] ):
                    kmers.append(t['k'])
                    values.append(np.mean(t['walk_lengths']))
            ax.plot(kmers, values, '-',
                    marker=d['plot_marker'], color=d['plot_color'])
            names.append(d['name'])
    ax.set_xlabel("k")
    ax.set_ylabel("Mean unique path length")
    ax.set_title("Mean length of unambiguous segments\nin k-de Bruijn graph")
    return _finish_plot(ax, names, legend_loc,
                        "No unique graph segment information")


def plot_de_bruijn_simulation_lengths(ax, data, legend_loc=0):
    names = []
    for d in data:
        kmers = list()
        values = list()
        if DE_BRUIJN_SIMULATION_NAME in d:
            for t in d[DE_BRUIJN_SIMULATION_NAME]:
                kmers.append(t['k'])
                values.append(n50(t['walk_lengths']))
            ax.plot(kmers, values, '-',
                    marker=d['plot_marker'], color=d['plot_color'])
            names.append(d['name'])
    ax.set_xlabel("k")
    ax.set_ylabel("Simulated contig length N50")
    ax.set_title("Simulated contig lengths vs k")#\nin the k-de Bruijn graph")
    #pl.ylim([0, 50000])
    return _finish_plot(ax, names, legend_loc,
                        "No de Bruin simulation information")


def plot_mean_quality_scores(ax, data, legend_loc=0, use_markers=False):
    # Plot mean quality
    names = []
    for d in data:
        if( QUALITY_SCORE_NAME in d and
            'mean_quality' in d[QUALITY_SCORE_NAME] ):
            names.append(d['name'])
            mean_quality = d[QUALITY_SCORE_NAME]['mean_quality']
            indices = range(0, len(mean_quality))
            marker = d['plot_marker'] if( use_markers ) else None
            ax.plot(indices, mean_quality, '-',
                    marker=marker, color=d['plot_color'])
    #pl.ylim([0, 40])
    ax.set_xlabel("Base position")
    ax.set_ylabel("Mean Phred Score")
    ax.set_title("Mean quality score by position")
    return _finish_plot(ax, names, legend_loc,
                        "No mean quality information")


def plot_q30_quality_scores(ax, data, legend_loc=0, use_markers=False):
    # Plot >q30 fraction
    names = []
    for d in data:
        if( QUALITY_SCORE_NAME in d and
            'fraction_q30' in d[QUALITY_SCORE_NAME] ):
            names.append(d['name'])
            q30_fraction = d[QUALITY_SCORE_NAME]['fraction_q30']
            indices = range(0, len(q30_fraction))
            marker = d['plot_marker'] if( use_markers ) else None
            ax.plot(indices, q30_fraction, '-',
                    marker=marker, color=d['plot_color'])
    ax.set_xlabel("Base position")
    ax.set_ylabel("Fraction at least Q30")
    ax.set_title("Fraction of bases at least Q30")
    return _finish_plot(ax, names, legend_loc,
                        "No >q30 quality information")


def plot_errors_per_base(ax, data, legend_loc=0, use_markers=False):
    names = []
    for d in data:
        if( ERRORS_PER_BASE_NAME in d and
            'base_count' in d[ERRORS_PER_BASE_NAME] and
            'error_count' in d[ERRORS_PER_BASE_NAME] ):
            names.append(d['name'])
            base_count = d[ERRORS_PER_BASE_NAME]['base_count']
            error_count = d[ERRORS_PER_BASE_NAME]['error_count']
            positions = np.arange(len(base_count))
            proportion_error = [ float(e) / b for b,e in zip(base_count, error_count) ]
            marker = d['plot_marker'] if( use_markers ) else None
            ax.plot(positions, proportion_error, '-',
                    marker=marker, color=d['plot_color'])
    ax.set_xlabel("Base position")
    ax.set_ylabel("Error rate")
    ax.set_title("Per-position error rate")
    return _finish_plot(ax, names, legend_loc,
                        "No per-base error count information")


def plot_pcr_duplicates(ax, data, is_legend=False):
    values = []
    for d in data:
        s = np.nan
        if( PCR_DUPLICATE_NAME in d and
            'num_duplicates' in d[PCR_DUPLICATE_NAME] and
            'num_pairs' in d[PCR_DUPLICATE_NAME] ):
            n_dups = d[PCR_DUPLICATE_NAME]['num_duplicates']
            n_pairs = d[PCR_DUPLICATE_NAME]['num_pairs']
            s = float(n_dups) / float(n_pairs)
        values.append(s)
    _fancy_barh(ax, values, data, val_fmt=None, is_legend=is_legend)
    ax.set_title("Est. PCR Duplicate Proportion")
    ax.set_xlabel("Duplicate proportion")
    if( not is_legend ):
        ax.set_ylabel("Sample")
    return True


def plot_fragment_sizes(ax, data, legend_loc=0, use_markers=False):
    max_threshold_x = 0
    # Trim outliers from the histograms
    min_freq = 1.0e-5 # clips distribution tail below min_freq
    names = []
    for d in data:
        if( FRAGMENT_SIZE_NAME in d and
            'sizes' in d[FRAGMENT_SIZE_NAME] ):
            names.append(d['name'])
            h = d[FRAGMENT_SIZE_NAME]['sizes']
            sizes = {}
            for i in h:
                if i not in sizes:
                    sizes[i] = 1
                else:
                    sizes[i] += 1
            n = len(h)
            x = list()
            y = list()
            for i,j in sorted(sizes.items()):
                f = float(j) / float(n)
                x.append(i)
                y.append(f)
                if( f >= min_freq and i > max_threshold_x ):
                    max_threshold_x = i
            marker = d['plot_marker'] if( use_markers ) else None
            ax.plot(x, y, '-',
                    marker=marker, color=d['plot_color'])
    ax.set_xlim([0, max_threshold_x])
    ax.set_xlabel("Fragment Size (bp)")
    ax.set_ylabel("Proportion")
    ax.set_title("Estimated Fragment Size Histogram")
    return _finish_plot(ax, names, legend_loc,
                        "No fragment size information")


def plot_first_error_position(ax, data, legend_loc=0, use_markers=False):
    names = []
    for d in data:
        if( FIRST_ERROR_NAME in d and
            'indices' in d[FIRST_ERROR_NAME] and
            'base_count' in d[FIRST_ERROR_NAME] and
            'error_count' in d[FIRST_ERROR_NAME] ):
            names.append(d['name'])
            indices = d[FIRST_ERROR_NAME]['indices']
            base_count = d[FIRST_ERROR_NAME]['base_count']
            error_count = d[FIRST_ERROR_NAME]['error_count']
            proportion_error = [ float(e) / float(b) for b,e in zip(base_count, error_count) ]
            marker = d['plot_marker'] if( use_markers ) else None
            ax.plot(indices, proportion_error, '-',
                        marker=marker, color=d['plot_color'])
    ax.set_xlabel("k-mer position")
    ax.set_ylabel("Proportion")
    ax.set_title("k-mer position of first error")
    return _finish_plot(ax, names, legend_loc,
                        "No error positon information")


def plot_kmer_distribution(ax, data, legend_loc=0, use_markers=False):
    k = 'k'
    CUTOFF = 0.90
    MIN_DELTA = 0.001
    names = []
    for d in data:
        if( KMER_DISTRIBUTION_NAME in d and
            'k' in d[KMER_DISTRIBUTION_NAME] ):
            names.append(d['name'])
            k = d[KMER_DISTRIBUTION_NAME]['k']
            x = list()
            y = list()
            m = list()
            for t in d[KMER_DISTRIBUTION_NAME]['distribution']:
                x.append(t['kmer-depth'])
                y.append(t['count'])
            # Normalize y and apply a cutoff
            s = sum(y)
            nx = list()
            ny = list()
            cumulative_sum = 0
            for a,b in zip(x,y):
                fb = float(b) / s
                cumulative_sum += fb
                # These checks cut off the long tail of the distribution
                if cumulative_sum < CUTOFF or fb > MIN_DELTA:
                    nx.append(a)
                    ny.append(fb)
            marker = d['plot_marker'] if( use_markers ) else None
            ax.plot(nx, ny, '-',
                        marker=marker, color=d['plot_color'])
    ax.set_xlabel(str(k) + "-mer count")
    ax.set_ylabel("Proportion")
    ax.set_title(str(k) + "-mer count distribution")
    return _finish_plot(ax, names, legend_loc,
                        "No kmer distribution information")


def plot_random_walk(ax, data, legend_loc=0, use_markers=True):
    names = []
    for d in data:
        if( RANDOM_WALK_LENGTH_NAME in d ):
            names.append(d['name'])
            kmers = list()
            values = list()
            for t in d[RANDOM_WALK_LENGTH_NAME]:
                kmers.append(t['k'])
                values.append(np.mean(t['walk_lengths']))
            marker = d['plot_marker'] if( use_markers ) else None
            ax.plot(kmers, values, '-',
                        marker=marker, color=d['plot_color'])
    ax.set_xlabel("k")
    ax.set_ylabel("Mean Random Walk Length")
    ax.set_title("Mean length of a random walk\nthrough the k-de Bruijn graph")
    #ax.set_ylim([0, 40000])
    return _finish_plot(ax, names, legend_loc,
                        "No random walk information")


def plotsample_gc_distribution(ax, d):
    """Plot the gc distribution vs k-mer coverage
        for a SINGLE SAMPLE, unlike the plot_<whatever> functions
        :params d: is single sample preqc data, and element of the data list
    """
    # @TCC make check for missing data
    # Plot the 2D histogram of coverage vs gc
    if GC_DISTRIBUTION_NAME in d:
        x = [ i * 100 for i in d[GC_DISTRIBUTION_NAME]['gc_samples'] ]
        y = d[GC_DISTRIBUTION_NAME]['cov_samples']
        # Use the median to determine the range to show and round
        # to nearest 100 to avoid aliasing artefacts
        m = np.median(y)
        y_limit = np.ceil( 2*m / 100.) * 100.
        hist,xedges,yedges = np.histogram2d(x,y, bins=[20, 50], range=[ [0, 100.0], [0, y_limit] ])
        # draw the plot
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
        im = ax.imshow(hist.T,extent=extent,interpolation='nearest',origin='lower', aspect='auto')
        # colormap, colorbar, labels, ect.
        im.set_cmap('gist_heat_r')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        ax.figure.colorbar(im, cax=cax)
        ax.set_title(d['name'] + ' GC Bias')
        ax.set_xlabel("GC %")
        ax.set_ylabel("k-mer coverage")
        return True
    else:
        # no data, draw an empty plot
        return _finish_plot(ax, [], 0, 'No GC data for %s' % d['name'])


## creating the figures and saving them ##

def load_preqc_datafiles(preqc_files):
    """load the preqc_files, assign names
        also assign colors and markers"""
    # load the data
    data = []
    deserial_fail = []
    for f in preqc_files:
        f = os.path.abspath(f)
        if( os.path.getsize(f) <= 0 ):
            print("Warning: empty file '%s' ... skipping"%f)
            continue
        if( f in [d['file'] for d in data] ):
            print("Warning: duplicate file '%s' ... skipping"%f)
            continue
        try:
            deserial = json.load(open(f, 'r'))
        except ValueError:
            deserial_fail.append(f)
            continue
        data.append(deserial)
        data[-1]['file'] = f
    # create unique names for each entry
    names = unique_names_from_filenames([d['file'] for d in data], splitext=True)
    # generate the right number of distinct colors
    plot_colors = get_distinct_colors(len(names)) # returns up to 12
    for i in range(len(names)):
        data[i]['name'] = names[i]
        data[i]['plot_color'] = plot_colors[i%len(plot_colors)]
        data[i]['plot_marker'] = PLOT_MARKERS[i%len(PLOT_MARKERS)]
    for failed in deserial_fail:
        print("Warning: failed to de-serialize file:", failed)
    return data


def make_report_with_subplots(output_pfx, data, save_png=False, pylab_show=False):
    # Configure the plots (rc level at least)
    MPL.rc('figure', figsize=(8,10.5)) # in inches
    MPL.rc('font', size=10)
    MPL.rc('xtick', labelsize=8)
    MPL.rc('ytick', labelsize=8)
    MPL.rc('legend', fontsize=10)
    MPL.rc('axes', titlesize=12)
    MPL.rc('axes', labelsize=10)
    MPL.rcParams['lines.linewidth'] = 1.5

    # open PdfPages magic
    pdfpages = PdfPages(output_pfx+'.pdf')

    # create empty figures and subplots
    figs = []
    subplots = []

    _create_fig_and_subplots('fig1', 3, 2, 'SGA Preqc Results : ', figs, subplots)
    _create_fig_and_subplots('fig2', 3, 2, 'SGA Preqc Results : ', figs, subplots)
    _create_fig_and_subplots('fig3', 3, 2, 'SGA Preqc Results : ', figs, subplots)

    ## the plots on each axes/subplot ##

    ## fig1
    #plot_legend(subplots[0][0], data) # genome_size has legend info too
    # Genome Characteristics
    plot_genome_size(subplots[0][0], data, is_legend=True) # doubles as legend
    # Graph topology plots
    #plot_graph_complexity(subplots[0][2], data, legend_loc=None) # @TCC untested due to lack of data
    #plot_mean_unipath_lengths(subplots[0][2], data, legend_loc=None)
    plot_de_bruijn_simulation_lengths(subplots[0][2], data, legend_loc=None)
    # de-Bruijn graph branches info
    plot_branch_classification_variant(subplots[0][1], data, legend_loc=None)
    plot_branch_classification_repeat(subplots[0][3], data, legend_loc=None)
    plot_branch_classification_error(subplots[0][5], data, legend_loc=None)

    ## fig2
    # Quality/Error rate plots
    plot_pcr_duplicates(subplots[1][0], data, is_legend=True) # doubles as legend
    plot_fragment_sizes(subplots[1][2], data, use_markers=False, legend_loc=None)
    plot_first_error_position(subplots[1][4], data, use_markers=False, legend_loc=None)
    plot_mean_quality_scores(subplots[1][1], data, use_markers=False, legend_loc=None)
    plot_q30_quality_scores(subplots[1][3], data, use_markers=False, legend_loc=None)
    plot_errors_per_base(subplots[1][5], data, use_markers=False, legend_loc=None)

    ## fig3
    plot_legend(subplots[2][0], data)
    # Coverage plots
    plot_kmer_distribution(subplots[2][2], data, use_markers=False, legend_loc=None)
    #plot_random_walk(subplots[2][4], data, use_markers=False, legend_loc=None) #@jts deprecated

    # use the rest of the subplots for gc
    gc_subplots = []
    gc_subplots.append(subplots[2][1])
    gc_subplots.append(subplots[2][3])
    gc_subplots.append(subplots[2][5])

    ## gc... make figs as needed
    for d in data:
        if( len(gc_subplots) <= 0 ):
            _create_fig_and_subplots('fig'+str(len(figs)+1), 3, 2, 'SGA Preqc Results : ',
                                     figs, gc_subplots, extend_subplots=True)
        ax = gc_subplots.pop(0)
        plotsample_gc_distribution(ax, d)

    ## finalize and save the figures
    for fig in figs :
        _final_fig_layout(fig)
        if( save_png ):
            fig.savefig(output_pfx+'_'+fig.figname+'.png')
        if( pdfpages is not None ):
            fig.savefig(pdfpages, format='pdf')

    # finalize the pdfpages if needed
    if( pdfpages is not None ):
        pdfpages.close()
    # (optional) interactively 'show' the pylab figs
    if( pylab_show ):
        pl.show()
    pl.close() # cleanup


def make_report_plot_per_page(plots, output_pfx, data, save_png=False, pylab_show=False):
    """outputs a plot per page
        :param plots: can either be 'all' or a list of plot names
    """
    # Configure the plots (rc level at least)
    #MPL.rc('figure', figsize=(8,10.5)) # in inches
    MPL.rc('font', size=12)
    MPL.rc('xtick', labelsize=10)
    MPL.rc('ytick', labelsize=10)
    MPL.rc('legend', fontsize=12)
    MPL.rc('axes', titlesize=12)
    MPL.rc('axes', labelsize=12)
    MPL.rcParams['lines.linewidth'] = 1.5

    # open PdfPages magic
    pdfpages = PdfPages(output_pfx+'.pdf')

    plot_funcs, plotsample_funcs = _get_available_plot_functions()

    if( plots == 'all' ):
        plots = sorted(plot_funcs.keys())
        plots.extend(sorted(plotsample_funcs.keys()))
        # remove some plot_funcs
        plots.remove('legend')

    for plotname in plots:
        # handle either all-samples in one plot or per-sample plots
        if( plotname in plot_funcs ): # not per-sample plot?
            func = plot_funcs[plotname]
            data_list = [data]
        elif( plotname in plotsample_funcs ): # per-sample plot?
            func = plotsample_funcs[plotname]
            data_list = data
        else: # unrecogized, skip
            print("Warning: plot name '"+plotname+"' not recognized ... skipping",
                  file=sys.stderr)
            continue

        # loop over each sample if needed, else d == data (data_list=[data])
        for d in data_list:
            figname = plotname
            # append sample name to figname if there are multiple samples in data_list
            if( len(data_list) > 1 ):
                figname += '_'+d['name'].replace('/','_').replace('\\','_')
            # make the fig
            print('plotting:',figname)
            fig, subplots = _create_fig_and_subplots(figname, 1, 1)
            # plot
            rv = func(subplots[0], d)
            # finalize and save the fig
            if( rv ): # only save if data plotted
                _final_fig_layout(fig)
                if( save_png ):
                    fig.savefig(output_pfx+'_'+fig.figname+'.png')
                if( pdfpages is not None ):
                    fig.savefig(pdfpages, format='pdf')
            if( not pylab_show ): # we're done with fig if not showing
                pl.close()

    # finalize the pdfpages if needed
    if( pdfpages is not None ):
        pdfpages.close()
    # (optional) interactively 'show' the pylab figs
    if( pylab_show ):
        pl.show()
    pl.close() # cleanup


###############
## Main loop hook... if run as script run main, else this is a module ##
if __name__ == "__main__":
    sys.exit(main(argv=None))
