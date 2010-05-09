#ifndef UNIEST_H
#define UNIEST_H

#include "Util.h"
#include "UniData.h"
#include "Contig.h"
#include "FragmentDistribution.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>

//
// Typedefs
//
typedef std::set<std::string> StringSet;
typedef std::map<ContigID, UniData> IDUniDataMap;
typedef std::vector<const UniData*> UniDataPVec;
//
// Functions
//
UniDataPVec filterByPercent(const IDUniDataMap& udMap, double percent);
UniDataPVec filterByUnique(const IDUniDataMap& udMap);

void fitParameters(const UniDataPVec& fitVec, const IntDist& fragDist, bool fitError, double& mean, double& error_rate);
IntDist generatePELikelihood(IntDist& fragDist, int length);

// Various parsers
void parseContigs(std::string file, IDUniDataMap& covMap);
void parseAligns(std::string file, IDUniDataMap& covMap);
void parsePairedAligns(std::string file, IDUniDataMap& covMap);
void parseGraph(std::string file, IDUniDataMap& udMap);
void parseOptions(int argc, char** argv);

//
// Getopt
//
#define PROGRAM "UniEst"
static const char *VERSION_MESSAGE =
PROGRAM "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION] ... CONTIGS\n"
"Estimate the uniqueness of all the sequences in the CONTIGS file using alignment and paired end data.\n"
"\n"
"  -a, --align=FILE        file to read alignments from\n"
"  -p, --paired=FILE       file to read paired alignments from (experimental)\n"
"  -h, --histogram=FILE    file to read the fragment size histogram from\n"
"  -g, --graph=FILE        file to read the graph (adjacency info) from (experimental)\n"
"  -k, --kmer=KMER_SIZE    k-mer size\n"
"  -l, --len_cutoff=SIZE   suppress calls for contigs below this size\n"
"  -t, --threshold         log-likelihood difference required to make a copy number call\n"
"      --no_depth          do not use depth information when making calls\n"
"  -v, --verbose           display verbose output\n"
"  -o, --outfile           name of file to write contigs to\n"
"      --help              display this help and exit\n";


namespace opt
{
    static unsigned int k;
    static unsigned int verbose;
    static double threshold;
    static unsigned int length_cutoff;
    static bool bUseDepth;
    static bool bUsePairs;
    static bool bUseGraph;
    static std::string outfile;
    static std::string alignFile;
    static std::string pairedFile;
    static std::string histFile;
    static std::string graphFile;
}

static const char* shortopts = "k:o:t:a:p:h:l:g:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NODEPTH };

static const struct option longopts[] = {
    { "kmer",        required_argument, NULL, 'k' },
    { "align",       required_argument, NULL, 'a' },
    { "paired",      required_argument, NULL, 'p' },    
    { "histogram",   required_argument, NULL, 'h' },    
    { "threshold",   required_argument, NULL, 't' },
    { "outfile",     required_argument, NULL, 'o' },
    { "len_cutoff",  required_argument, NULL, 'l' },
    { "graph",       required_argument, NULL, 'g' },
    { "no_depth",    no_argument,       NULL, OPT_NODEPTH },
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


#endif
