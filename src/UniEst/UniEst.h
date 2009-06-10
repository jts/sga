#ifndef UNIEST_H
#define UNIEST_H

#include "Util.h"
#include "UniData.h"
#include "Contig.h"
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
typedef std::map<ContigID, UniData> CoverageMap;


//
// Functions
//
double estimateGlobalMeanCoverage(CoverageMap& cov_map);
void estimateUniquenessByDepth(CoverageMap& cov_map, double mean_param);

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
"Usage: " PROGRAM " [OPTION] ... CONTIGS ALIGNMENTS\n"
"Estimate the uniqueness of all the sequences in the CONTIGS file.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -v, --verbose         display verbose output\n"
"  -t, --threshold       log-likelihood difference required to make a copy number call\n"
"  -o, --outfile         name of file to write contigs to\n"
"      --help            display this help and exit\n";


namespace opt
{
	static unsigned int k;
	static unsigned int verbose;
	static double threshold;
	static std::string outfile;
}

static const char* shortopts = "k:o:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "threshold",   required_argument, NULL, 't' },
	{ "outfile",     required_argument, NULL, 'o' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


#endif
