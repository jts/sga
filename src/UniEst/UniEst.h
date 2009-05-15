#ifndef UNIEST_H
#define UNIEST_H

#include "Util.h"
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
typedef std::vector<int> IntVec;
typedef std::vector<double> DoubleVec;
typedef std::set<std::string> StringSet;

struct LenMean
{
	size_t length;
	double mean;

	static bool sortByLen(const LenMean& lm1, const LenMean& lm2)
	{
		return lm1.length < lm2.length;
	}

	static bool sortByMean(const LenMean& lm1, const LenMean& lm2)
	{
		return lm1.mean < lm2.mean;
	}
	

	friend std::ostream& operator<< (std::ostream& out, LenMean& lm)
	{
		out << lm.length << "\t" << lm.mean;
		return out;
	}
};

struct UniData
{
	IntVec kmerCoverage;
	StringSet pairedCoverage[2];
};

typedef std::map<ContigID, UniData> CoverageMap;


//
// Functions
//
void addCoverage(CoverageMap& cov_map, const KAlignment ka);
void addOverhangCoverage(CoverageMap& covMap, const KAlignment& ka1, const KAlignment& ka2);
double getMeanCoverage(const IntVec& iVec);

double estimateGlobalMeanCoverage(CoverageMap& cov_map);
void estimateUniquenessByDepth(CoverageMap& cov_map, double mean_param);
void estimateUniquenessByOverhang(CoverageMap& covMap);

void parseOptions(int argc, char** argv);
void printCoverage(CoverageMap& cov_map, ContigID id);
double log_poisson(int k, double m);
double log_factorial(int k);


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
"      --help            display this help and exit\n";


namespace opt
{
	static unsigned int k;
	static unsigned int verbose;
}

static const char* shortopts = "k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


#endif
