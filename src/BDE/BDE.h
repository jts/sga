#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include "Util.h"
#include "Contig.h"
#include "Bigraph.h"
#include "PairStreamer.h"
#include "StatsCommon.h"
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
typedef std::map<ContigID, AlignPairVec> IDAPVecMap;
typedef std::map<ContigID, Contig> ContigMap;
typedef std::vector<int> IntVec;

//
// Structs
//
struct Profile
{
	Profile(ContigID i, int low, int high) : id(i), profile(low, high) {}
	ContigID id;
	IntDist profile;
};

typedef std::vector<Profile> ProfileVector;

//
// Functions
//
AlignPairVec processPairs(const AlignPairVec& apv, EdgeDir& direction, EdgeComp& orientation);
IDAPVecMap splitBlock(const AlignPairVec& block);

IntDist maxPost(int minDist, int maxDist, int sumContigLens, IntVec initFrags, const IntDist& dist, double& maxLP, int& argmax);
double intVecLogP(IntVec values, const IntDist& dist, double minP, int& discardCount);
IntDist makeConditional(const IntDist& original, int min, int max);

Contig& getContig(ContigMap& cm, ContigID id);
void writeDistEst(std::ostream& out, ContigID id0, ContigID id1, int dist, int numPairs);

// Parsing
void parseOptions(int argc, char** argv);

//
// Getopt
//
#define PROGRAM "BDE"
static const char *VERSION_MESSAGE =
PROGRAM "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION] ... CONTIGS PAIRED\n"
"Estimate distances between CONTIGS based on the PAIRED alignment file\n"
"\n"
"  -k, --kmer=KMER_SIZE    kmer-size\n"
"  -h, --histogram=FILE    file to read the fragment size histogram from\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n";


namespace opt
{
	unsigned int k;
	static std::string histFile;
	static unsigned int verbose;
}

static const char* shortopts = "o:h:k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "histogram",   required_argument, NULL, 'h' },	
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


#endif
