#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include "Util.h"
#include "Contig.h"
#include "ScaffoldData.h"
#include "Chain.h"
#include "SeqGraph.h"
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
// Structs
//


//
// Typedefs
//
typedef std::map<ContigID, ScaffoldData> SDMap;
typedef std::vector<Range> RangeVec;
typedef std::vector<ContigPosition> ContigPositionVector;

//
// Functions
//

void buildChains(SDMap& sdMap);
Contig& getContig(SDMap& sdMap, ContigID cid);

// Graph building
void addVertexToScaffoldGraph(SeqGraph& graph, ContigID id);
void addEdgeToScaffoldGraph(SeqGraph& graph, ContigID id1, ContigID id2, EdgeDir dir, EdgeComp comp, int dist);

// Parsing
void parseLinks(std::string filename, SDMap& sdMap);
void parseOptions(int argc, char** argv);

//
// Getopt
//
#define PROGRAM "Scaffold"
static const char *VERSION_MESSAGE =
PROGRAM "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION] ... CONTIGS DISTANCES\n"
"Build scaffolds from the CONTIGS based on the inferred DISTANCES.\n"
"\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n";


namespace opt
{
	static unsigned int verbose;
}

static const char* shortopts = "o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


#endif
