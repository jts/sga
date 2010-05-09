#ifndef RESOLVE_H
#define RESOLVE_H

#include "Bigraph.h"
#include "Contig.h"
#include "ContigGraph.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>

//
// Typedefs
//


//
// Functions
//
void parseOptions(int argc, char** argv);
std::string colorVertexByLength(Contig d);
std::string colorVertexByAnnotation(Contig d);

//
// Getopt
//
#define PROGRAM "Resolve"
static const char *VERSION_MESSAGE =
PROGRAM "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION] ... CONTIGS ADJACENCY\n"
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
