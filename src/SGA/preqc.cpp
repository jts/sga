//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// pre-QC - Perform pre-assembly quality checks on a set
//          of reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "preqc.h"
#include "Timer.h"
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SGACommon.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "preqc"
static const char *PREQC_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *PREQC_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Perform pre-assembly quality checks\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string prefix;
    static std::string readsFile;
}

static const char* shortopts = "p:d:t:o:k:n:b:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_RUNLENGTHS, OPT_KMERDIST, OPT_NOOVERLAP};

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "threads",            required_argument, NULL, 't' },
    { "prefix",             required_argument, NULL, 'p' },
    { "sample-rate",        required_argument, NULL, 'd' },
    { "kmer-size",          required_argument, NULL, 'k' },
    { "num-reads",          required_argument, NULL, 'n' },
    { "branch-cutoff",      required_argument, NULL, 'b' },
    { "kmer-distribution",  no_argument,       NULL, OPT_KMERDIST },
    { "no-overlap",         no_argument,       NULL, OPT_NOOVERLAP },
    { "run-lengths",        no_argument,       NULL, OPT_RUNLENGTHS },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int preQCMain(int argc, char** argv)
{
    parsePreQCOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    printf("Loading FM-index of %s\n", opt::readsFile.c_str());
    BWTIndexSet index_set;
    index_set.pBWT = new BWT(opt::prefix + BWT_EXT);

    delete index_set.pBWT;
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parsePreQCOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << PREQC_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << PREQC_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 1) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << PREQC_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
        opt::prefix = stripFilename(opt::readsFile);
}
