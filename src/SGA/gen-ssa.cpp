//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// gen-ssa - Build a sampled suffix array for a set of reads
// using its BWT
//
#include <iostream>
#include <fstream>
#include "SGACommon.h"
#include "Util.h"
#include "gen-ssa.h"
#include "SuffixArray.h"
#include "SeqReader.h"
#include "SACAInducedCopying.h"
#include "SampledSuffixArray.h"
#include "BWTDiskConstruction.h"
#include "BWT.h"
#include "Timer.h"

//
// Getopt
//
#define SUBPROGRAM "gen-ssa"

static const char *GENSSA_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *GENSSA_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Build a sampled suffix array for the reads in READSFILE using the BWT\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"  -t, --threads=NUM                    use NUM threads to construct the index (default: 1)\n"
"  -c, --check                          validate that the suffix array/bwt is correct\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string readsFile;
    static std::string prefix;
    static int numThreads = 1;
    static int sampleRate = 32;
    static bool validate = false;
}

static const char* shortopts = "t:cv";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "check",       no_argument,       NULL, 'c' },
    { "threads",     required_argument, NULL, 't' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int genSSAMain(int argc, char** argv)
{
    Timer t("sga gen-ssa");
    parseGenSSAOptions(argc, argv);
    
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    ReadInfoTable* pRIT = new ReadInfoTable(opt::readsFile, pBWT->getNumStrings(), RIO_NUMERICID);
    pBWT->printInfo();

    SampledSuffixArray* pSSA = new SampledSuffixArray();
    pSSA->build(pBWT, pRIT, opt::sampleRate);
    pSSA->printInfo();
    pSSA->writeSSA(opt::prefix + SSA_EXT);

    if(opt::validate)
        pSSA->validate(opt::readsFile, pBWT);

    delete pBWT;
    delete pRIT;
    delete pSSA;

    return 0;
}

// 
// Handle command line arguments
//
void parseGenSSAOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'c': opt::validate = true; break;
            case 't': arg >> opt::numThreads; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GENSSA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GENSSA_VERSION_MESSAGE;
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
        std::cout << "\n" << GENSSA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];
    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }
}
