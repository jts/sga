//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// haplotype-filter - perform quality checks
// on haplotypes and their associated variants
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "haplotype-filter.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGACommon.h"
#include "Timer.h"
#include "SequenceProcessFramework.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "haplotype-filter"
static const char *HAPLOTYPE_FILTER_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *HAPLOTYPE_FILTER_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Remove haplotypes and their associated variants from a data set.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -r, --reads=FILE                 load the FM-index of the reads in FILE\n"
"      -o, --out-prefix=STR             write the passed haplotypes and variants to STR.vcf and STR.fa\n" 
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string readsFile;
    static std::string outFile = "hapfilter";
}

static const char* shortopts = "i:o:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "threads",               required_argument, NULL, 't' },
    { "outfile",               required_argument, NULL, 'o' },
    { "reads",                 required_argument, NULL, 'r' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int haplotypeFilterMain(int argc, char** argv)
{
    parseHaplotypeFilterOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    std::string prefix = stripGzippedExtension(opt::readsFile);
    BWT* pBWT = new BWT(prefix + BWT_EXT);
    delete pBWT;

    // Cleanup
    delete pTimer;
    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseHaplotypeFilterOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'r': arg >> opt::readsFile; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << HAPLOTYPE_FILTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << HAPLOTYPE_FILTER_VERSION_MESSAGE;
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
        std::cout << "\n" << HAPLOTYPE_FILTER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
