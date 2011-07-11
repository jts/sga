//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// graph-diff - Find strings that are only present
// in one of two input graphs
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "SGACommon.h"
#include "GraphCompare.h"
#include "graph-diff.h"

// Defines to clarify awful template function calls
#define PROCESS_GDIFF_SERIAL SequenceProcessFramework::processSequencesSerial<SequenceWorkItem, GraphCompareResult, \
                                                                              GraphCompare, GraphCompareAggregateResults>

#define PROCESS_GDIFF_PARALLEL SequenceProcessFramework::processSequencesParallel<SequenceWorkItem, GraphCompareResult, \
                                                                                  GraphCompare, GraphCompareAggregateResults>

   
//
// Getopt
//
#define SUBPROGRAM "graph-diff"
static const char *GRAPH_DIFF_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *GRAPH_DIFF_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] --base BASE.fa --variant VARIANT.fa\n"
"Find and report strings only present in the graph of VARIANT when compared to BASE\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -b, --base=FILE                  the baseline reads are in FILE\n"
"      -r, --variant=FILE               the variant reads are in FILE\n"
"      -o, --outfile=FILE               write the strings found to FILE\n"
"      -k, --kmer=K                     use K as the k-mer size for variant discovery\n"
"      -x, --kmer-threshold=T           only used kmers seen at least T times\n"
"      -y, --max-branches=B             allow the search process to branch B times when \n"
"                                       searching for the completion of a bubble (default: 0)\n"
"      -t, --threads=NUM                use NUM computation threads\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int kmer = 55;
    static int kmerThreshold = 2;
    static int maxBranches = 0;
    static int sampleRate = 128;
    static std::string baseFile;
    static std::string variantFile;
    static std::string outFile = "variants.fa";
}

static const char* shortopts = "b:r:o:k:t:x:y:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "base",          required_argument, NULL, 'b' },
    { "variants",      required_argument, NULL, 'r' },
    { "kmer",          required_argument, NULL, 'k' },
    { "kmer-threshold",required_argument, NULL, 'x' },
    { "max-branches",  required_argument, NULL, 'y' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int graphDiffMain(int argc, char** argv)
{
    parseGraphDiffOptions(argc, argv);

    // Create BWTS
    std::string basePrefix = stripFilename(opt::baseFile);
    BWT* pBaseBWT = new BWT(basePrefix + BWT_EXT, opt::sampleRate);
    BWT* pBaseRBWT = new BWT(basePrefix + RBWT_EXT, opt::sampleRate);

    std::string variantPrefix = stripFilename(opt::variantFile);
    BWT* pVariantBWT = new BWT(variantPrefix + BWT_EXT, opt::sampleRate);
    BWT* pVariantRBWT = new BWT(variantPrefix + RBWT_EXT, opt::sampleRate);
    
    // Create the shared bit vector and shared results aggregator
    BitVector* pSharedBitVector = new BitVector(pVariantBWT->getBWLen());
    GraphCompareAggregateResults* pSharedResults = new GraphCompareAggregateResults();

    // Set the parameters shared between all threads
    GraphCompareParameters sharedParameters;
    sharedParameters.pVariantBWT = pVariantBWT;
    sharedParameters.pVariantRevBWT = pVariantRBWT;
    sharedParameters.pBaseBWT = pBaseBWT;
    sharedParameters.pBaseRevBWT = pBaseRBWT;
    sharedParameters.kmer = opt::kmer;
    sharedParameters.pBitVector = pSharedBitVector;
    sharedParameters.kmerThreshold = 3;
    sharedParameters.maxBranches = opt::maxBranches;

    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode graph diff\n", PROGRAM_IDENT);
        GraphCompare graphCompare(sharedParameters); 
        PROCESS_GDIFF_SERIAL(opt::variantFile, &graphCompare, pSharedResults);
        graphCompare.updateSharedStats(pSharedResults);
    }
    else
    {
        printf("[%s] starting parallel-mode graph diff with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        
        std::vector<GraphCompare*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            GraphCompare* pProcessor = new GraphCompare(sharedParameters);
            processorVector.push_back(pProcessor);
        }
        
        PROCESS_GDIFF_PARALLEL(opt::variantFile, processorVector, pSharedResults);
        
        for(size_t i = 0; i < processorVector.size(); ++i)
        {
            // Update the shared stats
            processorVector[i]->updateSharedStats(pSharedResults);

            delete processorVector[i];
            processorVector[i] = NULL;
        }


    }
    pSharedResults->printStats();

    // Cleanup
    delete pBaseBWT;
    delete pBaseRBWT;
    delete pVariantBWT;
    delete pVariantRBWT;
    delete pSharedBitVector;
    delete pSharedResults;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseGraphDiffOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'k': arg >> opt::kmer; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case 'b': arg >> opt::baseFile; break;
            case 'r': arg >> opt::variantFile; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'y': arg >> opt::maxBranches; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GRAPH_DIFF_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GRAPH_DIFF_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // Validate parameters
    if (argc - optind > 1) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if(opt::baseFile.empty() || opt::variantFile.empty())
    {
        std::cerr << SUBPROGRAM ": error a --base and --variant file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << GRAPH_DIFF_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}
