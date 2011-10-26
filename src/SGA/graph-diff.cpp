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
#include "SampledSuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "SGACommon.h"
#include "GraphCompare.h"
#include "VCFTester.h"
#include "DindelRealignWindow.h"
#include "graph-diff.h"

// Functions
void runVCFTester(GraphCompareParameters& parameters);
void runGraphDiff(GraphCompareParameters& parameters);

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
"          --reference=FILE             the reference FILE\n"
"      -o, --outfile=FILE               write the strings found to FILE\n"
"      -p, --prefix=NAME                prefix the output files with NAME\n"
"      -k, --kmer=K                     use K as the k-mer size for variant discovery\n"
"      -x, --kmer-threshold=T           only used kmers seen at least T times\n"
"      -y, --max-branches=B             allow the search process to branch B times when \n"
"                                       searching for the completion of a bubble (default: 0)\n"
"      -t, --threads=NUM                use NUM computation threads\n"
"          --test=VCF                   test the variants in the provided VCF file\n"
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
    static int cacheLength = 10;

    static std::string outPrefix = "graphdiff";
    static std::string referenceFile;
    static std::string baseFile;
    static std::string variantFile;
    static std::string inputVCFFile;
//    static std::string outFile = "variants.fa";
}

static const char* shortopts = "b:r:o:k:d:t:x:y:p:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_REFERENCE, OPT_TESTVCF };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "base",          required_argument, NULL, 'b' },
    { "variants",      required_argument, NULL, 'r' },
    { "outfile",       required_argument, NULL, 'o' },
    { "kmer",          required_argument, NULL, 'k' },
    { "kmer-threshold",required_argument, NULL, 'x' },
    { "max-branches",  required_argument, NULL, 'y' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "prefix",        required_argument, NULL, 'p' },
    { "references",    required_argument, NULL, OPT_REFERENCE },
    { "test"      ,    required_argument, NULL, OPT_TESTVCF },
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

    // Create indices for the base reads
    std::string basePrefix = stripFilename(opt::baseFile);
    BWT* pBaseBWT = new BWT(basePrefix + BWT_EXT, opt::sampleRate);
    SampledSuffixArray* pBaseSSA = new SampledSuffixArray(basePrefix + SAI_EXT, SSA_FT_SAI);

    // Create indices for the variant reads
    std::string variantPrefix = stripFilename(opt::variantFile);
    BWT* pVariantBWT = new BWT(variantPrefix + BWT_EXT, opt::sampleRate);
    SampledSuffixArray* pVariantSSA = new SampledSuffixArray(variantPrefix + SAI_EXT, SSA_FT_SAI);
    
    std::cout << "Variant index memory info\n";
    pVariantBWT->printInfo();
    pVariantSSA->printInfo();

    // Create indices for the reference
    std::string refPrefix = stripFilename(opt::referenceFile);
    BWT* pRefBWT = new BWT(refPrefix + BWT_EXT, opt::sampleRate);
    BWT* pRefRevBWT = new BWT(refPrefix + RBWT_EXT, opt::sampleRate);
    SampledSuffixArray* pRefSSA = new SampledSuffixArray(refPrefix + SSA_EXT);

    // Read in the reference 
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    // Create interval caches to speed up k-mer lookups
    BWTIntervalCache varBWTCache(opt::cacheLength, pVariantBWT);

    BWTIntervalCache baseBWTCache(opt::cacheLength, pBaseBWT);

    // Set the parameters shared between all threads
    GraphCompareParameters sharedParameters;

    sharedParameters.pBaseBWT = pBaseBWT;
    sharedParameters.pBaseBWTCache = &baseBWTCache;
    sharedParameters.pBaseSSA = pBaseSSA;

    sharedParameters.pVariantBWT = pVariantBWT;
    sharedParameters.pVariantBWTCache = &varBWTCache;
    sharedParameters.pVariantSSA = pVariantSSA;

    sharedParameters.pReferenceBWT = pRefBWT;
    sharedParameters.pReferenceRevBWT = pRefRevBWT;
    sharedParameters.pReferenceSSA = pRefSSA;
    sharedParameters.pRefTable = &refTable;

    sharedParameters.kmer = opt::kmer;
    sharedParameters.pBitVector = NULL;
    sharedParameters.kmerThreshold = opt::kmerThreshold;
    sharedParameters.maxBranches = opt::maxBranches;

    // If a VCF file was provided, just test the variants in the file without any discovery
    if(!opt::inputVCFFile.empty())
        runVCFTester(sharedParameters);
    else
        runGraphDiff(sharedParameters);

    // Cleanup
    delete pBaseBWT;
    delete pBaseSSA;

    delete pVariantBWT;
    delete pVariantSSA;

    delete pRefBWT;
    delete pRefSSA;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

void runVCFTester(GraphCompareParameters& parameters)
{
    (void)parameters;
    std::cout << "Testing variants in " << opt::inputVCFFile << "\n";
    std::string line;

    VCFTester tester(parameters);
    try
    {
        VCFFile input(opt::inputVCFFile, "r");
        VCFFile::VCFEntry record;

        while(1)
        {
            record = input.getNextEntry();
            if(record.isEmpty())
                break;
            else
                tester.process(record);
        }
    }
    catch(std::string e)
    {
        std::cerr << "Exception: " << e << "\n";
        exit(EXIT_FAILURE);
    }
}

void runGraphDiff(GraphCompareParameters& parameters)
{
    // Create the shared bit vector and shared results aggregator
    BitVector* pSharedBitVector = new BitVector(parameters.pVariantBWT->getBWLen());
    GraphCompareAggregateResults* pSharedResults = new GraphCompareAggregateResults(opt::outPrefix);
    
    // Set the bit vector
    parameters.pBitVector = pSharedBitVector;

    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode graph diff\n", PROGRAM_IDENT);
        GraphCompare graphCompare(parameters); 
        PROCESS_GDIFF_SERIAL(opt::variantFile, &graphCompare, pSharedResults);
        graphCompare.updateSharedStats(pSharedResults);
    }
    else
    {
        printf("[%s] starting parallel-mode graph diff with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        
        std::vector<GraphCompare*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            GraphCompare* pProcessor = new GraphCompare(parameters);
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
    
    delete pSharedBitVector;
    parameters.pBitVector = NULL;

    delete pSharedResults;
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
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
//            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'y': arg >> opt::maxBranches; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'p': arg >> opt::outPrefix; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_TESTVCF: arg >> opt::inputVCFFile; break;
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

    if(opt::referenceFile.empty())
    {
        std::cerr << SUBPROGRAM ": error, a --reference file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << GRAPH_DIFF_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}
