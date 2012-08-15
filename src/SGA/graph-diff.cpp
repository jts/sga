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
void runDebug(GraphCompareParameters& parameters);

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
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] --base BASE.fa --variant VARIANT.fa --ref REFERENCE.fa\n"
"Find and report strings only present in the graph of VARIANT when compared to BASE\n"
"\n"
"General options:\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=NAME                prefix the output files with NAME\n"
"      -t, --threads=NUM                use NUM computation threads\n"
//"          --test=VCF                   test the variants in the provided VCF file\n"
"\n"
"Index options:\n"
"      -r, --variant=FILE               call variants present in the read set in FILE\n"
"      -b, --base=FILE                  use the read set in FILE as the base line for comparison\n"
"                                       if this option is not given, reference-based calls will be made\n"
"          --reference=FILE             use the reference sequence in FILE\n"
"\n"
"Algorithm options:\n"
"      -k, --kmer=K                     use K-mers to discover variants\n"
"      -x, --min-discovery-count=T      require a variant k-mer to be seen at least T times\n"
"          --debruijn                   use the de Bruijn graph assembly algorithm (default: string graph)\n"
"      -m, --min-overlap=N              require at least N bp overlap when assembling using a string graph\n" 
"          --min-dbg-count=T            only use k-mers seen T times when assembling using a de Bruijn graph\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int cacheLength = 10;
    static int sampleRate = 128;
    
    static int kmer = 55;
    static int minDiscoveryCount = 2;
    static int maxDiscoveryCount = 400;
    static int minOverlap = 61;

    static bool deBruijnMode = false;
    static int minDBGCount = 2;
    static bool lowCoverage = false;

    static bool referenceMode = false;
    static std::string outPrefix = "graphdiff";
    static std::string indexPrefix;
    static std::string debugFile;
    static std::string referenceFile;
    static std::string baseFile;
    static std::string variantFile;
    static std::string inputVCFFile;
}

static const char* shortopts = "b:r:o:k:d:t:x:y:p:m:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_REFERENCE, OPT_TESTVCF, OPT_DEBUG, OPT_MIN_DBG_COUNT, OPT_INDEX, OPT_DEBRUIJN, OPT_LOWCOVERAGE };

static const struct option longopts[] = {
    { "verbose",              no_argument,       NULL, 'v' },
    { "threads",              required_argument, NULL, 't' },
    { "base",                 required_argument, NULL, 'b' },
    { "variant",              required_argument, NULL, 'r' },
    { "kmer",                 required_argument, NULL, 'k' },
    { "min-discovery-count",  required_argument, NULL, 'x' },
    { "sample-rate",          required_argument, NULL, 'd' },
    { "prefix",               required_argument, NULL, 'p' },
    { "min-overlap",          required_argument, NULL, 'm' },
    { "debruijn",             no_argument,       NULL, OPT_DEBRUIJN },
    { "low-coverage",         no_argument,       NULL, OPT_LOWCOVERAGE },
    { "index",                required_argument, NULL, OPT_INDEX },
    { "min-dbg-count",        required_argument, NULL, OPT_MIN_DBG_COUNT },
    { "debug",                required_argument, NULL, OPT_DEBUG },
    { "reference",            required_argument, NULL, OPT_REFERENCE },
    { "test"      ,           required_argument, NULL, OPT_TESTVCF },
    { "help",                 no_argument,       NULL, OPT_HELP },
    { "version",              no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int graphDiffMain(int argc, char** argv)
{
    parseGraphDiffOptions(argc, argv);

    // Create indices for the variant reads
    std::string variantPrefix = stripGzippedExtension(opt::variantFile);

    // Use debug index prefix if specified
    if(!opt::indexPrefix.empty())
        variantPrefix = opt::indexPrefix;

    //
    // Initialize indices
    //

    // Variant reads
    BWTIndexSet variantIndex;
    variantIndex.pBWT = new BWT(variantPrefix + BWT_EXT, opt::sampleRate);
    variantIndex.pSSA = new SampledSuffixArray(variantPrefix + SAI_EXT, SSA_FT_SAI);
    variantIndex.pCache = new BWTIntervalCache(opt::cacheLength, variantIndex.pBWT);
    if(opt::lowCoverage)
        variantIndex.pPopIdx = new PopulationIndex(variantPrefix + POPIDX_EXT);

    // Reference genome
    BWTIndexSet referenceIndex;
    std::string refPrefix = stripExtension(opt::referenceFile);

    referenceIndex.pBWT = new BWT(refPrefix + BWT_EXT, opt::sampleRate);
    referenceIndex.pSSA = new SampledSuffixArray(refPrefix + SSA_EXT);

    // Base reads, only if not in reference mode
    BWTIndexSet baseIndex;

    if(!opt::referenceMode)
    {
        std::string basePrefix = stripExtension(opt::baseFile);
        baseIndex.pBWT = new BWT(basePrefix + BWT_EXT, opt::sampleRate);
        baseIndex.pSSA = new SampledSuffixArray(basePrefix + SAI_EXT, SSA_FT_SAI);
        baseIndex.pCache = new BWTIntervalCache(opt::cacheLength, baseIndex.pBWT);
    }
    else
    {
        baseIndex.pBWT = referenceIndex.pBWT;
        baseIndex.pSSA = referenceIndex.pSSA;
        baseIndex.pCache = NULL;
    }

    // 
    std::cout << "Base index memory info\n";
    baseIndex.pBWT->printInfo();
    baseIndex.pSSA->printInfo();

    //
    std::cout << "Variant index memory info\n";
    variantIndex.pBWT->printInfo();
    variantIndex.pSSA->printInfo();

    //
    std::cout << "Reference index memory info\n";
    referenceIndex.pBWT->printInfo();
    referenceIndex.pSSA->printInfo();

    // Read in the reference 
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    // Set the parameters shared between all threads
    GraphCompareParameters sharedParameters;

    sharedParameters.baseIndex = baseIndex;
    sharedParameters.variantIndex = variantIndex;
    sharedParameters.referenceIndex = referenceIndex;
    sharedParameters.pRefTable = &refTable;
    sharedParameters.bReferenceMode = opt::referenceMode;

    sharedParameters.algorithm = opt::deBruijnMode ? GCA_DEBRUIJN_GRAPH : GCA_STRING_GRAPH;
    sharedParameters.kmer = opt::kmer;
    sharedParameters.pBitVector = NULL;
    sharedParameters.minDiscoveryCount = opt::minDiscoveryCount;
    sharedParameters.maxDiscoveryCount = opt::maxDiscoveryCount;
    sharedParameters.minDBGCount = opt::minDBGCount;
    sharedParameters.minOverlap = opt::minOverlap;

    if (opt::lowCoverage)
        sharedParameters.dindelRealignParameters.multiSample = 1;

    if(!opt::debugFile.empty())
    {
        runDebug(sharedParameters);
    }
    else
    {
        // If a VCF file was provided, just test the variants in the file without any discovery
        if(!opt::inputVCFFile.empty())
            runVCFTester(sharedParameters);
        else
            runGraphDiff(sharedParameters);
    }

    // Cleanup
    if(!opt::referenceMode)
    {
        delete baseIndex.pBWT;
        delete baseIndex.pSSA;
        delete baseIndex.pCache;
    }

    // Cleanup indices
    delete variantIndex.pBWT;
    delete variantIndex.pSSA;
    delete variantIndex.pCache;
    if(opt::lowCoverage)
        delete variantIndex.pPopIdx;

    delete referenceIndex.pBWT;
    delete referenceIndex.pSSA;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

void runVCFTester(GraphCompareParameters& parameters)
{
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
    BitVector* pSharedBitVector = new BitVector(parameters.variantIndex.pBWT->getBWLen());

    // This call can throw via dindel
    GraphCompareAggregateResults* pSharedResults;
    try {
        StringVector samples;

        // If in multi-sample mode, write sample names in the VCF header
        if(parameters.variantIndex.pPopIdx != NULL)
            samples = parameters.variantIndex.pPopIdx->getSamples();
        pSharedResults = new GraphCompareAggregateResults(opt::outPrefix, samples);
    }
    catch(std::string e)
    {
        std::cout << "Exception: " << e << "\n";
        exit(EXIT_FAILURE);
    }

    std::cout << "Bit vector capacity: " << pSharedBitVector->capacity() << "\n";

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

// Run in debug mode
void runDebug(GraphCompareParameters& parameters)
{
    // Create the shared bit vector and shared results aggregator
    BitVector* pSharedBitVector = new BitVector(parameters.variantIndex.pBWT->getBWLen());

    // Set the bit vector
    parameters.pBitVector = pSharedBitVector;

    GraphCompare graphCompare(parameters); 
    graphCompare.testKmersFromFile(opt::debugFile);
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
            case 'x': arg >> opt::minDiscoveryCount; break;
            case 'b': arg >> opt::baseFile; break;
            case 'r': arg >> opt::variantFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'p': arg >> opt::outPrefix; break;
            case 'm': arg >> opt::minOverlap; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case OPT_DEBRUIJN: opt::deBruijnMode = true; break;
            case OPT_LOWCOVERAGE: opt::lowCoverage = true; break;
            case OPT_MIN_DBG_COUNT: arg >> opt::minDBGCount; break;
            case OPT_DEBUG: arg >> opt::debugFile; break;
            case OPT_TESTVCF: arg >> opt::inputVCFFile; break;
            case OPT_INDEX: arg >> opt::indexPrefix; break;
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

    if(opt::variantFile.empty())
    {
        std::cerr << SUBPROGRAM ": error a --base and --variant file must be provided\n";
        die = true;
    }

    if(opt::baseFile.empty())
    {
        std::cerr << SUBPROGRAM ": reference-based calling enabled\n";
        opt::referenceMode = true;
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
