//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// hapgen - Generate candidate haplotypes from
// an assembly graph
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
#include "HapgenProcess.h"
#include "hapgen.h"

// Defines to clarify awful template function calls
#define PROCESS_HAPGEN_SERIAL SequenceProcessFramework::processSequencesSerial<SequenceWorkItem, GraphCompareResult, \
                                                                               GraphCompare, GraphCompareAggregateResults>

#define PROCESS_HAPGEN_PARALLEL SequenceProcessFramework::processSequencesParallel<SequenceWorkItem, GraphCompareResult, \
                                                                                  GraphCompare, GraphCompareAggregateResults>

   
//
// Getopt
//
#define SUBPROGRAM "hapgen"
static const char *HAPGEN_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *HAPGEN_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] --ref REF.fa --sites SITES.txt READS.fa\n"
"Generate candidate haplotypes around the SITES using the assembly graph of READS\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -r, --reference=FILE             the reference genome to use\n"
"      -s, --sites=FILE                 the coordinates on the reference to generate haplotypes for\n"
"      -o, --outfile=FILE               write results to VCF FILE\n"
"      -k, --kmer=K                     use K as the k-mer size for variant discovery\n"
"      -x, --kmer-threshold=T           only used kmers seen at least T times\n"
"      -t, --threads=NUM                use NUM computation threads\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

//static const char* PROGRAM_IDENT = PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int kmer = 55;
    static int kmerThreshold = 3;
    static int sampleRate = 128;
    static int cacheLength = 10;

    static std::string referenceFile;
    static std::string sitesFile;
    static std::string readsFile;
    static std::string outFile;
}

static const char* shortopts = "o:k:t:r:s:d:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "reference",     required_argument, NULL, 'r' },
    { "sites",         required_argument, NULL, 's' },
    { "outfile",       required_argument, NULL, 'o' },
    { "kmer",          required_argument, NULL, 'k' },
    { "kmer-threshold",required_argument, NULL, 'x' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int hapgenMain(int argc, char** argv)
{
    
    parseHapgenOptions(argc, argv);

    // In the BWTs and create interval caches
    std::string basePrefix = stripFilename(opt::readsFile);
    BWT* pBWT = new BWT(basePrefix + BWT_EXT, opt::sampleRate);
    BWT* pRevBWT = new BWT(basePrefix + RBWT_EXT, opt::sampleRate);
    SampledSuffixArray* pSSA = new SampledSuffixArray(basePrefix + SAI_EXT, SSA_FT_SAI);

    pBWT->printInfo();
    pSSA->printInfo();

    BWTIntervalCache* pBWTCache = new BWTIntervalCache(opt::cacheLength, pBWT);
    BWTIntervalCache* pRevBWTCache = new BWTIntervalCache(opt::cacheLength, pRevBWT);

    // Read in a table of the reference genome
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    HapgenParameters parameters;
    parameters.pBWT = pBWT;
    parameters.pRevBWT = pRevBWT;
    parameters.pSSA = pSSA;
    parameters.pBWTCache = pBWTCache;
    parameters.pRevBWTCache = pRevBWTCache;
    parameters.kmer = opt::kmer;
    parameters.kmerThreshold = opt::kmerThreshold;
    parameters.pRefTable = &refTable;
    parameters.verbose = opt::verbose;
    parameters.vcfOutfile = opt::outFile;

    HapgenProcess processor(parameters);

    std::cout << "Parsing file\n";

    std::istream* pReader = createReader(opt::sitesFile);
    std::string line;

    while(getline(*pReader, line))
    {
        std::stringstream parser(line);
        std::string refName;
        std::string comment;
        size_t start;
        size_t end;
        parser >> refName >> start >> end >> comment;
        try {
            processor.processSite(refName, start, end, comment);
        } catch (std::string errorString)
        {
            std::cerr << "String exception: " << errorString << std::endl;
        }
    }


    delete pReader;

    // Cleanup
    delete pBWT;
    delete pRevBWT;
    delete pBWTCache;
    delete pRevBWTCache;
    delete pSSA;

        

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseHapgenOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'k': arg >> opt::kmer; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case 's': arg >> opt::sitesFile; break;
            case 'r': arg >> opt::referenceFile; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << HAPGEN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << HAPGEN_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // Validate parameters
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

    if(opt::referenceFile.empty() || opt::sitesFile.empty() || opt::outFile.empty())
    {
        std::cerr << SUBPROGRAM ": error a --reference and --sites  and --outfile file must be provided\n";
        die = true;
    }

    

    if (die) 
    {
        std::cout << "\n" << HAPGEN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::readsFile = argv[optind++];
}
