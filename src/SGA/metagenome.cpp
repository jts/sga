//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// metagenome - Assemble metagenomics data
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
#include "BitVector.h"
#include "MetAssembleProcess.h"
#include "metagenome.h"

// Defines to clarify awful template function calls
#define PROCESS_METASSEMBLE_SERIAL SequenceProcessFramework::processSequencesSerial<SequenceWorkItem, MetAssembleResult, \
                                                                                    MetAssemble, MetAssembleAggregateResults>
 
#define PROCESS_METASSEMBLE_PARALLEL SequenceProcessFramework::processSequencesParallel<SequenceWorkItem, MetAssembleResult, \
                                                                                        MetAssemble, MetAssembleAggregateResults>

   
//
// Getopt
//
#define SUBPROGRAM "metagenome"
static const char *METAGENOME_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *METAGENOME_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] READS\n"
"Assemble contigs from the metagenomics data in READS\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -k, --kmer=K                     use a k-mer size of size K\n"
"      -x, --kmer-threshold=T           only use kmers seen at least T times\n"
"      -t, --threads=NUM                use NUM computation threads\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"      -o, --outfile=FILE               write contigs to FILE (default: contigs.fa)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int kmer = 51;
    static int kmerThreshold = 4;
    static int sampleRate = 128;
    static int cacheLength = 10;
    static size_t minLength = 200;

    static std::string inFile;
    static std::string outFile = "contigs.fa";
}

static const char* shortopts = "o:d:k:t:x:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
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
int metagenomeMain(int argc, char** argv)
{
    Timer progTimer(PROGRAM_IDENT);

    parseMetagenomeOptions(argc, argv);

    // Create BWTS
    std::string prefix = stripFilename(opt::inFile);
    BWT* pBWT = new BWT(prefix + BWT_EXT, opt::sampleRate);
    BWT* pRevBWT = new BWT(prefix + RBWT_EXT, opt::sampleRate);
    pBWT->printInfo();

    // Create the shared bit vector and shared results aggregator
    BitVector* pSharedBitVector = new BitVector(pBWT->getBWLen());
    
    // Create interval caches to speed up k-mer lookups
    BWTIntervalCache* pBWTCache = new BWTIntervalCache(opt::cacheLength, pBWT);
    BWTIntervalCache* pRevBWTCache = new BWTIntervalCache(opt::cacheLength, pRevBWT);

    MetAssembleParameters sharedParameters;
    sharedParameters.pBWT = pBWT;
    sharedParameters.pRevBWT = pRevBWT;
    sharedParameters.pBWTCache = pBWTCache;
    sharedParameters.pRevBWTCache = pRevBWTCache;
    sharedParameters.kmer = opt::kmer;
    sharedParameters.kmerThreshold = opt::kmerThreshold;
    sharedParameters.pBitVector = pSharedBitVector;
    sharedParameters.minLength = opt::minLength;

    MetAssembleAggregateResults resultsProcess(opt::outFile);

    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode assembly\n", PROGRAM_IDENT);
        MetAssemble assembleProcess(sharedParameters); 
        PROCESS_METASSEMBLE_SERIAL(opt::inFile, &assembleProcess, &resultsProcess);
    }
    else
    {
        printf("[%s] starting parallel-mode assembly with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        
        std::vector<MetAssemble*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            MetAssemble* pProcessor = new MetAssemble(sharedParameters);
            processorVector.push_back(pProcessor);
        }
        
        PROCESS_METASSEMBLE_PARALLEL(opt::inFile, processorVector, &resultsProcess);
        
        for(size_t i = 0; i < processorVector.size(); ++i)
        {
            delete processorVector[i];
            processorVector[i] = NULL;
        }
    }

    // Cleanup
    delete pBWT;
    delete pRevBWT;
    delete pBWTCache;
    delete pRevBWTCache;
    delete pSharedBitVector;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseMetagenomeOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'k': arg >> opt::kmer; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << METAGENOME_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << METAGENOME_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // Validate parameters
    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }

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

    if (die) 
    {
        std::cout << "\n" << METAGENOME_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::inFile = argv[optind++];
}
