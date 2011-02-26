//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// cluster - extract connected components from an asqg file
//
#include <iostream>
#include <fstream>
#include "cluster.h"
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGPairedAlgorithms.h"
#include "SGDebugAlgorithms.h"
#include "SGVisitors.h"
#include "SGSearch.h"
#include "Timer.h"
#include "OverlapCommon.h"
#include "BitVector.h"
#include "ClusterProcess.h"

//
// Getopt
//
#define SUBPROGRAM "cluster"

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

static const char *CLUSTER_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *CLUSTER_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] READS\n"
"Extract clusters of reads belonging to the same connected components.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the clusters to FILE (default: clusters.txt)\n"
"      -c, --cluster-size=N             only write clusters with at least N reads (default: 2)\n"
"      -m, --min-overlap=N              require an overlap of at least N bases between reads (default: 45)\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences aligned (default: exact matches only)\n"
"      -t, --threads=NUM                use NUM worker threads to compute the overlaps (default: no threading)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static size_t minSize = 2;
    static std::string outFile;
    static std::string readsFile;
    static std::string prefix;
    static int seedLength = 16;
    static int seedStride = seedLength;
    static int numThreads = 1;
    static double errorRate = 0.0f;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
}

static const char* shortopts = "o:m:c:t:e:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "out",            required_argument, NULL, 'o' },
    { "cluster-size",   required_argument, NULL, 'c' },
    { "min-overlap",    required_argument, NULL, 'm' },
    { "error-rate",     required_argument, NULL, 'e' },
    { "threads",        required_argument, NULL, 't' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int clusterMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga cluster");
    parseClusterOptions(argc, argv);
    cluster();
    delete pTimer;

    return 0;
}

void cluster()
{
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT,opt::errorRate, opt::seedLength, opt::seedStride, true);

    BitVector markedReads(pBWT->getNumStrings());

    std::string preclustersFile = opt::prefix + ".preclusters";
    std::ostream* pPreWriter = createWriter(preclustersFile);
    ClusterPostProcess postProcessor(pPreWriter, opt::minSize, &markedReads);
    
    // Make pre-clusters from the reads
    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode read clustering\n", PROGRAM_IDENT);
        ClusterProcess processor(pOverlapper, opt::minOverlap, &markedReads);
        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                         ClusterResult, 
                                                         ClusterProcess, 
                                                         ClusterPostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        printf("[%s] starting parallel-mode read clustering computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        
        std::vector<ClusterProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            ClusterProcess* pProcessor = new ClusterProcess(pOverlapper, opt::minOverlap, &markedReads);
            processorVector.push_back(pProcessor);
        }

        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                         ClusterResult, 
                                                         ClusterProcess, 
                                                         ClusterPostProcess>(opt::readsFile, processorVector, &postProcessor);
        
        for(size_t i = 0; i < processorVector.size(); ++i)
        {
            delete processorVector[i];
            processorVector[i] = NULL;
        }
    }
    delete pPreWriter;
    delete pBWT;
    delete pRBWT;
    delete pOverlapper;

    // Open the preclusters file and convert them to read names
    SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
    ReadInfoTable* pRIT = new ReadInfoTable(opt::readsFile, pFwdSAI->getNumStrings());

    std::istream* pPreReader = createReader(preclustersFile);
    std::ostream* pClusterWriter = createWriter(opt::outFile);
    std::string line;
    while(getline(*pPreReader,line))
    {
        std::stringstream parser(line);
        std::string clusterName;
        std::string readSequence;
        size_t clusterSize;
        int64_t lowIdx;
        int64_t highIdx;
        parser >> clusterName >> clusterSize >> readSequence >> lowIdx >> highIdx;
        assert(lowIdx <= highIdx);

        for(int64_t i = lowIdx; i <= highIdx; ++i)
        {
            const ReadInfo& targetInfo = pRIT->getReadInfo(pFwdSAI->get(i).getID());
            std::string readName = targetInfo.id;
            *pClusterWriter << clusterName << "\t" << clusterSize << "\t" << readName << "\t" << readSequence << "\n";
        }
    }
    unlink(preclustersFile.c_str());

    delete pFwdSAI;
    delete pRIT;
    delete pPreReader;
    delete pClusterWriter;
}

// 
// Handle command line arguments
//
void parseClusterOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'c': arg >> opt::minSize; break;
            case 'e': arg >> opt::errorRate; break;
            case 'm': arg >> opt::minOverlap; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << CLUSTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CLUSTER_VERSION_MESSAGE;
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

    if (die) 
    {
        std::cerr << "Try `" << SUBPROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filename
    opt::readsFile = argv[optind++];
    opt::prefix = stripFilename(opt::readsFile);

    if(opt::outFile.empty())
        opt::outFile = opt::prefix + ".clusters";
}
