//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// connect - Determine the complete sequence of a 
// paired end fragment by finding a walk that
// connects the ends.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "connect.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "ConnectProcess.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "connect"
static const char *CONNECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *CONNECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Resolve the complete sequence of a paired end fragment by finding a walk through the graph connecting the ends\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -e, --error-rate                 the maximum error rate allowed between two sequences to consider them aligned (default: exact matches only)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 45)\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -o, --outfile=FILE               write the connected reads to FILE\n"
"      -l, --seed-length=LEN            force the seed length to be LEN. By default, the seed length in the overlap step\n"
"                                       is calculated to guarantee all overlaps with --error-rate differences are found.\n"
"                                       This option removes the guarantee but will be (much) faster. As SGA can tolerate some\n"
"                                       missing edges, this option may be preferable for some data sets.\n"
"      -s, --seed-stride=LEN            force the seed stride to be LEN. This parameter will be ignored unless --seed-length\n"
"                                       is specified (see above). This parameter defaults to the same value as --seed-length\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string prefix;
    static std::string readsFile;
    static std::string outFile;

    static double errorRate = 0.0f;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
    static int seedLength = 0;
    static int seedStride = 0;
}

static const char* shortopts = "p:m:e:t:l:s:o:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "threads",     required_argument, NULL, 't' },
    { "min-overlap", required_argument, NULL, 'm' },
    { "outfile",     required_argument, NULL, 'o' },
    { "prefix",      required_argument, NULL, 'p' },
    { "error-rate",  required_argument, NULL, 'e' },
    { "seed-length", required_argument, NULL, 'l' },
    { "seed-stride", required_argument, NULL, 's' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int connectMain(int argc, char** argv)
{
    parseConnectOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, 
                                                         opt::errorRate, opt::seedLength, 
                                                         opt::seedStride, true);
    
    std::ostream* pWriter = createWriter(opt::outFile);

    ConnectPostProcess postProcessor(pWriter);

    if(opt::numThreads <= 1)
    {
        // Serial mode
        ConnectProcess processor(pOverlapper, opt::minOverlap);
        SequenceProcessFramework::processSequencesSerial<SequenceWorkItemPair,
                                                         ConnectResult, 
                                                         ConnectProcess, 
                                                         ConnectPostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        // Parallel mode
        std::vector<ConnectProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            ConnectProcess* pProcessor = new ConnectProcess(pOverlapper, opt::minOverlap);
            processorVector.push_back(pProcessor);
        }
        
        SequenceProcessFramework::processSequencesParallel<SequenceWorkItemPair,
                                                           ConnectResult, 
                                                           ConnectProcess, 
                                                           ConnectPostProcess>(opt::readsFile, processorVector, &postProcessor); 
        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
    }

    delete pBWT;
    delete pRBWT;
    delete pOverlapper;
    delete pWriter;
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseConnectOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'm': arg >> opt::minOverlap; break;
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 'e': arg >> opt::errorRate; break;
            case 't': arg >> opt::numThreads; break;
            case 'l': arg >> opt::seedLength; break;
            case 's': arg >> opt::seedStride; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << CONNECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CONNECT_VERSION_MESSAGE;
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
        std::cout << "\n" << CONNECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Validate parameters
    if(opt::errorRate <= 0)
        opt::errorRate = 0.0f;

    if(opt::errorRate > 1.0f)
    {
        std::cerr << "Invalid error-rate parameter: " << opt::errorRate << "\n";
        exit(EXIT_FAILURE);
    }

    if(opt::seedLength < 0)
        opt::seedLength = 0;

    if(opt::seedLength > 0 && opt::seedStride <= 0)
        opt::seedStride = opt::seedLength;

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }

    if(opt::outFile.empty())
    {
        opt::outFile = opt::prefix + ".connect.fa";
    }
}
