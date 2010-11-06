//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct - Correct sequencing errors in reads using the FM-index
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "correct.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "ErrorCorrectProcess.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "correct"
static const char *CORRECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *CORRECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Correct sequencing errors in all the reads in READSFILE\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"          --discard                    detect and discard low-quality reads\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"      -a, --algorithm=STR              specify the correction algorithm to use. STR must be one of hybrid, kmer, overlap.\n"
"                                       The default algorithm is hybrid which first attempts kmer correction, then performs\n"
"                                       overlap correction on the remaining uncorrected reads.\n"
"          --metrics=FILE               collect error correction metrics (error rate by position in read, etc) and write\n"
"                                       them to FILE\n"
"\nOverlap correction parameters:\n"
"      -e, --error-rate                 the maximum error rate allowed between two sequences to consider them overlapped (default: 0.04)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 45)\n"
"      -c, --conflict=INT               use INT as the threshold to detect a conflicted base in the multi-overlap (default: 5)\n"
"      -l, --seed-length=LEN            force the seed length to be LEN. By default, the seed length in the overlap step\n"
"                                       is calculated to guarantee all overlaps with --error-rate differences are found.\n"
"                                       This option removes the guarantee but will be (much) faster. As SGA can tolerate some\n"
"                                       missing edges, this option may be preferable for some data sets.\n"
"      -s, --seed-stride=LEN            force the seed stride to be LEN. This parameter will be ignored unless --seed-length\n"
"                                       is specified (see above). This parameter defaults to the same value as --seed-length\n"
"      -b, --branch-cutoff=N            stop the overlap search at N branches. This parameter is used to control the search time for\n"
"                                       highly-repetitive reads. If the number of branches exceeds N, the search stops and the read\n"
"                                       will not be corrected. This is not enabled by default.\n"
"      -r, --rounds=NUM                 iteratively correct reads up to a maximum of NUM rounds (default: 1)\n"
"\nKmer correction parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 41)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int numRounds = 1;
    static std::string prefix;
    static std::string readsFile;
    static std::string outFile;
    static std::string discardFile;
    static std::string metricsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE;
    
    static double errorRate = 0.04;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
    static int seedLength = 0;
    static int seedStride = 0;
    static int conflictCutoff = 5;
    static int branchCutoff = -1;

    static int kmerLength = 41;
    static int kmerThreshold = 3;

    static ErrorCorrectAlgorithm algorithm = ECA_HYBRID;
}

static const char* shortopts = "p:m:d:e:t:l:s:o:r:b:a:c:k:x:vi";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_DISCARD };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "min-overlap",   required_argument, NULL, 'm' },
    { "rounds",        required_argument, NULL, 'r' },
    { "outfile",       required_argument, NULL, 'o' },
    { "prefix",        required_argument, NULL, 'p' },
    { "error-rate",    required_argument, NULL, 'e' },
    { "seed-length",   required_argument, NULL, 'l' },
    { "seed-stride",   required_argument, NULL, 's' },
    { "algorithm",     required_argument, NULL, 'a' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "conflict",      required_argument, NULL, 'c' },
    { "branch-cutoff", required_argument, NULL, 'b' },
    { "kmer-size",     required_argument, NULL, 'k' },
    { "kmer-threshold",required_argument, NULL, 'x' },
    { "discard",       no_argument,       NULL, OPT_DISCARD },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "metrics",       required_argument, NULL, OPT_METRICS },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int correctMain(int argc, char** argv)
{
    parseCorrectOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);


    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, 
                                                         opt::errorRate, opt::seedLength, 
                                                         opt::seedStride, false, opt::branchCutoff);
    
    std::ostream* pWriter = createWriter(opt::outFile);
    std::ostream* pDiscardWriter = (!opt::discardFile.empty() ? createWriter(opt::discardFile) : NULL);

    bool bCollectMetrics = !opt::metricsFile.empty();
    pBWT->printInfo();

    ErrorCorrectPostProcess postProcessor(pWriter, pDiscardWriter, bCollectMetrics);

    if(opt::numThreads <= 1)
    {
        // Serial mode
        ErrorCorrectProcess processor(pOverlapper, opt::minOverlap, opt::numRounds, opt::conflictCutoff, 
                                      opt::kmerLength, opt::kmerThreshold, opt::algorithm, opt::verbose > 1);

        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                         ErrorCorrectResult, 
                                                         ErrorCorrectProcess, 
                                                         ErrorCorrectPostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        // Parallel mode
        std::vector<ErrorCorrectProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            ErrorCorrectProcess* pProcessor = new ErrorCorrectProcess(pOverlapper, opt::minOverlap, 
                                                                      opt::numRounds, opt::conflictCutoff, 
                                                                      opt::kmerLength, opt::kmerThreshold, 
                                                                      opt::algorithm, opt::verbose > 1);
            processorVector.push_back(pProcessor);
        }
        
        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                           ErrorCorrectResult, 
                                                           ErrorCorrectProcess, 
                                                           ErrorCorrectPostProcess>(opt::readsFile, processorVector, &postProcessor);

        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
    }

    if(bCollectMetrics)
    {
        std::ostream* pMetricsWriter = createWriter(opt::metricsFile);
        postProcessor.writeMetrics(pMetricsWriter);
        delete pMetricsWriter;
    }

    delete pBWT;
    delete pRBWT;
    delete pOverlapper;
    delete pTimer;
    
    delete pWriter;
    if(pDiscardWriter != NULL)
        delete pDiscardWriter;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseCorrectOptions(int argc, char** argv)
{
    std::string algo_str;
    bool bDiscardReads = false;
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
            case 'r': arg >> opt::numRounds; break;
            case 'a': arg >> algo_str; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'c': arg >> opt::conflictCutoff; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'b': arg >> opt::branchCutoff; break;
            case OPT_DISCARD: bDiscardReads = true; break;
            case OPT_METRICS: arg >> opt::metricsFile; break;
            case OPT_HELP:
                std::cout << CORRECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CORRECT_VERSION_MESSAGE;
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

    if(opt::numRounds <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of rounds: " << opt::numRounds << ", must be at least 1\n";
        die = true;
    }
    
    if(opt::kmerLength <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
        die = true;
    }

    if(opt::kmerThreshold <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::kmerThreshold << ", must be greater than zero\n";
        die = true;
    }

    // Determine the correction algorithm to use
    if(!algo_str.empty())
    {
        if(algo_str == "hybrid")
            opt::algorithm = ECA_HYBRID;
        else if(algo_str == "kmer")
            opt::algorithm = ECA_KMER;
        else if(algo_str == "overlap")
            opt::algorithm = ECA_OVERLAP;
        else
        {
            std::cerr << SUBPROGRAM << ": unrecognized -a,--algorithm parameter: " << algo_str << "\n";
            die = true;
        }
    }

    if (die) 
    {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
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

    std::string out_prefix = stripFilename(opt::readsFile);
    if(opt::outFile.empty())
    {
        opt::outFile = out_prefix + ".ec.fa";
    }

    if(bDiscardReads)
    {
        opt::discardFile = out_prefix + ".discard.fa";
    }
    else
    {
        opt::discardFile.clear();
    }
}
