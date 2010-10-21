//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// qc - Perform a quality check on a set of reads, discarding low quality
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "qc.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "QCProcess.h"
#include "BWTDiskConstruction.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "qc"
static const char *QC_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *QC_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Perform a quality check on a set of reads, discarding low quality reads.\n"
"By default, the quality check looks for a tiled set of high-coverage k-mers across the reads.\n"
"This check is fast and can detect chimeric reads or reads with internal uncorrected mismatches.\n"
"Automatically rebuilds the FM-index without the discarded reads.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -o, --outfile=FILE               write the qc-passed reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
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
    static std::string discardFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE;
    
    static int kmerLength = 27;
    static int kmerThreshold = 3;
    static int gapArrayStorage = 4;
}

static const char* shortopts = "p:d:t:o:k:x:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_DISCARD };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "outfile",       required_argument, NULL, 'o' },
    { "prefix",        required_argument, NULL, 'p' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "kmer-size",     required_argument, NULL, 'k' },
    { "kmer-threshold",required_argument, NULL, 'x' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "metrics",       required_argument, NULL, OPT_METRICS },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int qcMain(int argc, char** argv)
{
    parseQCOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);


    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
    pBWT->printInfo();
    
    std::ostream* pWriter = createWriter(opt::outFile);
    std::ostream* pDiscardWriter = createWriter(opt::discardFile);

    QCPostProcess postProcessor(pWriter, pDiscardWriter);
    if(opt::numThreads <= 1)
    {
        // Serial mode
        QCProcess processor(pBWT, pRBWT, opt::kmerLength, opt::kmerThreshold);

        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                         QCResult, 
                                                         QCProcess, 
                                                         QCPostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        // Parallel mode
        std::vector<QCProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            QCProcess* pProcessor = new QCProcess(pBWT, pRBWT, opt::kmerLength, opt::kmerThreshold);
            processorVector.push_back(pProcessor);
        }
        
        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                           QCResult, 
                                                           QCProcess, 
                                                           QCPostProcess>(opt::readsFile, processorVector, &postProcessor);

        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
    }

    // close filehandles
    delete pWriter;
    delete pDiscardWriter;

    // Rebuild the FM-index without the discarded reads
    std::string out_prefix = stripFilename(opt::outFile);
    removeReadsFromIndices(opt::prefix, opt::discardFile, out_prefix, BWT_EXT, SAI_EXT, false, opt::numThreads, opt::gapArrayStorage);
    removeReadsFromIndices(opt::prefix, opt::discardFile, out_prefix, RBWT_EXT, RSAI_EXT, true, opt::numThreads, opt::gapArrayStorage);

    delete pBWT;
    delete pRBWT;
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseQCOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << QC_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << QC_VERSION_MESSAGE;
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

    if (die) 
    {
        std::cout << "\n" << QC_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }

    if(opt::outFile.empty())
    {
        opt::outFile = opt::prefix + ".qcpass.fa";
    }

    opt::discardFile = opt::prefix + ".discard.fa";
}
