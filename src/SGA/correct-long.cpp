//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct-long - Correct sequencing errors in long reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "correct-long.h"
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
#include "CorrectionThresholds.h"
#include "KmerDistribution.h"
#include "LRAlignment.h"

//
// Getopt
//
#define SUBPROGRAM "correct-long"
static const char *CORRECT_LONG_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *CORRECT_LONG_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Correct sequencing errors in the long reads in READSFILE\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
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
    static std::string metricsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
}

static const char* shortopts = "p:m:d:t:o:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_DISCARD, OPT_LEARN };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "min-overlap",   required_argument, NULL, 'm' },
    { "outfile",       required_argument, NULL, 'o' },
    { "prefix",        required_argument, NULL, 'p' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int correctLongMain(int argc, char** argv)
{
    parseCorrectLongOptions(argc, argv);

    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    //BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
    SampledSuffixArray* pSSA = new SampledSuffixArray(opt::prefix + SSA_EXT);

    SeqRecord record;
    SeqReader reader(opt::readsFile);
    while(reader.get(record))
    {
        std::cout << "Aligning sequence " << record.id << "\n";
        LRAlignment::bwaswAlignment(record.seq.toString(), pBWT, pSSA);
    }

    // example 454 read from E. coli K12
    //std::string query = "GGCGTCTTTTATAAAGATGAGCCCATCAAAGAACTGGAGTCGGCGCTGGTGGCGCAAGGCTTTCAGATTATCTGGCCACAAAACAGCGTTGATTTGCTGAAATTTATCGAGCATAACCCTCGAATTTGCGGCGTGATTTTTGACTGGGATGAGTACAGTCTCGATTTATGTAGCGATATCAATCAGCTTAATGAATATCTCCCGCTTTATGCCTTCATCAACACCCACTCGA";
    //std::string query = "TCAAAGAACTGGAGTCGGCGCTGGTGGCGCAAGGCTTTCAGATTAT";
    //LRAlignment::bwaswAlignment(query, pBWT, pSSA);

    delete pBWT;
    //delete pRBWT;
    delete pSSA;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseCorrectLongOptions(int argc, char** argv)
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
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << CORRECT_LONG_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CORRECT_LONG_VERSION_MESSAGE;
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

    if (die) 
    {
        std::cout << "\n" << CORRECT_LONG_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
        opt::prefix = stripFilename(opt::readsFile);

    std::string out_prefix = stripFilename(opt::readsFile);
    if(opt::outFile.empty())
        opt::outFile = out_prefix + ".ec.fa";
}
