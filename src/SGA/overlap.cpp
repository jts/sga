//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// overlap - compute pairwise overlaps between reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "LCPArray.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "AssembleExact.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "OverlapProcess.h"

//
enum OutputType
{
    OT_ASQG,
    OT_RAW
};

// Functions
size_t computeHitsSerial(const std::string& prefix, const std::string& readsFile, 
                         const OverlapAlgorithm* pOverlapper, int minOverlap, 
                         StringVector& filenameVec, std::ostream* pASQGWriter);

size_t computeHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                           const OverlapAlgorithm* pOverlapper, int minOverlap, 
                           StringVector& filenameVec, std::ostream* pASQGWriter);

//
// Getopt
//
#define SUBPROGRAM "overlap"
static const char *OVERLAP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *OVERLAP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Compute pairwise overlap between all the sequences in READS\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences aligned\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -i, --irreducible                only output the irreducible edges for each node\n"
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
    static OutputType outputType = OT_ASQG;
    static std::string prefix;
    static std::string readsFile;
    
    static double errorRate;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
    static int seedLength = 0;
    static int seedStride = 0;
    static bool bIrreducibleOnly;
}

static const char* shortopts = "p:m:d:e:t:l:s:vi";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "threads",     required_argument, NULL, 't' },
    { "min-overlap", required_argument, NULL, 'm' },
    { "max-diff",    required_argument, NULL, 'd' },
    { "prefix",      required_argument, NULL, 'p' },
    { "error-rate",  required_argument, NULL, 'e' },
    { "seed-length", required_argument, NULL, 'l' },
    { "seed-stride", required_argument, NULL, 's' },
    { "irreducible", no_argument,       NULL, 'i' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int overlapMain(int argc, char** argv)
{
    parseOverlapOptions(argc, argv);

    // Prepare the output ASQG file
    assert(opt::outputType == OT_ASQG);
    (void)opt::outputType;

    // Open output file
    std::string asqgFilename = opt::prefix + ASQG_EXT + GZIP_EXT;
    std::ostream* pASQGWriter = createWriter(asqgFilename);

    // Build and write the ASQG header
    ASQG::HeaderRecord headerRecord;
    headerRecord.setOverlapTag(opt::minOverlap);
    headerRecord.setErrorRateTag(opt::errorRate);
    headerRecord.setInputFileTag(opt::readsFile);
    headerRecord.setContainmentTag(true); // containments are always present
    headerRecord.setTransitiveTag(!opt::bIrreducibleOnly);
    headerRecord.write(*pASQGWriter);

    // Compute the overlap hits
    StringVector hitsFilenames;
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, 
                                                         opt::errorRate, opt::seedLength, 
                                                         opt::seedStride, opt::bIrreducibleOnly);
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    size_t count;
    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode overlap computation\n", PROGRAM_IDENT);
        count = computeHitsSerial(opt::prefix, opt::readsFile, pOverlapper, opt::minOverlap, hitsFilenames, pASQGWriter);
    }
    else
    {
        printf("[%s] starting parallel-mode overlap computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        count = computeHitsParallel(opt::numThreads, opt::prefix, opt::readsFile, pOverlapper, opt::minOverlap, hitsFilenames, pASQGWriter);
    }

    delete pOverlapper;
    delete pBWT; 
    delete pRBWT;

    // Parse the hits files and write the overlaps to the ASQG file
    convertHitsToASQG(hitsFilenames, pASQGWriter);

    // Cleanup
    delete pASQGWriter;
    delete pTimer;
    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// Compute the hits for each read in the input file without threading
// Return the number of reads processed
size_t computeHitsSerial(const std::string& prefix, const std::string& readsFile, 
                         const OverlapAlgorithm* pOverlapper, int minOverlap, 
                         StringVector& filenameVec, std::ostream* pASQGWriter)
{
    std::string filename = prefix + HITS_EXT + GZIP_EXT;
    filenameVec.push_back(filename);

    OverlapProcess processor(filename, pOverlapper, minOverlap);
    OverlapPostProcess postProcessor(pASQGWriter, pOverlapper);

    size_t numProcessed = 
           SequenceProcessFramework::processSequencesSerial<OverlapResult, 
                                                            OverlapProcess, 
                                                            OverlapPostProcess>(readsFile, &processor, &postProcessor);
    return numProcessed;
}

// Compute the hits for each read in the SeqReader file with threading
// The way this works is we create a vector of numThreads OverlapProcess pointers and 
// pass this to the SequenceProcessFragmework which wraps the processes
// in threads and distributes the reads to each thread.
// The number of reads processsed is returned
size_t computeHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                           const OverlapAlgorithm* pOverlapper, int minOverlap, 
                           StringVector& filenameVec, std::ostream* pASQGWriter)
{
    std::string filename = prefix + HITS_EXT + GZIP_EXT;

    std::vector<OverlapProcess*> processorVector;
    for(int i = 0; i < numThreads; ++i)
    {
        std::stringstream ss;
        ss << prefix << "-thread" << i << HITS_EXT << GZIP_EXT;
        std::string outfile = ss.str();
        filenameVec.push_back(outfile);
        OverlapProcess* pProcessor = new OverlapProcess(outfile, pOverlapper, minOverlap);
        processorVector.push_back(pProcessor);
    }

    // The post processing is performed serially so only one post processor is created
    OverlapPostProcess postProcessor(pASQGWriter, pOverlapper);
    
    size_t numProcessed = 
           SequenceProcessFramework::processSequencesParallel<OverlapResult, 
                                                              OverlapProcess, 
                                                              OverlapPostProcess>(readsFile, processorVector, &postProcessor);
    for(int i = 0; i < numThreads; ++i)
        delete processorVector[i];
    return numProcessed;
}

//
void convertHitsToASQG(const StringVector& hitsFilenames, std::ostream* pASQGWriter)
{
    // Load the suffix array index and the reverse suffix array index
    // Note these are not the full suffix arrays
    SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
    SuffixArray* pRevSAI = new SuffixArray(opt::prefix + RSAI_EXT);

    // Load the read table and output the initial vertex set, consisting of all the reads
    ReadTable* pFwdRT = new ReadTable(opt::readsFile);

    // Convert the hits to overlaps and write them to the asqg file as initial edges
    for(StringVector::const_iterator iter = hitsFilenames.begin(); iter != hitsFilenames.end(); ++iter)
    {
        printf("[%s] parsing file %s\n", PROGRAM_IDENT, iter->c_str());
        std::istream* pReader = createReader(*iter);
    
        // Read each hit sequentially, converting it to an overlap
        std::string line;
        while(getline(*pReader, line))
        {
            size_t readIdx;
            bool isSubstring;
            OverlapVector ov;
            OverlapCommon::parseHitsString(line, pFwdRT, pFwdSAI, pRevSAI, readIdx, ov, isSubstring);
            for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter)
            {
                ASQG::EdgeRecord edgeRecord(*iter);
                edgeRecord.write(*pASQGWriter);
            }
        }

        delete pReader;
    }

    // Delete allocated data
    delete pFwdSAI;
    delete pRevSAI;
    delete pFwdRT;
}

// Before sanity checks on the overlaps and write them out
void writeOverlap(Overlap& ovr, std::ofstream& containHandle, std::ofstream& overlapHandle)
{
    // Ensure that the overlap is not a containment
    if(ovr.match.coord[0].isContained() || ovr.match.coord[1].isContained())
    {
        containHandle << ovr << "\n";
        return;
    }

    // Unless both overlaps are extreme, skip
    if(!ovr.match.coord[0].isExtreme() || !ovr.match.coord[1].isExtreme())
    {
        std::cerr << "Skipping non-extreme overlap: " << ovr << "\n";
        return;
    }

    bool sameStrand = !ovr.match.isRC();
    bool proper = false;
    
    if(sameStrand)
    {
        proper = (ovr.match.coord[0].isLeftExtreme() != ovr.match.coord[1].isLeftExtreme() && 
                  ovr.match.coord[0].isRightExtreme() != ovr.match.coord[1].isRightExtreme());
    }
    else
    {
        proper = (ovr.match.coord[0].isLeftExtreme() == ovr.match.coord[1].isLeftExtreme() && 
                  ovr.match.coord[0].isRightExtreme() == ovr.match.coord[1].isRightExtreme());
    }
    
    if(!proper)
    {
        std::cerr << "Skipping improper overlap: " << ovr << "\n";
        return;
    }

    // All checks passed, output the overlap
    overlapHandle << ovr << "\n";
}

// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'm': arg >> opt::minOverlap; break;
            case 'p': arg >> opt::prefix; break;
            case 'e': arg >> opt::errorRate; break;
            case 't': arg >> opt::numThreads; break;
            case 'l': arg >> opt::seedLength; break;
            case 's': arg >> opt::seedStride; break;
            case 'i': opt::bIrreducibleOnly = true; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << OVERLAP_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << OVERLAP_VERSION_MESSAGE;
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
        std::cout << "\n" << OVERLAP_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Validate parameters
    if(opt::errorRate <= 0)
        opt::errorRate = 0.0f;

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
}
