//-----------------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
// Released under the GPL
//-----------------------------------------------------
//
// overlap-long - compute overlaps for long reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "overlap-long.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "OverlapProcess.h"
#include "ReadInfoTable.h"
#include "KmerOverlaps.h"

// Functions
size_t computeHitsSerial(const std::string& prefix, const std::string& readsFile, 
                         const OverlapAlgorithm* pOverlapper, int minOverlap, 
                         StringVector& filenameVec, std::ostream* pASQGWriter);

size_t computeHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                           const OverlapAlgorithm* pOverlapper, int minOverlap, 
                           StringVector& filenameVec, std::ostream* pASQGWriter);

//
void convertHitsToASQG(const std::string& indexPrefix, const StringVector& hitsFilenames, std::ostream* pASQGWriter);


//
// Getopt
//
#define SUBPROGRAM "overlap-long"
static const char *OVERLAP_LONG_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Ontario Institute for Cancer Research\n";

static const char *OVERLAP_LONG_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Compute pairwise overlap between all the long reads in READS\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM worker threads to compute the overlaps (default: no threading)\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences aligned (default: exact matches only)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 45)\n"
"      -f, --target-file=FILE           perform the overlap queries against the reads in FILE\n"
"      -l, --seed-length=LEN            force the seed length to be LEN. By default, the seed length in the overlap step\n"
"                                       is calculated to guarantee all overlaps with --error-rate differences are found.\n"
"                                       This option removes the guarantee but will be (much) faster. As SGA can tolerate some\n"
"                                       missing edges, this option may be preferable for some data sets.\n"
"      -s, --seed-stride=LEN            force the seed stride to be LEN. This parameter will be ignored unless --seed-length\n"
"                                       is specified (see above). This parameter defaults to the same value as --seed-length\n"
"      -d, --sample-rate=N              sample the symbol counts every N symbols in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string readsFile;
    static std::string targetFile;
    static std::string outFile;
    
    static double errorRate = 0.0f;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
    static int seedLength = 100;
    static int seedStride = 0;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    static bool bIrreducibleOnly = true;
    static bool bExactIrreducible = false;
}

static const char* shortopts = "m:d:e:t:l:s:o:f:vix";

enum { OPT_HELP = 1, OPT_VERSION, OPT_EXACT };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "threads",     required_argument, NULL, 't' },
    { "min-overlap", required_argument, NULL, 'm' },
    { "sample-rate", required_argument, NULL, 'd' },
    { "outfile",     required_argument, NULL, 'o' },
    { "target-file", required_argument, NULL, 'f' },
    { "prefix",      required_argument, NULL, 'p' },
    { "error-rate",  required_argument, NULL, 'e' },
    { "seed-length", required_argument, NULL, 'l' },
    { "seed-stride", required_argument, NULL, 's' },
    { "exhaustive",  no_argument,       NULL, 'x' },
    { "exact",       no_argument,       NULL, OPT_EXACT },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Write an ascii pictogram of an overlap
// O = overhang of length block_size (region of no overlap)
// P = perfect match of length block_size
// I = matched block with mismatches
std::string ascii_overlap(const std::string& s0, const std::string& s1,
                          SequenceOverlap& ovr, int block_size = 50)
{
    std::string out;

    int pre_overhang = ovr.match[0].start / block_size;
    int post_overhang = (ovr.length[0] - ovr.match[0].end - 1) / block_size;

    std::string p0;
    std::string p1;

    ovr.makePaddedMatches(s0, s1, &p0, &p1);

    out.append(pre_overhang, '-');

    for(size_t i = 0; i < p0.size(); i += block_size)
    {
        if(p0.substr(i, block_size) == p1.substr(i, block_size))
            out.append(1, '=');
        else
            out.append(1, 'x');
    }
    out.append(post_overhang, '-');
    return out;
}

//
// Main
//
int overlapLongMain(int argc, char** argv)
{
    parseOverlapLongOptions(argc, argv);

    // Open output file
    std::ostream* pASQGWriter = createWriter(opt::outFile);

    // Build and write the ASQG header
    ASQG::HeaderRecord headerRecord;
    headerRecord.setOverlapTag(opt::minOverlap);
    headerRecord.setErrorRateTag(opt::errorRate);
    headerRecord.setInputFileTag(opt::readsFile);
    headerRecord.setTransitiveTag(true);
    headerRecord.write(*pASQGWriter);

    // Determine which index files to use. If a target file was provided,
    // use the index of the target reads
    std::string indexPrefix;
    if(!opt::targetFile.empty())
        indexPrefix = stripFilename(opt::targetFile);
    else
        indexPrefix = stripFilename(opt::readsFile);

    BWT* pBWT = new BWT(indexPrefix + BWT_EXT, opt::sampleRate);
    SampledSuffixArray* pSSA = new SampledSuffixArray(indexPrefix + SAI_EXT, SSA_FT_SAI);
    
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    pBWT->printInfo();

    // Read the sequence file and write vertex records for each
    // Also store the read names in a vector of strings
    ReadTable reads;
    
    SeqReader* pReader = new SeqReader(opt::readsFile, SRF_NO_VALIDATION);
    SeqRecord record;
    while(pReader->get(record))
    {
        reads.addRead(record.toSeqItem());
        ASQG::VertexRecord vr(record.id, record.seq.toString());
        vr.write(*pASQGWriter);

        if(reads.getCount() % 100000 == 0)
            printf("Read %zu sequences\n", reads.getCount());
    }

    delete pReader;
    pReader = NULL;

    BWTIndexSet index;
    index.pBWT = pBWT;
    index.pSSA = pSSA;
    index.pReadTable = &reads;

    // Make a prefix for the temporary hits files
    size_t n_reads = reads.getCount();

#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
    #pragma omp parallel for
#endif
    for(size_t read_idx = 0; read_idx < n_reads; ++read_idx)
    {
        const SeqItem& curr_read = reads.getRead(read_idx);

        printf("read %s %zubp\n", curr_read.id.c_str(), curr_read.seq.length());
        SequenceOverlapPairVector sopv = 
            KmerOverlaps::retrieveMatches(curr_read.seq.toString(),
                                          opt::seedLength,
                                          opt::minOverlap,
                                          1 - opt::errorRate,
                                          100,
                                          index);

        printf("Found %zu matches\n", sopv.size());
        for(size_t i = 0; i < sopv.size(); ++i)
        {
            std::string match_id = reads.getRead(sopv[i].match_idx).id;

            // We only want to output each edge once so skip this overlap
            // if the matched read has a lexicographically lower ID
            if(curr_read.id > match_id)
                continue;

            std::string ao = ascii_overlap(sopv[i].sequence[0], sopv[i].sequence[1], sopv[i].overlap, 50);
            printf("\t%s\t[%d %d] ID=%s OL=%d PI:%.2lf C=%s\n", ao.c_str(),
                                                                sopv[i].overlap.match[0].start,
                                                                sopv[i].overlap.match[0].end,
                                                                match_id.c_str(),
                                                                sopv[i].overlap.getOverlapLength(),
                                                                sopv[i].overlap.getPercentIdentity(),
                                                                sopv[i].overlap.cigar.c_str());

            // Convert to ASQG
            SeqCoord sc1(sopv[i].overlap.match[0].start, sopv[i].overlap.match[0].end, sopv[i].overlap.length[0]);
            SeqCoord sc2(sopv[i].overlap.match[1].start, sopv[i].overlap.match[1].end, sopv[i].overlap.length[1]);
            
            // KmerOverlaps returns the coordinates of the overlap after flipping the reads
            // to ensure the strand matches. The ASQG file wants the coordinate of the original
            // sequencing strand. Flip here if necessary
            if(sopv[i].is_reversed)
                sc2.flip();

            // Convert the SequenceOverlap the ASQG's overlap format
            Overlap ovr(curr_read.id, sc1, match_id,  sc2, sopv[i].is_reversed, -1);

            ASQG::EdgeRecord er(ovr);
            er.setCigarTag(sopv[i].overlap.cigar);
            er.setPercentIdentityTag(sopv[i].overlap.getPercentIdentity());

#pragma omp critical
            {
                er.write(*pASQGWriter);
            }
        }
    }

    // Cleanup
    delete pReader;
    delete pBWT; 
    delete pSSA;
    
    delete pASQGWriter;
    delete pTimer;
    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

/*
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
           SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                            OverlapResult, 
                                                            OverlapProcess, 
                                                            OverlapPostProcess>(readsFile, &processor, &postProcessor);
    return numProcessed;
}

// Compute the hits for each read in the SeqReader file with threading
// The way this works is we create a vector of numThreads OverlapProcess pointers and 
// pass this to the SequenceProcessFramework which wraps the processes
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
           SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                              OverlapResult, 
                                                              OverlapProcess, 
                                                              OverlapPostProcess>(readsFile, processorVector, &postProcessor);
    for(int i = 0; i < numThreads; ++i)
        delete processorVector[i];
    return numProcessed;
}
*/

// 
// Handle command line arguments
//
void parseOverlapLongOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'm': arg >> opt::minOverlap; break;
            case 'o': arg >> opt::outFile; break;
            case 'e': arg >> opt::errorRate; break;
            case 't': arg >> opt::numThreads; break;
            case 'l': arg >> opt::seedLength; break;
            case 's': arg >> opt::seedStride; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'f': arg >> opt::targetFile; break;
            case OPT_EXACT: opt::bExactIrreducible = true; break;
            case 'x': opt::bIrreducibleOnly = false; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << OVERLAP_LONG_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << OVERLAP_LONG_VERSION_MESSAGE;
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

    if(!IS_POWER_OF_2(opt::sampleRate))
    {
        std::cerr << SUBPROGRAM ": invalid parameter to -d/--sample-rate, must be power of 2. got: " << opt::sampleRate << "\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << OVERLAP_LONG_USAGE_MESSAGE;
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

    if(opt::outFile.empty())
    {
        std::string prefix = stripFilename(opt::readsFile);
        if(!opt::targetFile.empty())
        {
            prefix.append(1,'.');
            prefix.append(stripFilename(opt::targetFile));
        }
        opt::outFile = prefix + ASQG_EXT + GZIP_EXT;
    }
}
