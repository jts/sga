//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// filterBAM - remove erroneous paired end connections
// from a bam file.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "filterBAM.h"
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
#include "SGUtil.h"
#include "gmap.h"
#include "SGSearch.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"

// Structs

//
// Getopt
//
#define SUBPROGRAM "filterBAM"
static const char *FILTERBAM_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *FILTERBAM_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE BAMFILE\n"
"Discard mate-pair alignments from a BAM file that are potentially erroneous\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
//"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -m, --max-distance=LEN           search the graph for a path completing the mate-pair fragment. If the path is less than LEN\n"
"                                       then the pair will be discarded.\n"
"      -o, --out-bam=FILE               write the filtered reads to FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;

    static int maxDistance = 700;
    
    static std::string outFile;
    static std::string asqgFile;
    static std::string bamFile;
}

static const char* shortopts = "d:t:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",         no_argument,       NULL, 'v' },
    { "threads",         required_argument, NULL, 't' },
    { "max-distance",    required_argument, NULL, 'd' },
    { "outfile",         required_argument, NULL, 'o' },
    { "help",            no_argument,       NULL, OPT_HELP },
    { "version",         no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int filterBAMMain(int argc, char** argv)
{
    parseFilterBAMOptions(argc, argv);

    // Read the graph and compute walks
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, false);

    Timer* pTimer = new Timer(PROGRAM_IDENT);    

    // 
    int numPairsTotal = 0;
    int numPairsFiltered = 0;
    int numPairsDiffContig = 0;

    // Open the bam file for reading
    BamTools::BamReader* pBamReader = new BamTools::BamReader;
    pBamReader->Open(opt::bamFile);

    // Open the bam file for writing
    BamTools::BamWriter* pBamWriter = new BamTools::BamWriter;
    pBamWriter->Open(opt::outFile, pBamReader->GetHeaderText(), pBamReader->GetReferenceData());

    BamTools::BamAlignment record1;
    BamTools::BamAlignment record2;
    bool done = false;

    const BamTools::RefVector& referenceVector = pBamReader->GetReferenceData();
    while(!done)
    {
        numPairsTotal += 1;

        // Read a pair from the BAM
        // Read record 1. Skip secondary alignments of the previous pair
        do
        {
            if(!pBamReader->GetNextAlignment(record1))
            {
                done = true;
                break;
            }
        } while(!record1.IsPrimaryAlignment());

        // Stop if no alignment could be parsed from the stream
        if(done)
            break;

        // Read record 2. Skip any 
        do
        {
            if(!pBamReader->GetNextAlignment(record2))
            {
                done = true;
                break;
            }
        } while(!record2.IsPrimaryAlignment());

        // If this read failed, there is a mismatch between the pairing
        if(done)
        {
            std::cout << "Could not read pair for read: " << record1.Name << "\n";
        }
    
        if(!record1.IsMapped() || !record2.IsMapped())
            continue;

        // Ensure the pairing is correct
        assert(record1.Name == record2.Name);
        
        std::string vertexID1 = referenceVector[record1.RefID].RefName;
        std::string vertexID2 = referenceVector[record2.RefID].RefName;

        // Get the vertices for this pair using the mapped IDs
        Vertex* pX = pGraph->getVertex(vertexID1);
        Vertex* pY = pGraph->getVertex(vertexID2);

        // Ensure that the vertices are found
        assert(pX != NULL && pY != NULL);

#ifdef DEBUG_CONNECT
        std::cout << "Finding path from " << vertexID1 << " to " << vertexID2 << "\n";
#endif

        EdgeDir walkDirectionXOut = ED_SENSE;
        EdgeDir walkDirectionYIn = ED_SENSE;

        // Flip walk directions if the alignment is to the reverse strand
        if(record1.IsReverseStrand())
            walkDirectionXOut = !walkDirectionXOut;
        
        if(record2.IsReverseStrand())
            walkDirectionYIn = !walkDirectionYIn;

        int fromX = walkDirectionXOut == ED_SENSE ? record1.Position : record1.GetEndPosition();
        int toY = walkDirectionYIn == ED_SENSE ? record2.Position : record2.GetEndPosition();

        // Calculate the amount of contig X that already covers the fragment
        // Using this number, we calculate how far we should search
        int coveredX = walkDirectionXOut == ED_SENSE ? pX->getSeqLen() - fromX : fromX;
        int maxWalkDistance = opt::maxDistance - coveredX;

        bool bShortInsertPair = false;
        if(pX == pY)
        {
            if(abs(record1.InsertSize) < opt::maxDistance)
                bShortInsertPair = true;
        }
        else
        {

            SGWalkVector walks;
            SGSearch::findWalks(pX, pY, walkDirectionXOut, maxWalkDistance, 10000, walks);

            if(!walks.empty())
            {
                for(size_t i = 0; i < walks.size(); ++i)
                {
                    std::string fragment = walks[i].getFragmentString(pX, 
                                                                      pY, 
                                                                      fromX,
                                                                      toY,
                                                                      walkDirectionXOut,
                                                                      walkDirectionYIn);
                    if((int)fragment.size() < opt::maxDistance)
                    {
                        bShortInsertPair = true;
                        numPairsDiffContig += 1;
                        //std::cout << "Found completing fragment (" << pX->getID() << " -> " << pY->getID() << ": " << fragment.size() << "\n";
                        break;
                    }
                }
            }
        }
        
        if(bShortInsertPair)
        {
            numPairsFiltered += 1;
        }
        else
        {
            pBamWriter->SaveAlignment(record1);
            pBamWriter->SaveAlignment(record2);
        }
        if(numPairsTotal++ % 50000 == 0)
            printf("[sga filterBAM] Processed %d pairs\n", numPairsTotal);
    }

    std::cout << "Total pairs: " << numPairsTotal << "\n";
    std::cout << "Total filtered out: " << numPairsFiltered << "\n";
    std::cout << "Total filtered with different contigs: " << numPairsDiffContig << "\n";

    delete pTimer;
    delete pGraph;
    delete pBamReader;
    delete pBamWriter;
    return 0;
}

// 
// Handle command line arguments
//
void parseFilterBAMOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'd': arg >> opt::maxDistance; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << FILTERBAM_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FILTERBAM_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 2) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 2) 
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
        std::cout << "\n" << FILTERBAM_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::asqgFile = argv[optind++];
    opt::bamFile = argv[optind++];

    if(opt::outFile.empty())
    {
        opt::outFile = "filtered." + opt::bamFile;
    }
}
