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
#include "SGUtil.h"
#include "gmap.h"
#include "SGSearch.h"
#include "api/BamReader.h"

//#define DEBUG_CONNECT 1

// Struct
// Functions
void markWalkVertices(SGWalk& walk, GraphColor color);

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
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE BAMFILE\n"
"Resolve the complete sequence of a paired end fragment by finding a walk through the graph connecting the ends.\n"
"The graph is read from the ASQGFILE. The reads alignments (to the contigs that make up the graph graph) are\n"
"read from the BAMFILE.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
//"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -m, --max-distance=LEN           maximum expected distance between the PE reads (start to end). This option specifies\n"
"                                       how long the search should proceed for. Default: 250\n"
"      -o, --outfile=FILE               write the connected reads to FILE\n"
"          --connected-only             only write the connected read pairs, not the unresolved vertices in the graph\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;

    static int minDistance = 150;
    static int maxDistance = 400;
    
    static bool bWriteCovered = true;
    static bool bWriteUnresolved = false;

    // In hetSVMode, we will output up to two paths between the pairs
    static bool hetSVMode = false;

    static std::string outFile;
    static std::string unconnectedFile = "unconnected.fa";
    static std::string asqgFile;
    static std::string bamFile;
}

static const char* shortopts = "p:m:e:t:l:s:o:d:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_HETSV, OPT_CONNECTED_ONLY };

static const struct option longopts[] = {
    { "verbose",         no_argument,       NULL, 'v' },
    { "threads",         required_argument, NULL, 't' },
    { "max-distance",    required_argument, NULL, 'd' },
    { "outfile",         required_argument, NULL, 'o' },
    { "connected-only",  no_argument,       NULL, OPT_CONNECTED_ONLY },
    { "het-sv",          no_argument,       NULL, OPT_HETSV },
    { "help",            no_argument,       NULL, OPT_HELP },
    { "version",         no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Visit every vertex in the graph and write out
// the vertices that were covered by a path but not
// resolved
struct CoveredVertexVisitor
{

    CoveredVertexVisitor(std::ostream* pWriter) : m_pWriter(pWriter), resolved(0), internal(0), not_used(0) {}
    ~CoveredVertexVisitor()
    {
        std::cout << "Resolved: " << resolved << "\n";
        std::cout << "Internal: " << internal << "\n";
        std::cout << "Not used: " << not_used << "\n";
    }

    // not used
    void previsit(StringGraph*) {}
    void postvisit(StringGraph*) {}

    //
    bool visit(StringGraph*, Vertex* pX)
    {
        if(pX->getColor() == GC_BLACK)
        {
            resolved += 1;
        }

        if(pX->getColor() == GC_WHITE)
        {
            not_used += 1;
        }

        if(pX->getColor() == GC_RED)
        {
            SeqRecord record;
            record.id = pX->getID();
            record.seq = pX->getSeq().toString();
            record.write(*m_pWriter);
            internal += 1;
        }

        return false;
    }

    int getNumWrote() const { return internal; }

    //
    std::ostream* m_pWriter;
    int resolved;
    int internal;
    int not_used;

};

//
// Main
//
int connectMain(int argc, char** argv)
{
    parseConnectOptions(argc, argv);

    // Read the graph and compute walks
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, false);

    Timer* pTimer = new Timer(PROGRAM_IDENT);    

    // Open the bam file for reading
    BamTools::BamReader* pBamReader = new BamTools::BamReader;
    pBamReader->Open(opt::bamFile);

    std::ostream* pWriter = createWriter(opt::outFile);
    std::cout << "NUM REFS: " << pBamReader->GetReferenceCount() << "\n";
    
    int numPairsAttempted = 0;
    int numPairsResolved = 0;
    int numUnresolvedWrote = 0;
    int numCoveredWrote = 0;


    // In heterozygous SV mode, write up to 2 paths
    size_t maxPaths = (opt::hetSVMode ? 2 : 1);

    BamTools::BamAlignment record1;
    BamTools::BamAlignment record2;
    bool done = false;

    const BamTools::RefVector& referenceVector = pBamReader->GetReferenceData();
    while(!done)
    {
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

        // Skip the pair if either vertex is not found
        if(pX == NULL || pY == NULL)
            continue;

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

        // Calculate the coordinates of the sequenced part of the contigs
        // that these reads lie on
        SeqCoord sc1(0,0, pX->getSeqLen());
        sc1.interval.start = record1.Position;
        sc1.interval.end = record1.GetEndPosition();

        SeqCoord sc2(0,0, pY->getSeqLen());
        sc2.interval.start = record2.Position;
        sc2.interval.end = record2.GetEndPosition();

        int skipX = 0;
        if(!record1.IsReverseStrand())
            skipX = record1.Position;
        else
            skipX = pX->getSeqLen() - record1.GetEndPosition() - 1;

        int skipY = 0;
        if(!record2.IsReverseStrand())
            skipY = record2.Position;
        else
            skipY = pY->getSeqLen() - record2.GetEndPosition() - 1;

        SGWalkVector walks;
        SGSearch::findWalks(pX, pY, walkDirectionXOut, opt::maxDistance, 1000, walks);

        // Mark used vertices in the graph
        // If the entire path was resolved, mark black
        // otherwise mark as red
        GraphColor used_color = (walks.size() <= maxPaths) ? GC_BLACK : GC_RED;

        for(size_t i = 0; i < walks.size(); i +=1 )
        {
            markWalkVertices(walks[i], used_color);
        }

        if(walks.size() <= maxPaths)
        {
            for(size_t i = 0; i < walks.size(); ++i)
            {
                std::stringstream idSS;
                idSS << getPairBasename(record1.Name);
                size_t numReads = walks[i].getNumVertices();

                // Calculate the alignment positions of the reads on the merged path string
                // First, if the two contigs are not from the same strand flip the coordinate
                // of the second read
                bool needFlip = !walks[i].areEndpointsFromSameStrand();
                if(needFlip)
                    sc2.flip();

                // Next, translate the position of the right-most coordinate
                // (wrt to the path string) by the added length
                if(pX != pY)
                {
                    // We need to translate one of the two coordinates
                    // If the extension was antisense, s1 is translated
                    // otherwise s2 is.
                    int overlap_distance = walks[i].getEndToStartDistance();
                    if(walkDirectionXOut == ED_ANTISENSE)
                    {
                        int translate = pY->getSeqLen() + overlap_distance;
                        sc1.interval.start += translate;
                        sc1.interval.end += translate;
                    }
                    else
                    {
                        int translate = pX->getSeqLen() + overlap_distance;
                        sc2.interval.start += translate;
                        sc2.interval.end += translate;
                    }
                }


                // Update the length of the coordinates
                int totalLength = walks[i].getStartToEndDistance();
                sc1.seqlen = totalLength;
                sc2.seqlen = totalLength;

                assert(sc1.isValid());
                assert(sc2.isValid());

                // Validate that the path is as expected
                // This has 2 conditions:
                // 1) The inferred fragment is orientated correctly
                // 2) The fragment size is within the expected range
                bool correctOrientation = false;
                if(walkDirectionXOut == ED_SENSE)
                {
                    if(sc1.interval.start <= sc2.interval.start)
                        correctOrientation = true;
                }
                else
                {
                    if(sc1.interval.start >= sc2.interval.start)
                        correctOrientation = true;
                }

#ifdef DEBUG_CONNECT
                std::cout << "FINAL SC1: " << sc1 << "\n";
                std::cout << "FINAL SC2: " << sc2 << "\n";
                std::cout << "OrientationCheck: " << correctOrientation << "\n";
#endif

                // Calculate the seqcoord on the path string representing the paired end fragment
                SeqCoord pe(0,0, totalLength);
                pe.interval.start = std::min(sc1.interval.start, sc2.interval.start);
                pe.interval.end = std::max(sc1.interval.end, sc2.interval.end);
                int fragSize = pe.length();
                bool correctSize = false;
                if(fragSize >= opt::minDistance && fragSize <= opt::maxDistance)
                    correctSize = true;

                int expectedLength = totalLength - (skipX + skipY);

#ifdef DEBUG_CONNECT
                std::cout << "FRAGSIZE: " << fragSize << "\n";
                std::cout << "EXPECTED: " << expectedLength << "\n";
                std::cout << "SIZE CHECK: " << correctSize << "\n";
#endif
                
                if(correctOrientation && correctSize)
                {
                    // All checks passed, output the fragment
                    //std::string pathString = walks[i].getString(SGWT_START_TO_END);
                    //std::string fragment2 = pathString.substr(pe.interval.start, pe.length());

                    int fromX = walkDirectionXOut == ED_SENSE ? record1.Position : record1.GetEndPosition();
                    int toY = walkDirectionYIn == ED_SENSE ? record2.Position : record2.GetEndPosition();
                    
                    //std::cout << "fromX: " << fromX << "\n";
                    ///std::cout << "toY: " << toY << "\n";

                    std::string fragment = walks[i].getFragmentString(pX, 
                                                                      pY, 
                                                                      fromX,
                                                                      toY,
                                                                      walkDirectionXOut,
                                                                      walkDirectionYIn);
                    
                    /*
                    std::cout << "X: " << pX->getID() << " Y: " << pY->getID() << "\n";
                    std::cout << "Fragsize  " << fragment.size() << "\n";
                    std::cout << "Fragptp  " << pointToPoint.size() << "\n";

                    std::cout << "Fragment: " << fragment << "\n";
                    std::cout << "PTP:      " << pointToPoint << "\n";
                    */
                    
                    //assert(fragment == fragment2);
                    if(opt::hetSVMode)
                    {
                        idSS << "-walk:" << i;
                    }

                    idSS << " numVertices:" << numReads;

                    SeqRecord resolved;
                    resolved.id = idSS.str();
                    resolved.seq = fragment;
                    resolved.write(*pWriter);
                
                    // Mark all the vertices in this walk as resolved
                    markWalkVertices(walks[i], GC_BLACK);
                }
            }
            numPairsResolved += 1;
        }
        else if(opt::bWriteUnresolved)
        {
            // Write the unconnected reads
            SeqRecord unresolved1;
            unresolved1.id = record1.Name;
            unresolved1.seq = record1.QueryBases;
            
            SeqRecord unresolved2;
            unresolved2.id = record2.Name;
            unresolved2.seq = record2.QueryBases;

            unresolved1.write(*pWriter);
            unresolved2.write(*pWriter);

            numUnresolvedWrote += 2;
        }
        numPairsAttempted += 1;

        if(numPairsAttempted % 50000 == 0)
            printf("[sga connect] Processed %d pairs\n", numPairsAttempted);
    }

    //
    if(opt::bWriteCovered)
    {
        CoveredVertexVisitor cvv(pWriter);
        pGraph->visit(cvv);
        numCoveredWrote += cvv.getNumWrote();
    }

    double proc_time_secs = pTimer->getElapsedWallTime();
    printf("connect: Resolved %d out of %d pairs (%lf) in %lfs (%lf pairs/s)\n",
            numPairsResolved, numPairsAttempted, 
            (double)numPairsResolved / numPairsAttempted,
            proc_time_secs,
            numPairsAttempted / proc_time_secs);
    
    printf("Wrote %d unconnected pairs\n", numUnresolvedWrote);
    printf("Wrote %d vertices that were covered by a path but not full resolved\n", numCoveredWrote);

    delete pTimer;
    delete pGraph;
    delete pWriter;
    delete pBamReader;
    return 0;
}

// Mark all the vertices in the walk as color 
// unless they are already black (we never change from black->red)
void markWalkVertices(SGWalk& walk, GraphColor color)
{
    VertexPtrVec verts = walk.getVertices();
    for(size_t i = 0; i < verts.size(); ++i)
    {
        if(verts[i]->getColor() != GC_BLACK)
            verts[i]->setColor(color);
    }
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
            case 'o': arg >> opt::outFile; break;
            case 'd': arg >> opt::maxDistance; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_CONNECTED_ONLY: 
                opt::bWriteCovered = false; 
                opt::bWriteUnresolved = false; 
                break;
            case OPT_HETSV: opt::hetSVMode = true; break;
            case OPT_HELP:
                std::cout << CONNECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CONNECT_VERSION_MESSAGE;
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
        std::cout << "\n" << CONNECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::asqgFile = argv[optind++];
    opt::bamFile = argv[optind++];

    if(opt::outFile.empty())
    {
        std::string prefix = stripFilename(opt::bamFile);
        opt::outFile = prefix + ".connect.fa";
        opt::unconnectedFile = prefix + ".single.fa";
    }
}
