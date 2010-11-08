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
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE GMAPFILE\n"
"Resolve the complete sequence of a paired end fragment by finding a walk through the graph connecting the ends\n"
"The read adjacency information is given in ASQGFILE, which is the direct output from the sga-overlap step.\n"
"It should not contain duplicate reads. The GMAPFILE specifies the vertices to walk between, read pairs\n"
"are assumed to be on consecutive lines.\n"
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
    static int maxDistance = 250;
    
    static bool bWriteCovered = true;
    static bool bWriteUnresolved = false;

    // In hetSVMode, we will output up to two paths between the pairs
    static bool hetSVMode = false;

    static std::string outFile;
    static std::string unconnectedFile = "unconnected.fa";
    static std::string asqgFile;
    static std::string gmapFile;
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
    std::istream* pReader = createReader(opt::gmapFile);
    std::ostream* pWriter = createWriter(opt::outFile);

    GmapRecord record1;
    GmapRecord record2;
    
    int numPairsAttempted = 0;
    int numPairsResolved = 0;
    int numUnresolvedWrote = 0;
    int numCoveredWrote = 0;


    // In heterozygous SV mode, write up to 2 paths
    size_t maxPaths = (opt::hetSVMode ? 2 : 1);

    while(*pReader >> record1)
    {
        bool good = *pReader >> record2;
        assert(good);
    
        if(!record1.isMapped() || !record2.isMapped())
            continue;

        // Ensure the pairing is correct
        assert(getPairID(record1.readID) == record2.readID);

        // Get the vertices for this pair using the mapped IDs
        Vertex* pX = pGraph->getVertex(record1.mappedID);
        Vertex* pY = pGraph->getVertex(record2.mappedID);

        // Skip the pair if either vertex is not found
        if(pX == NULL || pY == NULL)
            continue;

        EdgeDir walkDirection = ED_SENSE;
        if(record1.isRC)
            walkDirection = !walkDirection;

        SGWalkVector walks;
        SGSearch::findWalks(pX, pY, walkDirection, opt::maxDistance, 10000, walks);

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
                idSS << getPairBasename(record1.readID);
                if(opt::hetSVMode)
                {
                    idSS << "-walk:" << i;
                }

                SeqRecord resolved;
                resolved.id = idSS.str();
                resolved.seq = walks[i].getString(SGWT_START_TO_END);
                resolved.write(*pWriter);
            
                // Mark all the vertices in this walk as resolved
                markWalkVertices(walks[i], GC_BLACK);
            }
            numPairsResolved += 1;
        }
        else if(opt::bWriteUnresolved)
        {
            // Write the unconnected reads
            SeqRecord unresolved1;
            unresolved1.id = record1.readID;
            unresolved1.seq = record1.readSeq;
            
            SeqRecord unresolved2;
            unresolved2.id = record2.readID;
            unresolved2.seq = record2.readSeq;

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
    printf("Wrote %d reads that were covered by a path but not full resolved\n", numCoveredWrote);

    delete pTimer;
    delete pGraph;
    delete pReader;
    delete pWriter;
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
    opt::gmapFile = argv[optind++];

    if(opt::outFile.empty())
    {
        std::string prefix = stripFilename(opt::gmapFile);
        opt::outFile = prefix + ".connect.fa";
        opt::unconnectedFile = prefix + ".single.fa";
    }
}
