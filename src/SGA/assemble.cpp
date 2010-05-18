//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// assemble - convert read overlaps into contigs
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "assemble.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGPairedAlgorithms.h"
#include "SGDebugAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"
#include "SeqDAVG.h"
#include "SuffixTree.h"

//
// Getopt
//
#define SUBPROGRAM "assemble"
static const char *ASSEMBLE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *ASSEMBLE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE\n"
"Create contigs from the assembly graph ASQGFILE.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the contigs to FILE (default: contigs.fa)\n"
"      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
"                                       the overlap set so that the overlap step only needs to be run once.\n"
"      -b, --bubble                     perform bubble removal\n"
"      -t, --trim                       trim terminal branches\n"
"      -c, --correct                    error correct reads and write to correctedReads.fa\n"
"      --edge-stats                     print out the distribution of overlap lengths and number of errors\n"
"                                       for edges found in the overlap step.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;
    static std::string prefix;
    static std::string outFile;
    static std::string debugFile;
    static unsigned int minOverlap;
    static bool bEdgeStats;
    static bool bCorrectReads;
    static bool bRemodelGraph;
    static bool bTrim;
    static bool bBubble;
    static bool bValidate;
}

static const char* shortopts = "p:o:m:d:vbtc";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VALIDATE };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "prefix",         required_argument, NULL, 'p' },
    { "out",            required_argument, NULL, 'o' },
    { "min-overlap",    required_argument, NULL, 'm' },
    { "debug-file",     required_argument, NULL, 'd' },
    { "bubble",         no_argument,       NULL, 'b' },
    { "trim",           no_argument,       NULL, 't' },
    { "correct",        no_argument,       NULL, 'c' },    
    { "remodel",        no_argument,       NULL, 'r' },
    { "edge-stats",     no_argument,       NULL, 'x' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { "validate",       no_argument,       NULL, OPT_VALIDATE},
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int assembleMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga assemble");
    parseAssembleOptions(argc, argv);
    assemble();
    delete pTimer;

    return 0;
}

void assemble()
{
    Timer t("sga assemble");
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, opt::minOverlap, true);
    pGraph->printMemSize();

    // Visitor functors
    SGTransitiveReductionVisitor trVisit;
    SGGraphStatsVisitor statsVisit;
    SGRemodelVisitor remodelVisit;
    SGEdgeStatsVisitor edgeStatsVisit;
    SGTrimVisitor trimVisit;
    SGBubbleVisitor bubbleVisit;
    SGContainRemoveVisitor containVisit;
    SGErrorCorrectVisitor errorCorrectVisit;
    SGValidateStructureVisitor validationVisit;
    SGPairedPathResolveVisitor peResolveVisit;

    if(!opt::debugFile.empty())
    {
        // Pre-assembly graph stats
        std::cout << "Initial graph stats\n";
        pGraph->visit(statsVisit);

        SGDebugGraphCompareVisitor* pDebugGraphVisit = new SGDebugGraphCompareVisitor(opt::debugFile);
        
        /*
        pGraph->visit(*pDebugGraphVisit);
        while(pGraph->visit(realignVisitor))
            pGraph->visit(*pDebugGraphVisit);
        SGOverlapWriterVisitor overlapWriter("final-overlaps.ovr");
        pGraph->visit(overlapWriter);
        */
        //pDebugGraphVisit->m_showMissing = true;
        pGraph->visit(*pDebugGraphVisit);
        pGraph->visit(statsVisit);
        delete pDebugGraphVisit;
        //return;
    }

    if(opt::bEdgeStats)
    {
        std::cout << "Computing edge stats\n";
        pGraph->visit(edgeStatsVisit);
    }

    // Pre-assembly graph stats
    std::cout << "Initial graph stats\n";
    pGraph->visit(statsVisit);    

    // Remove containments from the graph
    std::cout << "Removing contained vertices\n";
    pGraph->visit(containVisit);
    pGraph->writeASQG("afterCR.asqg.gz");
    // Pre-assembly graph stats
    std::cout << "Post-contain graph stats\n";
    pGraph->visit(statsVisit);    

    // Remove transitive edges from the graph
    std::cout << "Removing transitive edges\n";
    pGraph->visit(trVisit);

    // Resolve PE paths
    //pGraph->visit(peResolveVisit);

    if(opt::bValidate)
    {
        std::cout << "Validating graph structure\n";
        pGraph->visit(validationVisit);
    }

    std::cout << "Writing graph file\n";
    pGraph->writeASQG("afterTR.asqg.gz");

    std::cout << "Pre-remodelling graph stats\n";
    pGraph->visit(statsVisit);

    if(opt::bCorrectReads)
    {
        std::cout << "Correcting reads\n";
        pGraph->visit(errorCorrectVisit);

        std::cout << "Writing corrected reads\n";
        SGFastaVisitor correctedVisitor("correctedReads.fa");
        pGraph->visit(correctedVisitor);
        pGraph->writeASQG("afterEC.asqg.gz");
    }

    if(opt::bRemodelGraph)
    {
        // Remodel graph
        std::cout << "Remodelling graph\n";
        pGraph->visit(remodelVisit);
        pGraph->writeASQG("afterRM.asqg.gz");

        while(pGraph->hasContainment())
        {
            std::cout << "Removing contained reads\n";
            pGraph->visit(containVisit);
        }
        pGraph->visit(trVisit);
        std::cout << "After remodel graph stats: \n";
        pGraph->visit(statsVisit);
        //pGraph->writeASQG("afterRM.asqg.gz");
    }

    if(opt::bTrim)
    {
        WARN_ONCE("USING NAIVE TRIMMING");
        std::cout << "Trimming bad vertices\n"; 
        pGraph->visit(trimVisit);
        pGraph->visit(trimVisit);
    }

    // Simplify the graph by compacting edges
    std::cout << "Pre-simplify graph stats\n";
    pGraph->visit(statsVisit);
    pGraph->simplify();

    if(opt::bBubble)
    {
        std::cout << "\nPerforming bubble removal\n";
        // Bubble removal
        pGraph->visit(bubbleVisit);
        pGraph->simplify();
    }    

    std::cout << "\nFinal graph stats\n";
    pGraph->visit(statsVisit);

#ifdef VALIDATE
    VALIDATION_WARNING("SGA/assemble")
    pGraph->validate();
#endif

    // Write the results
    pGraph->writeDot("final.dot");
    SGFastaVisitor av(opt::outFile);
    pGraph->visit(av);

    delete pGraph;
}

// 
// Handle command line arguments
//
void parseAssembleOptions(int argc, char** argv)
{
    // Set defaults
    opt::outFile = "contigs.fa";
    opt::minOverlap = 0;

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 'm': arg >> opt::minOverlap; break;
            case 'd': arg >> opt::debugFile; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'b': opt::bBubble = true; break;
            case 't': opt::bTrim = true; break;
            case 'c': opt::bCorrectReads = true; break;
            case 'r': opt::bRemodelGraph = true; break;
            case 'x': opt::bEdgeStats = true; break;
            case OPT_VALIDATE: opt::bValidate = true; break;
            case OPT_HELP:
                std::cout << ASSEMBLE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << ASSEMBLE_VERSION_MESSAGE;
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
    opt::asqgFile = argv[optind++];
}
