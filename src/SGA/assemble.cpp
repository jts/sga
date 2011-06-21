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
#include "EncodedString.h"

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
"      -o, --out-prefix=NAME            use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
"      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
"                                       the overlap set so that the overlap step only needs to be run once.\n"

"\nBubble/Variation removal parameters:\n"
"      -b, --bubble=N                   perform N bubble removal steps (default: 3)\n"
"      -d, --max-divergence=F           only remove variation if the divergence between sequences is less than F (default: 0.05)\n"
"      -g, --max-gap-divergence=F       only remove variation if the divergence between sequences when only counting indels is less than F (default: 0.01)\n"
"                                       Setting this to 0.0 will suppress removing indel variation\n"
"          --max-indel=D                do not remove variation that is an indel of length greater than D (default: 20)\n"
"\n"
"\nTrimming parameters:\n"
"      -x, --cut-terminal=N             cut off terminal branches in N rounds (default: 10)\n"
"      -l, --min-branch-length=LEN      remove terminal branches only if they are less than LEN bases in length (default: 150)\n"

"\nSmall repeat resolution parameters:\n"
"      -r,--resolve-small=LEN           resolve small repeats using spanning overlaps when the difference between the shortest\n"
"                                       and longest overlap is greater than LEN (default: not performed)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;
    static std::string outContigsFile;
    static std::string outVariantsFile;
    static std::string outGraphFile;

    static unsigned int minOverlap;
    static bool bEdgeStats = false;
    static bool bSmoothGraph = false;
    static int resolveSmallRepeatLen = -1;

    // Trim parameters
    static int numTrimRounds = 10;
    static size_t trimLengthThreshold = 150;
    
    // Bubble parameters
    static int numBubbleRounds = 3;
    static double maxBubbleDivergence = 0.05f;
    static double maxBubbleGapDivergence = 0.01f;
    static int maxIndelLength = 20;

    // 
    static int coverageCutoff = 0;
    static bool bValidate;
    static bool bExact = true;
}

static const char* shortopts = "p:o:m:d:g:b:a:c:r:x:l:sv";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VALIDATE, OPT_EDGESTATS, OPT_EXACT, OPT_MAXINDEL };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "out-prefix",         required_argument, NULL, 'o' },
    { "min-overlap",        required_argument, NULL, 'm' },
    { "bubble",             required_argument, NULL, 'b' },
    { "cut-terminal",       required_argument, NULL, 'x' },
    { "min-branch-length",  required_argument, NULL, 'l' },
    { "resolve-small",      required_argument, NULL, 'r' },
    { "coverage",           required_argument, NULL, 'c' },    
    { "max-divergence",     required_argument, NULL, 'd' },
    { "max-gap-divergence", required_argument, NULL, 'g' },
    { "max-indel",          required_argument, NULL, OPT_MAXINDEL },
    { "smooth",             no_argument,       NULL, 's' },
    { "edge-stats",         no_argument,       NULL, OPT_EDGESTATS },
    { "exact",              no_argument,       NULL, OPT_EXACT },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { "validate",           no_argument,       NULL, OPT_VALIDATE},
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
    if(opt::bExact)
        pGraph->setExactMode(true);
    pGraph->printMemSize();

    // Visitor functors
    SGTransitiveReductionVisitor trVisit;
    SGGraphStatsVisitor statsVisit;
    SGRemodelVisitor remodelVisit;
    SGEdgeStatsVisitor edgeStatsVisit;
    SGTrimVisitor trimVisit(opt::trimLengthThreshold);
    SGBubbleVisitor bubbleVisit;
    SGBubbleEdgeVisitor bubbleEdgeVisit;

    SGContainRemoveVisitor containVisit;
    SGErrorCorrectVisitor errorCorrectVisit;
    SGValidateStructureVisitor validationVisit;

    // Pre-assembly graph stats
    std::cout << "Initial graph stats\n";
    pGraph->visit(statsVisit);    

    // Remove containments from the graph
    std::cout << "Removing contained vertices\n";
    while(pGraph->hasContainment())
        pGraph->visit(containVisit);

    // Pre-assembly graph stats
    std::cout << "Post-contain removal graph stats\n";
    pGraph->visit(statsVisit);    

    // Remove any extraneous transitive edges that may remain in the graph
    std::cout << "Removing transitive edges\n";
    pGraph->visit(trVisit);

    // Compact together unbranched chains of vertices
    pGraph->simplify();
    
    if(opt::bValidate)
    {
        std::cout << "Validating graph structure\n";
        pGraph->visit(validationVisit);
    }

    //
    std::cout << "Pre-remodelling graph stats\n";
    pGraph->visit(statsVisit);

    // Remove dead-end branches from the graph
    if(opt::numTrimRounds > 0)
    {
        std::cout << "Trimming bad vertices\n"; 
        int numTrims = opt::numTrimRounds;
        while(numTrims-- > 0)
           pGraph->visit(trimVisit);
        std::cout << "After trimming stats\n";
        pGraph->visit(statsVisit);
    }

    // Resolve small repeats
    if(opt::resolveSmallRepeatLen > 0)
    {
        SGSmallRepeatResolveVisitor smallRepeatVisit(opt::resolveSmallRepeatLen);
        std::cout << "Resolving small repeats\n";

        while(pGraph->visit(smallRepeatVisit)) {}
        
        std::cout << "After small repeat resolve graph stats\n";
        pGraph->visit(statsVisit);
    }

    //
    if(opt::coverageCutoff > 0)
    {
        std::cout << "Coverage visit\n";
        SGCoverageVisitor coverageVisit(opt::coverageCutoff);
        pGraph->visit(coverageVisit);
        pGraph->visit(trimVisit);
        pGraph->visit(trimVisit);
        pGraph->visit(trimVisit);
    }

    // Peform another round of simplification
    pGraph->simplify();
    
    if(opt::numBubbleRounds > 0)
    {
        std::cout << "\nPerforming variation smoothing\n";
        SGSmoothingVisitor smoothingVisit(opt::outVariantsFile, opt::maxBubbleGapDivergence, opt::maxBubbleDivergence, opt::maxIndelLength);
        int numSmooth = opt::numBubbleRounds;
        while(numSmooth-- > 0)
            pGraph->visit(smoothingVisit);
        pGraph->simplify();
    }
    
    pGraph->renameVertices("contig-");

    std::cout << "\nFinal graph stats\n";
    pGraph->visit(statsVisit);

    // Rename the vertices to have contig IDs instead of read IDs
    //pGraph->renameVertices("contig-");

    // Write the results
    SGFastaVisitor av(opt::outContigsFile);
    pGraph->visit(av);

    pGraph->writeASQG(opt::outGraphFile);

    delete pGraph;
}

// 
// Handle command line arguments
//
void parseAssembleOptions(int argc, char** argv)
{
    // Set defaults
    opt::minOverlap = 0;
    std::string prefix = "default";
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> prefix; break;
            case 'm': arg >> opt::minOverlap; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'l': arg >> opt::trimLengthThreshold; break;
            case 'b': arg >> opt::numBubbleRounds; break;
            case 'd': arg >> opt::maxBubbleDivergence; break;
            case 'g': arg >> opt::maxBubbleGapDivergence; break;
            case 's': opt::bSmoothGraph = true; break;
            case 'x': arg >> opt::numTrimRounds; break;
            case 'c': arg >> opt::coverageCutoff; break;
            case 'r': arg >> opt::resolveSmallRepeatLen; break;
            case OPT_MAXINDEL: arg >> opt::maxIndelLength; break;
            case OPT_EXACT: opt::bExact = true; break;
            case OPT_EDGESTATS: opt::bEdgeStats = true; break;
            case OPT_VALIDATE: opt::bValidate = true; break;
            case OPT_HELP:
                std::cout << ASSEMBLE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << ASSEMBLE_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
                
        }
    }

    // Build the output names
    opt::outContigsFile = prefix + "-contigs.fa";
    opt::outVariantsFile = prefix + "-variants.fa";
    opt::outGraphFile = prefix + "-graph.asqg.gz";

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
