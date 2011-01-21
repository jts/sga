//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// scaffold - Scaffold a set of contigs using
// distances estimates and/or reference information.
//
#include <getopt.h>
#include "Util.h"
#include "config.h"
#include "ScaffoldGraph.h"
#include "ScaffoldVisitors.h"

#define SUBPROGRAM "scaffold"
static const char *SCAFFOLD_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SCAFFOLD_USAGE_MESSAGE =
"Usage: " SUBPROGRAM " [OPTION] ... CONTIGSFILE DISTANCE-EST\n"
"Construct scaffolds from CONTIGSFILE using distance estimates in DISTANCE-EST\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -m, --min-length=N               only use contigs at least N bp in length to build scaffolds (default: no minimun).\n"
"      -a, --astatistic-file=FILE       load Myers' A-statistic values from FILE. This is used to\n"
"                                       determine unique and repetitive contigs with the -u/--unique-astat\n"
"                                       and -r/--repeat-astat parameters (required)\n"
"      -u, --unique-astat=FLOAT         Contigs with an a-statitic value about FLOAT will be considered unique (default: 20.0)\n"
"      -r, --repeat-astat=FLOAT         Contigs with an a-statistic below FLOAT will be considered repetitive (default: 5.0)\n"
"                                       Contigs with an a-statistic between these thresholds will not be\n"
"                                       classified as unique or repetitive\n"
"      -o, --outfile=FILE               write the scaffolds to FILE (default: CONTIGSFILE.scaf\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string contigsFile;
    static std::string distanceEstFile;
    static std::string astatFile;
    static std::string outFile;
    static double uniqueAstatThreshold = 20.0f;
    static double repeatAstatThreshold = 5.0f;
    static int minContigLength = 0;
}

static const char* shortopts = "vm:a:u:r:o:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "min-length",      required_argument, NULL, 'm' },
    { "astatistic-file",required_argument, NULL, 'a' },
    { "unique-astat",   required_argument, NULL, 'u' },
    { "repeat-astat",   required_argument, NULL, 'r' },
    { "outfile",        required_argument, NULL, 'o' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
void parseScaffoldOptions(int argc, char** argv);

//
int scaffoldMain(int argc, char** argv)
{
    parseScaffoldOptions(argc, argv);
    std::cout << "Building scaffolds from " << opt::contigsFile << " using " << opt::distanceEstFile << "\n";

    int maxOverlap = 100;

    ScaffoldStatsVisitor statsVisitor;
    ScaffoldGraph graph;
    
    graph.loadVertices(opt::contigsFile, opt::minContigLength);
    graph.loadDistanceEstimateEdges(opt::distanceEstFile);

    if(!opt::astatFile.empty())
    {
        graph.loadAStatistic(opt::astatFile);
        ScaffoldAStatisticVisitor astatVisitor(opt::uniqueAstatThreshold, 
                                               opt::repeatAstatThreshold);
        graph.visit(astatVisitor);
    }


    // Create chains of vertices from the links
    std::cout << "[sga-scaffold] Removing non-unique vertices from scaffold graph\n";
    graph.writeDot("pregraph.dot");
    graph.deleteVertices(SVC_REPEAT);
    graph.writeDot("scaffold.dot");

    ScaffoldLinkValidator linkValidator(maxOverlap, 0.2f);   
    graph.visit(linkValidator);
    exit(1);
    graph.visit(statsVisitor);


    ScaffoldChainVisitor chainVisitor(maxOverlap);

    // Collapse all chains
    while(graph.visit(chainVisitor)) {}
    graph.writeDot("finalChain.dot");
    graph.visit(statsVisitor);

    ScaffoldMultiEdgeRemoveVisitor cutVisitor;
    graph.visit(cutVisitor);

    ScaffoldWriterVisitor writer(opt::outFile);
    graph.visit(writer);

    return 0;
}

//
void parseScaffoldOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'm': arg >> opt::minContigLength; break;
            case 'a': arg >> opt::astatFile; break;
            case 'u': arg >> opt::uniqueAstatThreshold; break;
            case 'r': arg >> opt::repeatAstatThreshold; break;
            case 'o': arg >> opt::outFile; break;
            case OPT_HELP:
                std::cout << SCAFFOLD_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SCAFFOLD_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 2) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 3) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if (die) 
    {
        std::cerr << "Try `" << SUBPROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    // 
    opt::contigsFile = argv[optind++];
    opt::distanceEstFile = argv[optind++];

    if(opt::contigsFile.empty())
    {
        std::cerr << SUBPROGRAM ": a contigs file must be provided\n";
        exit(1);
    }

    if(opt::distanceEstFile.empty())
    {
        std::cerr << SUBPROGRAM ": a distance estimation file must be provided\n";
        exit(1);
    }

    if(opt::astatFile.empty())
    {
        std::cerr << SUBPROGRAM ": an a-statistic file must be provided\n";
        exit(1);
    }

    if(opt::uniqueAstatThreshold < opt::repeatAstatThreshold)
    {
        std::cerr << SUBPROGRAM ": the unique a-stat threshold must be greater than the repeat a-stat threshold\n";
        std::cerr << "Found unique value: " << opt::uniqueAstatThreshold << " repeat value: " << opt::repeatAstatThreshold << "\n";
        exit(1);
    }

    if(opt::outFile.empty())
    {
        opt::outFile = opt::contigsFile + ".scaf";
    }

    
}

