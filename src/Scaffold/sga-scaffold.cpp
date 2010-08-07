//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// sga-scaffold - Post-assembly scaffolder
//
#include <getopt.h>
#include "Util.h"
#include "config.h"
#include "ScaffoldGraph.h"
#include "ScaffoldVisitors.h"

#define PROGRAM "sga-scaffold"
#define AUTHOR "Jared Simpson"

#define SUBPROGRAM "overlap"
static const char *SCAFFOLD_VERSION_MESSAGE =
PROGRAM " Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR ".\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SCAFFOLD_USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION] ... CONTIGSFILE DISTANCE-EST\n"
"Construct scaffolds from CONTIGSFILE using distance estimates in DISTANCE-EST\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -m, --minContingLength=N         only use contigs at least N bp in length to build scaffolds.\n"
"      -a, --astatistic-file=FILE            load Myers' A-statistic values from FILE. This is used to\n"
"                                       determine unique and repetitive contigs with the -u/--unique-astat\n"
"                                       and -r/--repeat-astat parameters\n"
"      -u, --unique-astat=FLOAT         Contigs with an a-statitic value about FLOAT will be considered unique\n"
"      -r, --repeat-astat=FLOAT         Contigs with an a-statistic below FLOAT will be considered repetitive\n"
"                                       Contigs with an a-statistic between these thresholds will not be\n"
"                                       classified as unique or repetitive\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string contigsFile;
    static std::string distanceEstFile;
    static std::string astatFile;
    static double uniqueAstatThreshold = 30.0f;
    static double repeatAstatThreshold = -5.0f;
    static int minContigLength = 0;
}

static const char* shortopts = "vm:a:u:r:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "minLength",      required_argument, NULL, 'm' },
    { "astatistic-file",required_argument, NULL, 'a' },
    { "unique-astat",   required_argument, NULL, 'u' },
    { "repeat-astat",   required_argument, NULL, 'r' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
void parseScaffoldOptions(int argc, char** argv);

//
int main(int argc, char** argv)
{
    parseScaffoldOptions(argc, argv);
    std::cout << "Building scaffolds from " << opt::contigsFile << " using " << opt::distanceEstFile << "\n";

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
    graph.visit(statsVisitor);
    graph.writeDot("scaffold.dot");
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
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 3) 
    {
        std::cerr << PROGRAM ": too many arguments\n";
        die = true;
    }

    if (die) 
    {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    // 
    opt::contigsFile = argv[optind++];
    opt::distanceEstFile = argv[optind++];

    if(opt::contigsFile.empty())
    {
        std::cerr << PROGRAM ": a contigs file must be provided\n";
        exit(1);
    }

    if(opt::distanceEstFile.empty())
    {
        std::cerr << PROGRAM ": a distance estimation file must be provided\n";
        exit(1);
    }

    if(opt::uniqueAstatThreshold < opt::repeatAstatThreshold)
    {
        std::cerr << PROGRAM ": the unique a-stat threshold must be greater than the repeat a-stat threshold\n";
        std::cerr << "Found unique value: " << opt::uniqueAstatThreshold << " repeat value: " << opt::repeatAstatThreshold << "\n";
        exit(1);
    }

    
}

