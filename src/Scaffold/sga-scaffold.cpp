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
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string contigsFile;
    static std::string distanceEstFile;
}

static const char* shortopts = "v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
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
    graph.loadVertices(opt::contigsFile);
    graph.visit(statsVisitor);
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
        die = true;
    }

    if(opt::distanceEstFile.empty())
    {
        std::cerr << PROGRAM ": a distance estimation file must be provided\n";
        die = true;
    }

    
}

