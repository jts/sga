//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// cluster - extract connected components from an asqg file
//
#include <iostream>
#include <fstream>
#include "cluster.h"
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGPairedAlgorithms.h"
#include "SGDebugAlgorithms.h"
#include "SGVisitors.h"
#include "SGSearch.h"
#include "Timer.h"

//
// Getopt
//
#define SUBPROGRAM "cluster"
static const char *CLUSTER_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *CLUSTER_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ASQGFILE\n"
"Extract clusters of reads belonging to the same connected components.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the clusters to FILE (default: clusters.txt)\n"
"      -m, --min-size=N                 only write clusters with at least N reads (default: 2)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static size_t minSize = 2;
    static std::string asqgFile;
    static std::string outFile;
}

static const char* shortopts = "o:m:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "out",            required_argument, NULL, 'o' },
    { "min-size",       required_argument, NULL, 'm' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int clusterMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga cluster");
    parseClusterOptions(argc, argv);
    cluster();
    delete pTimer;

    return 0;
}

void cluster()
{
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
    pGraph->printMemSize();

    std::ostream* pWriter = createWriter(opt::outFile);

    typedef std::vector<VertexPtrVec> ConnectedComponents;
    
    VertexPtrVec allVertices = pGraph->getAllVertices();
    
    ConnectedComponents components;
    SGSearchTree::connectedComponents(allVertices, components);
    
    size_t numClusters = 0;
    for(size_t i = 0; i < components.size(); ++i)
    {
        VertexPtrVec& currComponent = components[i];

        if(currComponent.size() < opt::minSize)
            continue;

        std::stringstream clusterName;
        clusterName << "cluster-" << numClusters++;
        //std::cout << clusterName.str() << " has " << currComponent.size() << " vertices\n";
        for(size_t j = 0; j < currComponent.size(); ++j)
            *pWriter << clusterName.str() << "\t" << currComponent.size() << "\t" << currComponent[j]->getID() << "\t" << currComponent[j]->getSeq() << "\n";
    }

    std::cout << "Wrote " << numClusters << " clusters with at least " << opt::minSize << " reads\n";
    delete pWriter;
    delete pGraph;
}

// 
// Handle command line arguments
//
void parseClusterOptions(int argc, char** argv)
{
    // Set defaults
    opt::outFile = "clusters.txt";

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'm': arg >> opt::minSize; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << CLUSTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CLUSTER_VERSION_MESSAGE;
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
