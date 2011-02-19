//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// walk - extract a walks between vertices from the graph
//
#include <iostream>
#include <fstream>
#include "walk.h"
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGPairedAlgorithms.h"
#include "SGDebugAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"
#include "SGSearch.h"

// functions

//
// Getopt
//
#define SUBPROGRAM "walk"
static const char *WALK_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *WALK_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE\n"
"Extract walks from ID1 to ID2 from the graph.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -d, --distance=NUM               the maximum length of the walks to find.\n"
"      -s,--startt ID                   start the walk at vertex with ID\n"
"      -e,--end ID                      end the walk at vertex with ID\n"
"      -w,walk-str STR                  extract the defined walk given by STR\n"
"                                       STR should be a comma delimited list of IDs\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;
    static std::string walkStr;
    static std::string id1;
    static std::string id2;
    static std::string outFile;
    static int maxDistance = 500;
}

static const char* shortopts = "o:d:s:e:w:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "out",            required_argument, NULL, 'o' },
    { "distance",       required_argument, NULL, 'd' },
    { "start",          required_argument, NULL, 's' },
    { "end",            required_argument, NULL, 'e' },
    { "walk-str",       required_argument, NULL, 'w' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int walkMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga walk");
    parseWalkOptions(argc, argv);
    walk();
    delete pTimer;

    return 0;
}

void walk()
{
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
    pGraph->printMemSize();
    std::ostream* pWriter = createWriter(opt::outFile);


    SGWalkVector walkVector;

    if(!opt::walkStr.empty())
    {
        // Build the explicitly defined walk.
        std::cout << "Building walk from string: " << opt::walkStr << "\n";
        StringVector vertIDs = split(opt::walkStr, ',');
        assert(vertIDs.size() > 0);
        Vertex* pCurr = pGraph->getVertex(vertIDs.front());
        SGWalk out(pCurr);
        
        for(size_t i = 1; i < vertIDs.size(); ++i)
        {
            Vertex* pNext = pGraph->getVertex(vertIDs[i]);
            EdgePtrVec epv = pCurr->findEdgesTo(pNext->getID());
            if(epv.size() == 0)
            {
                std::cerr << "[sga walk] Error: edge between " << pCurr->getID() << " and " << pNext->getID() << " not found\n";
                exit(EXIT_FAILURE);
            }
            out.addEdge(epv.front());
            pCurr = pNext;
        }

        walkVector.push_back(out);
    }
    else
    {
        // Search the the walk between id1 and id2
        Vertex* pX = pGraph->getVertex(opt::id1);
        Vertex* pY = pGraph->getVertex(opt::id2);
        SGSearch::findWalks(pX, pY, ED_SENSE, opt::maxDistance, 1000, false, walkVector);
        SGSearch::findWalks(pX, pY, ED_ANTISENSE, opt::maxDistance, 1000, false, walkVector);
    }

    for(size_t i = 0; i < walkVector.size(); ++i)
    {
        SGWalk& walk = walkVector[i];
        std::cout << "Walk " << i << "\n";
        walk.print();
        std::string str = walk.getString(SGWT_START_TO_END);
        std::cout << "Str: " << str << "\n";
        std::cout << "Len: " << str.length() << "\n";

        *pWriter << ">test" << i << " " << walk.pathSignature() << "\n";
        *pWriter << str << "\n";
    }

    delete pGraph;
    delete pWriter;
}

// 
// Handle command line arguments
//
void parseWalkOptions(int argc, char** argv)
{
    // Set defaults
    opt::outFile = "walks.fa";

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'd': arg >> opt::maxDistance; break;
            case 's': arg >> opt::id1; break;
            case 'e': arg >> opt::id2; break;
            case 'w': arg >> opt::walkStr; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << WALK_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << WALK_VERSION_MESSAGE;
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

    if(opt::walkStr.empty() && (opt::id1.empty() || opt::id2.empty()))
    {
        std::cerr << SUBPROGRAM ": a walk definition string or start/end vertices must be provided\n";
        exit(EXIT_FAILURE);
    }

    // Parse the input filename
    opt::asqgFile = argv[optind++];
}
