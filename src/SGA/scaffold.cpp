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
#include "ScaffoldSearch.h"

#define SUBPROGRAM "scaffold"
static const char *SCAFFOLD_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SCAFFOLD_USAGE_MESSAGE =
"Usage: " SUBPROGRAM " [OPTION] ... CONTIGSFILE DISTANCE-EST\n"
"Construct scaffolds from CONTIGSFILE using distance estimates. \n"
"The distance estimates are read from the --pe and --matepair parameters\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"          --pe=FILE                    load links derived from paired-end (short insert) libraries from FILE\n"
"          --mate-pair=FILE             load links derived from mate-pair (long insert) libraries from FILE\n"
"      -m, --min-length=N               only use contigs at least N bp in length to build scaffolds (default: no minimun).\n"
"      -g, --asqg-file=FILE             optionally load the sequence graph from FILE\n"
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
    static std::string peDistanceEstFile;
    static std::string mateDistanceEstFile;
    static std::string astatFile;
    static std::string outFile;
    static std::string asqgFile;
    static double uniqueAstatThreshold = 20.0f;
    static double repeatAstatThreshold = 5.0f;
    static int minContigLength = 0;
}

static const char* shortopts = "vm:a:u:r:o:g:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PE, OPT_MATEPAIR };

static const struct option longopts[] = {
    { "verbose",         no_argument,       NULL, 'v' },
    { "min-length",      required_argument, NULL, 'm' },
    { "asgq-file",       required_argument, NULL, 'g' }, 
    { "astatistic-file", required_argument, NULL, 'a' },
    { "unique-astat",    required_argument, NULL, 'u' },
    { "repeat-astat",    required_argument, NULL, 'r' },
    { "outfile",         required_argument, NULL, 'o' },
    { "pe",              required_argument, NULL, OPT_PE },
    { "mate-pair",       required_argument, NULL, OPT_MATEPAIR },
    { "help",            no_argument,       NULL, OPT_HELP },
    { "version",         no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
void parseScaffoldOptions(int argc, char** argv);

//
int scaffoldMain(int argc, char** argv)
{
    parseScaffoldOptions(argc, argv);

    int maxOverlap = 100;

    ScaffoldGraph graph;
    
    graph.loadVertices(opt::contigsFile, opt::minContigLength);

    if(!opt::peDistanceEstFile.empty())
        graph.loadDistanceEstimateEdges(opt::peDistanceEstFile, false);
    
    if(!opt::mateDistanceEstFile.empty())
        graph.loadDistanceEstimateEdges(opt::mateDistanceEstFile, true);

    assert(!opt::astatFile.empty());
    
    // Load the a-stat data and mark vertices as unique and repeat
    graph.loadAStatistic(opt::astatFile);
    ScaffoldAStatisticVisitor astatVisitor(opt::uniqueAstatThreshold, 
                                           opt::repeatAstatThreshold);
    graph.visit(astatVisitor);

    std::cout << "[sga-scaffold] Removing non-unique vertices from scaffold graph\n";
    graph.writeDot("pregraph.dot");
    graph.deleteVertices(SVC_REPEAT);
    graph.writeDot("scaffold.dot");
    
    // Remove polymorphic nodes from the graph
    ScaffoldPolymorphismVisitor polyVisitor(maxOverlap);
    while(graph.visit(polyVisitor)) {}

    graph.writeDot("nopoly-scaffold.dot");

    /*
    ScaffoldWalkVector outWalks;
    ScaffoldSearch::findPrimaryWalks(graph.getVertex("contig-14709"), graph.getVertex("contig-18695"), ED_SENSE, 10000, 1000, outWalks);
    //ScaffoldSearch::printWalks(outWalks);

    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        std::cout << "walk " << i << " " << outWalks[i].getContigLengthSum() << " " << outWalks[i].getGapSum() << "\n";
    }

    ScaffoldWalk& selectedWalk = outWalks[34];
    ScaffoldVertexPtrVector walkVertices = selectedWalk.getVertices();

    ScaffoldSearch::connectedComponents(&graph);

    for(size_t i = 0; i < walkVertices.size(); ++i)
    {
        walkVertices[i]->setClassification(SVC_REPEAT);
    }
    selectedWalk.printDot(std::cout);

    graph.writeDot("testPath.dot");

    exit(1);
    */
    /*
    // Compute the layout of the scaffolds
    ScaffoldLayoutVisitor layoutVisitor;
    graph.visit(layoutVisitor);
    graph.writeDot("afterLayout.dot");
    graph.visit(layoutVisitor);
    */
    // Break up any remaining scaffolds
    ScaffoldMultiEdgeRemoveVisitor cutVisitor;
    graph.visit(cutVisitor);

    // Write out the scaffolds
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
            case 'g': arg >> opt::asqgFile; break;
            case 'u': arg >> opt::uniqueAstatThreshold; break;
            case 'r': arg >> opt::repeatAstatThreshold; break;
            case 'o': arg >> opt::outFile; break;
            case OPT_PE: 
                arg >> opt::peDistanceEstFile; 
                break;
            case OPT_MATEPAIR: 
                arg >> opt::mateDistanceEstFile; 
                break;
            case OPT_HELP:
                std::cout << SCAFFOLD_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SCAFFOLD_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 2) 
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

    if(opt::contigsFile.empty())
    {
        std::cerr << SUBPROGRAM ": a contigs file must be provided\n";
        exit(1);
    }

    if(opt::peDistanceEstFile.empty() && opt::mateDistanceEstFile.empty())
    {
        std::cerr << SUBPROGRAM ": a distance estimation file must be provided using the --pe and/or --mate-pair options\n";
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

