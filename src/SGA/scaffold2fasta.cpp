//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// scaffold2fasta - Convert a scaffold file into a fasta
// file representing the scaffold with gaps
//
#include <getopt.h>
#include "Util.h"
#include "config.h"
#include "ScaffoldGraph.h"
#include "ScaffoldVisitors.h"
#include "SGUtil.h"

#define SUBPROGRAM "scaffold2fasta"
static const char *SCAFFOLD2FASTA_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SCAFFOLD2FASTA_USAGE_MESSAGE =
"Usage: " SUBPROGRAM " [OPTION] ... [-f CONTIGSFILE | -a ASQGFILE] SCAFFOLDFILE\n"
"Write out a fasta file for the scaffolds indicated in SCAFFOLDFILE\n"
"One of -f CONTIGSFILE or -a ASQGFILE must be provided. If an asqg file is provided,\n"
"the program will attempt to determine the sequence linking the scaffold components by\n"
"walking the graph/\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -f, --contig-file=FILE           read the contig sequences from FILE\n"
"      -a, --asqg-file=FILE             read the sequence graph from FILE. This supercedes --contig-file\n"
"          --no-singletons              do not output scaffolds that consist of a single contig\n"
"      -o, --outfile=FILE               write the scaffolds to FILE (default: scaffolds.fa)\n"
"      -m, --min-length=N               only output scaffolds longer than N bases\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string contigFile;
    static std::string asqgFile;
    static std::string outFile = "scaffolds.fa";
    static std::string scafFile;
    static bool bNoSingletons = false;
    static int minScaffoldLength = 0;
}

static const char* shortopts = "vm:o:f:a:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NOSINGLETON };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "min-length",     required_argument, NULL, 'm' },
    { "outfile",        required_argument, NULL, 'o' },
    { "contig-file",    required_argument, NULL, 'f' },
    { "asqg-file",      required_argument, NULL, 'a' },
    { "no-singleton",   no_argument,       NULL, OPT_NOSINGLETON },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
void parseScaffold2fastaOptions(int argc, char** argv);

//
int scaffold2fastaMain(int argc, char** argv)
{
    parseScaffold2fastaOptions(argc, argv);

    if(opt::asqgFile.empty())
        assert(false && "only asqg file is implemented atm");

    std::cout << "Reading graph from " << opt::asqgFile << "\n";
    std::cout << "Reading scaffolds from " << opt::scafFile << "\n";

    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
    delete pGraph;
    return 0;
}

//
void parseScaffold2fastaOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'm': arg >> opt::minScaffoldLength; break;
            case 'f': arg >> opt::contigFile; break;
            case 'a': arg >> opt::asqgFile; break;
            case 'o': arg >> opt::outFile; break;
            case OPT_NOSINGLETON: opt::bNoSingletons = false; break;
            case OPT_HELP:
                std::cout << SCAFFOLD2FASTA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SCAFFOLD2FASTA_VERSION_MESSAGE;
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
    opt::scafFile = argv[optind++];

    if(opt::contigFile.empty() && opt::asqgFile.empty())
    {
        std::cerr << SUBPROGRAM ": a contigs file or asqg file must be provided\n";
        exit(1);
    }
}

