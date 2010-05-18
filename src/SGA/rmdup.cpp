//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// rmdup - remove duplicated reads from the data set
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "rmdup.h"
#include "overlap.h"
#include "Timer.h"
#include "SGACommon.h"

//
// Getopt
//
#define SUBPROGRAM "rmdup"
static const char *RMDUP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *RMDUP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READFILE\n"
"Remove duplicate reads from the data set.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the output to FILE (default: READFILE.rmdup.fa)\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences identical\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string prefix;
    static std::string outFile;
    static std::string readsFile;
    static double errorRate;
}

static const char* shortopts = "p:o:e:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VALIDATE };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "prefix",         required_argument, NULL, 'p' },
    { "out",            required_argument, NULL, 'o' },
    { "error-rate",     required_argument, NULL, 'e' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int rmdupMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga rmdup");
    parseRmdupOptions(argc, argv);
    rmdup();
    delete pTimer;

    return 0;
}

void rmdup()
{
#if 0
    StringVector hitsFilenames;
    computeHitsBWT(OM_FULLREAD, NULL, hitsFilenames);
    parseDupHits(hitsFilenames);
#endif
}

void parseDupHits(const StringVector& hitsFilenames)
{
    // Load the suffix array index and the reverse suffix array index
    // Note these are not the full suffix arrays
    SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
    SuffixArray* pRevSAI = new SuffixArray(opt::prefix + RSAI_EXT);

    // Load the read table and output the initial vertex set, consisting of all the reads
    ReadTable* pFwdRT = new ReadTable(opt::readsFile);
    ReadTable* pRevRT = new ReadTable();
    pRevRT->initializeReverse(pFwdRT);

    // Convert the hits to overlaps and write them to the asqg file as initial edges
    for(StringVector::const_iterator iter = hitsFilenames.begin(); iter != hitsFilenames.end(); ++iter)
    {
        printf("[%s] parsing file %s\n", PROGRAM_IDENT, iter->c_str());
        std::istream* pReader = createReader(*iter);
    
        // Read each hit sequentially, converting it to an overlap
        std::string line;
        while(getline(*pReader, line))
        {
            OverlapVector ov = hitStringToOverlaps(line, pFwdRT, pRevRT, pFwdSAI, pRevSAI);
            (void)ov;
        }
        delete pReader;
    }

    // Delete allocated data
    delete pFwdSAI;
    delete pRevSAI;
    delete pFwdRT;
    delete pRevRT;
}

// 
// Handle command line arguments
//
void parseRmdupOptions(int argc, char** argv)
{
    // Set defaults
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 'e': arg >> opt::errorRate; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << RMDUP_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << RMDUP_VERSION_MESSAGE;
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

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }
}

