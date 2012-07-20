//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// beetl-convert - Utility to convert a beetl-constructed
// bwt index into SGA's format
//
#include <iostream>
#include <fstream>
#include "SGACommon.h"
#include "Util.h"
#include "convert-beetl.h"
#include "BWT.h"
#include "BWTWriterBinary.h"
#include "Timer.h"
#include "ReadInfoTable.h"
#include "SampledSuffixArray.h"

//
// Getopt
//
#define SUBPROGRAM "convert-beetl"

static const char *CONVERT_BEETL_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *CONVERT_BEETL_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... BWTFILE READSFILE\n"
"Convert the BWTFILE constructed from READSFILE into the SGA index format\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      --prefix=STR                     use STR as the prefix of the output files\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string inBWTFile;
    static std::string readsFile;
    static std::string prefix;
}

static const char* shortopts = "p:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NO_REVERSE };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "prefix",      required_argument, NULL, 'p' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int convertBeetlMain(int argc, char** argv)
{
    Timer t("sga beetl-convert");
    parseConvertBeetlOptions(argc, argv);
    std::cout << "Converting " << opt::inBWTFile << "\n";

    // Read the ASCII bwt beetl file and convert it to the SGA
    // run length encoded binary format

    // To write the header of the file, we need to count the number of strings
    // and symbols in the BWT. We do this in an initial pas
    std::istream* pReader = createReader(opt::inBWTFile);
    size_t numSymbols = 0;
    size_t numStrings = 0;
    char c;
    while(*pReader >> c) 
    {
        numSymbols += 1;
        if(c == 'N')
        {
            std::cerr << "Error: ambiguous character found in BWT\n";
            std::cerr << "sga preprocess must be run on the data\n";
            exit(EXIT_FAILURE);
        }

        if(c == '$')
            numStrings += 1;
    }
    delete pReader;

    printf("Read %zu symbols and %zu strings from the beetl bwt\n", numSymbols, numStrings);

    //
    std::string outBWTName = opt::prefix + BWT_EXT;
    BWTWriterBinary* pWriter = new BWTWriterBinary(outBWTName);
    pWriter->writeHeader(numStrings, numSymbols, BWF_NOFMI);
    
    // Re-read the file, writing out the bw chars
    pReader = createReader(opt::inBWTFile);
    while(*pReader >> c) 
        pWriter->writeBWChar(c);
    pWriter->finalize();
    delete pWriter;
    delete pReader;

    // Create the suffix array index files using the sampled suffix array machinery
    std::cout << "Generating lexicographic index (.sai)\n";
    BWT* pBWT = new BWT(outBWTName, 512); 
    SampledSuffixArray* pSSA = new SampledSuffixArray();

    pSSA->buildLexicoIndex(pBWT, opt::readsFile);
    pSSA->writeLexicoIndex(opt::prefix + SAI_EXT);

    delete pBWT;
    delete pSSA;
    return 0;
}

// 
// Handle command line arguments
//
void parseConvertBeetlOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << CONVERT_BEETL_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CONVERT_BEETL_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 2) 
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

    // Parse the input filenames
    opt::inBWTFile = argv[optind++];
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
        opt::prefix = stripFilename(opt::readsFile);
}
