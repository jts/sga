//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// merge - merge read files and their indices
//
#include <iostream>
#include <fstream>
#include "SGACommon.h"
#include "Util.h"
#include "merge.h"
#include "SACAInducedCopying.h"
#include "BWTDiskConstruction.h"
#include "BWT.h"

//
void removeFiles(const std::string& inFile);

//
// Getopt
//
#define SUBPROGRAM "merge"

static const char *MERGE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *MERGE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READS1 READS2 ...\n"
"Merge the sequence files READS1, READS2, ... READSN and their associated indices into a single file/index\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"  -p, --prefix=PREFIX                  write final index to file using PREFIX (the default is to concatenate the input filenames)\n"
"  -r, --remove                         remove the original BWT, SAI and reads files after the merge\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string prefix;
    static bool bRemove;
}

static const char* shortopts = "p:m:vr";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "prefix",      required_argument, NULL, 'p' },
    { "remove",      no_argument,       NULL, 'r' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int mergeMain(int argc, char** argv)
{
    parseMergeOptions(argc, argv);
    StringVector inFiles;
    while(optind < argc)
    {
        inFiles.push_back(argv[optind++]);
    }
    assert(inFiles.size() == 2);
    if(inFiles[0] == inFiles[1])
        return 0; // avoid self-merge

    std::string prefix1 = stripFilename(inFiles[0]);
    std::string prefix2 = stripFilename(inFiles[1]);

    if(opt::prefix.empty())
    {
        std::string prefix1 = stripFilename(inFiles[0]);
        std::string prefix2 = stripFilename(inFiles[1]);
        opt::prefix = prefix1 + "." + prefix2;
    }

    // Merge the forward and reverse indices
    mergeIndependentIndices(inFiles[0], inFiles[1], opt::prefix, BWT_EXT, SAI_EXT, false);
    mergeIndependentIndices(inFiles[0], inFiles[1], opt::prefix, RBWT_EXT, RSAI_EXT, true);

    // Merge the read files
    mergeReadFiles(inFiles[0], inFiles[1], opt::prefix);

    if(opt::bRemove)
    {
        // Delete the original reads, bwt and sai files
        removeFiles(inFiles[0]);
        removeFiles(inFiles[1]);
    }
    return 0;
}

//
void removeFiles(const std::string& inFile)
{
    std::string prefix = stripFilename(inFile);
    std::string bwt_file = prefix + BWT_EXT;
    std::string rbwt_file = prefix + RBWT_EXT;

    std::string sai_file = prefix + SAI_EXT;
    std::string rsai_file = prefix + RSAI_EXT;

    unlink(bwt_file.c_str());
    unlink(rbwt_file.c_str());
    unlink(sai_file.c_str());
    unlink(rsai_file.c_str());
    unlink(inFile.c_str());
}

// 
// Handle command line arguments
//
void parseMergeOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 'r': opt::bRemove = true; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << MERGE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << MERGE_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 

    if (die) 
    {
        std::cerr << "Try `" << SUBPROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
}
