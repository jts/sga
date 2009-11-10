//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// exact - test algorithms for assembling exact data
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "exact.h"
#include "SGUtil.h"
#include "AssembleExact.h"

//
// Getopt
//
#define SUBPROGRAM "exact"
static const char *EXACT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *EXACT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Assemble the data in READSFILE assuming the reads are exact representations of substrings of the source genome\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -p, --prefix=FILE                use PREFIX instead of the basename of READSFILE\n"
"      -m, --min-overlap=OVERLAP_LEN    minimum overlap required between two reads\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string readsFile;
	static std::string prefix;
	static unsigned int minOverlap;
}

static const char* shortopts = "p:m:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int exactMain(int argc, char** argv)
{
	parseExactOptions(argc, argv);
	exact();
	return 0;
}

void exact()
{
	BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
	BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
	ReadTable* pRT = new ReadTable(opt::readsFile);


	delete pBWT;
	delete pRBWT;
	delete pRT;
}

// 
// Handle command line arguments
//
void parseExactOptions(int argc, char** argv)
{
	opt::minOverlap = DEFAULT_MIN_OVERLAP;
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'p': arg >> opt::prefix; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case 'm': arg >> opt::minOverlap; break;
			case OPT_HELP:
				std::cout << EXACT_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << EXACT_VERSION_MESSAGE;
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
