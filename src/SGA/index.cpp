//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// index - index reads using a suffix tree
//
#include <iostream>
#include <fstream>
#include "SGACommon.h"
#include "Util.h"
#include "index.h"
#include "SuffixArray.h"
#include "SeqReader.h"
#include "saca.h"
#include "BWT.h"

//
// Getopt
//
#define SUBPROGRAM "index"

static const char *INDEX_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *INDEX_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Index the reads in READSFILE using a suffixarray/bwt\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"  -c, --check                          validate that the suffix array/bwt is correct\n"
"  -p, --prefix=PREFIX                  write index to file using PREFIX instead of prefix of READSFILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string readsFile;
	static std::string prefix;
	static bool validate;
}

static const char* shortopts = "p:m:cv";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "check",       no_argument,       NULL, 'c' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

int indexMain(int argc, char** argv)
{
	parseIndexOptions(argc, argv);
	std::cout << "Building index for " << opt::readsFile << "\n";

	// Parse the initial read table
	ReadTable* pRT = new ReadTable(opt::readsFile);
	
	// Create and write the suffix array for the forward reads
	buildIndex(opt::prefix, pRT, false);
	
	// Reverse all the reads
	pRT->reverseAll();

	// Build the reverse suffix array
	buildIndex(opt::prefix, pRT, true);
	
	delete pRT;
	return 0;
}

void buildIndex(std::string prefix, const ReadTable* pRT, bool isReverse)
{
	// Create suffix array from read table
	SuffixArray* pSA = new SuffixArray(pRT);
	BWT* pBWT = new BWT(pSA, pRT);

	if(opt::validate)
	{
		std::cout << "Validating suffix array\n";
		pSA->validate(pRT);
	}

	if(opt::verbose > 1)
	{
		//std::cout << "SuffixArray:\n";
		//pSA->print(pRT);
		std::cout << "BWT:\n";
		pBWT->print(pRT, pSA);
	}

	//std::string sa_filename = prefix + (!isReverse ? SA_EXT : RSA_EXT);
	//pSA->write(sa_filename);

	std::string bwt_filename = prefix + (!isReverse ? BWT_EXT : RBWT_EXT);
	pBWT->write(bwt_filename);

	std::string sufidx_filename = prefix + (!isReverse ? SAI_EXT : RSAI_EXT);
	pSA->writeIndex(sufidx_filename);

	delete pSA;
	pSA = NULL;
	delete pBWT;
	pBWT = NULL;
}

// 
// Handle command line arguments
//
void parseIndexOptions(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'p': arg >> opt::prefix; break;
			case '?': die = true; break;
			case 'c': opt::validate = true; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cout << INDEX_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << INDEX_VERSION_MESSAGE;
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
