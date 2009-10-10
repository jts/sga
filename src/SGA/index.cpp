//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// index - index reads using a suffix tree
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "index.h"
#include "SuffixArray.h"
#include "SeqReader.h"

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
"      -p, --prefix=PREFIX              write index to file using PREFIX instead of prefix of READSFILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string readsFile;
	static std::string prefix;
}

static const char* shortopts = "p:m:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "prefix",     required_argument, NULL, 'p' },
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
	buildIndex(opt::prefix + ".sa", pRT);
	
	return 0;
	/*
	// Create the reverse read table
	ReadTable* pRevRT = new ReadTable();
	pRevRT->initializeReverse(pRT);
	delete pRT; // done with the initial reads
	
	// Build the reverse suffix array
	buildIndex(opt::prefix + ".rsa", pRevRT);
	delete pRevRT;
	return 0;
	*/
}

void buildIndex(std::string outfile, const ReadTable* pRT)
{
	// Make initial suffix arrays
	SuffixArray* pSA = new SuffixArray();
	pSA->initialize(*pRT);
	pSA->sort(pRT);
	pSA->validate(pRT);
	if(opt::verbose > 0)
		pSA->print(pRT);
	writeSA(outfile, pSA);
	delete pSA;
	pSA = NULL;
}

void writeSA(std::string filename, const SuffixArray* pSA)
{
	std::ofstream out(filename.c_str());
	out << *pSA;
	out.close();
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
