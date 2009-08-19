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
"      -o, --outfile=FILE               write overlaps to FILE [basename(READSFILE).sa]\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string readsFile;
	static std::string outFile;
}

static const char* shortopts = "o:m:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "outfile",     required_argument, NULL, 'o' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

int indexMain(int argc, char** argv)
{
	parseIndexOptions(argc, argv);
	std::cout << "Building index for " << opt::readsFile << "\n";
	buildIndex(opt::readsFile);
	return 0;
}

void buildIndex(std::string filename)
{
	ReadTable rt(filename);

	// Make initial suffix arrays
	SuffixArray sa;
	sa.initialize(rt);
	sa.sort(&rt);
	sa.validate(&rt);

	std::ofstream out(opt::outFile.c_str());
	out << sa;
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
			case 'o': arg >> opt::outFile; break;
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
	if(opt::outFile.empty())
	{
		opt::outFile = stripFilename(opt::readsFile) + ".sa";
	}
}
