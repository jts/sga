//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// assemble - convert read overlaps into contigs
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "assemble.h"
#include "StringGraph.h"

//
// Getopt
//
#define SUBPROGRAM "assemble"
static const char *ASSEMBLE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *ASSEMBLE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Create contigs for the reads in READSFILE. Overlaps are read from PREFIX.ovr. PREFIX defaults to the basename of READSFILE\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -p, --prefix=FILE                use PREFIX instead of the basename of READSFILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string readsFile;
	static std::string prefix;
}

static const char* shortopts = "p:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int assembleMain(int argc, char** argv)
{
	parseAssembleOptions(argc, argv);
	assemble();
	return 0;
}

void assemble()
{
	StringGraph* pGraph = createStringGraph(opt::readsFile, opt::prefix + ".ovr");

	SGFastaVisitor* bv = new SGFastaVisitor("before.fa");
	pGraph->visit(*bv);
	delete bv;
	
	pGraph->validate();
	pGraph->simplify();
	//pGraph->mergeVertices("tiled_1", "tiled_2");
	//pGraph->mergeVertices("tiled_0", "tiled_1");
	pGraph->validate();
	//pGraph->simplify();
	pGraph->writeDot("postSimplify.fa");

	SGFastaVisitor av("after.fa");
	pGraph->visit(av);
	delete pGraph;
}

// 
// Handle command line arguments
//
void parseAssembleOptions(int argc, char** argv)
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
				std::cout << ASSEMBLE_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << ASSEMBLE_VERSION_MESSAGE;
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
