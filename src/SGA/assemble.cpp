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
#include "SGUtil.h"
#include "SGAlgorithms.h"

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
"      -b, --bubble                     perform bubble removal\n"
"      -t, --trim                       trim terminal branches\n"
"      -c, --close                      perform transitive closure\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string readsFile;
	static std::string prefix;
	static bool bTransClose;
	static bool bTrim;
	static bool bBubble;
}

static const char* shortopts = "p:vbtc";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "bubble",      no_argument,       NULL, 'b' },
	{ "trim",        no_argument,       NULL, 't' },
	{ "close",       no_argument,       NULL, 'c' },	
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
	StringGraph* pGraph = loadStringGraph(opt::readsFile, opt::prefix + ".ovr", opt::prefix + ".ctn");

	//pGraph->validate();
	pGraph->writeDot("before.dot");
	
	// Visitor functors
	SGTrimVisitor trimVisit;
	SGIslandVisitor islandVisit;
	SGTransRedVisitor trVisit;
	SGVariantVisitor varVisit;
	SGTCVisitor tcVisit;
	SGBubbleVisitor bubbleVisit;
	SGGraphStatsVisitor statsVisit;
	SGEdgeClassVisitor edgeClassVisit;

	// Pre-assembly graph stats
	std::cout << "Initial graph stats\n";
	pGraph->visit(statsVisit);

/*
	if(opt::bTrim)
	{
		std::cout << "Performing island/trim reduction\n";
		pGraph->visit(islandVisit);
		pGraph->visit(trimVisit);
		pGraph->visit(statsVisit);
	}
*/
	if(opt::bTransClose)
	{
		pGraph->visit(edgeClassVisit);
		std::cout << "\nPerforming transitive closure\n";
		int max_tc = 10;
		while(pGraph->visit(tcVisit) && --max_tc > 0);
		pGraph->visit(statsVisit);
		pGraph->visit(edgeClassVisit);
	}

	if(opt::bTrim)
	{
		std::cout << "\nPerforming island/trim reduction\n";
		pGraph->visit(islandVisit);
		pGraph->visit(trimVisit);
		pGraph->visit(trimVisit);
		pGraph->visit(trimVisit);
		pGraph->visit(statsVisit);
	}


	std::cout << "\nPerforming transitive reduction\n";
	pGraph->visit(trVisit);
	pGraph->simplify();
	pGraph->visit(statsVisit);

	//pGraph->visit(trimVisit);
	//pGraph->simplify();
	
	if(opt::bBubble)
	{
		std::cout << "\nPerforming bubble removal\n";
		// Bubble removal
		while(pGraph->visit(bubbleVisit))
			pGraph->simplify();
		pGraph->visit(statsVisit);
	}

	pGraph->simplify();

	std::cout << "\nFinal graph stats\n";
	pGraph->visit(statsVisit);

#ifdef VALIDATE
	VALIDATION_WARNING("SGA/assemble")
	pGraph->validate();
#endif

	// Write the results
	pGraph->writeDot("final.dot");
	SGFastaVisitor av("contigs.fa");
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
			case 'b': opt::bBubble = true; break;
            case 't': opt::bTrim = true; break;
			case 'c': opt::bTransClose = true; break;
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
