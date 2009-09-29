//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// oview - view overlap alignments
//
#include <iostream>
#include <fstream>
#include <iterator>
#include "Util.h"
#include "oview.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "LCPArray.h"
#include "SGUtil.h"

static const char* DEFAULT_PADDING = "        ";

//
// Getopt
//
#define SUBPROGRAM "oview"
static const char *OVIEW_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *OVIEW_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Compute pairwise overlap between all the sequences in READS\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static std::string prefix;
	static std::string readsFile;
}

static const char* shortopts = "p:m:ve";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "exact",       no_argument,       NULL, 'e' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int oviewMain(int argc, char** argv)
{
	parseOviewOptions(argc, argv);

	// parse reads and index them
	ReadTable* pRT = new ReadTable(opt::readsFile);
	pRT->indexReadsByID();

	// read the contain map
	std::string containFile = opt::prefix + ".ctn";
	ContainMap containMap(containFile);

	// read the overlaps and draw
	std::string overlapFile = opt::prefix + ".ovr";
	std::ifstream overlapReader(overlapFile.c_str());
	Overlap o;

	while(overlapReader >> o)
	{
		// skip contained reads
		if(containMap.isContained(o.read[0].id) || containMap.isContained(o.read[1].id))
			continue;

		std::cout << "Overlap string: " << o << "\n\n";
		
		// Determine the left and right sequence
		SeqCoord* pLeftSC;
		SeqCoord* pRightSC;
		size_t leftIdx = (o.read[0].interval.start > o.read[1].interval.start) ? 0 : 1;
		pLeftSC = &o.read[leftIdx];
		pRightSC = &o.read[1 - leftIdx];

		std::string leftSeq = pRT->getRead(pLeftSC->id).seq;
		std::string rightSeq = pRT->getRead(pRightSC->id).seq;

		// Reverse the sequences if necessary
		if(pLeftSC->isReverse())
			leftSeq = reverseComplement(leftSeq);
		if(pRightSC->isReverse())
			rightSeq = reverseComplement(rightSeq);

		// Draw the left sequence
		std::cout << DEFAULT_PADDING << leftSeq << "\n";

		// Setup the padding string for the right sequence
		int offset = std::min(pLeftSC->interval.start, pLeftSC->interval.end);
		assert(offset >= 0 && offset < pLeftSC->seqlen);
		char* d_padding = new char[offset + 1];
		memset(d_padding, 32, offset);
		d_padding[offset] = 0;

		// Draw the right sequence
		std::cout << DEFAULT_PADDING << d_padding << rightSeq << "\n";
		std::cout << "\n";
		delete [] d_padding;
	}

	delete pRT;
	return 0;
}


// 
// Handle command line arguments
//
void parseOviewOptions(int argc, char** argv)
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
				std::cout << OVIEW_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << OVIEW_VERSION_MESSAGE;
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
		std::cout << "\n" << OVIEW_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::prefix.empty())
	{
		opt::prefix = stripFilename(opt::readsFile);
	}
}
