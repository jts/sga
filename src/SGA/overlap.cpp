//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// overlap - compute pairwise overlaps between reads
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"
#include "BWT.h"

//
// Getopt
//
#define SUBPROGRAM "overlap"
static const char *OVERLAP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *OVERLAP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... INDEX READSFILE\n"
"Compute pairwise overlap between all the sequences in READS using the index INDEX\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -m, --min-overlap=OVERLAP_LEN    minimum overlap required between two reads [30]\n"
"      -o, --outfile=FILE               write overlaps to FILE [basename(READSFILE).ovr]\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static unsigned int minOverlap;
	static std::string indexFile;
	static std::string readsFile;
	static std::string outFile;
}

static const char* shortopts = "o:m:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "outfile",     required_argument, NULL, 'o' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int overlapMain(int argc, char** argv)
{
	parseOverlapOptions(argc, argv);
	computeOverlaps();
	return 0;
}

void computeOverlaps()
{
	ReadTable rt(opt::readsFile);

	// Load suffix array
	std::ifstream inSA(opt::indexFile.c_str());
	checkFileHandle(inSA, opt::indexFile);

	if(!inSA.is_open())
	{
		std::cerr << "Error: could not open " << opt::indexFile << " for read\n";
		exit(1);
	}

	SuffixArray sa;
	inSA >> sa;
	sa.validate(&rt);

	// Convert SA to a BWT
	BWT b(&sa, &rt);
	b.print(&rt);

	// Open the writer
	std::ofstream outHandle(opt::outFile.c_str());
	assert(outHandle.is_open());

	// Compute overlaps
	for(size_t i = 0; i < rt.getCount(); ++i)
	{
		HitVector hitVec = b.getOverlaps(rt.getRead(i).seq, 30);
		for(size_t j = 0; j < hitVec.size(); ++j)
		{
			// Convert the hit to an overlap
			Hit& hit = hitVec[j];
			
			// Skip self alignments
			if(hit.said.getID() != i)
			{
				// Get the read names for the strings
				std::string rn1 = rt.getRead(i).id;
				std::string rn2 = rt.getRead(hit.said.getID()).id;

				// Compute the endpoints of the overlap
				int s1 = hit.qstart;
				int e1 = s1 + hit.len - 1;

				int s2 = hit.said.getPos();
				int e2 = s2 + hit.len - 1;

				Overlap o(rn1, s1, e1, rn2, s2, e2);
				outHandle << o << "\n";
			}
		}
	}
	outHandle.close();
	inSA.close();
}

// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'm': arg >> opt::minOverlap; break;
			case 'o': arg >> opt::outFile; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cout << OVERLAP_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << OVERLAP_VERSION_MESSAGE;
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
	opt::indexFile = argv[optind++];
	opt::readsFile = argv[optind++];

	if(opt::outFile.empty())
	{
		opt::outFile = stripFilename(opt::readsFile) + ".ovr";
	}
}
