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
#include <iterator>
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
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Compute pairwise overlap between all the sequences in READS\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -m, --min-overlap=OVERLAP_LEN    minimum overlap required between two reads [30]\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static unsigned int minOverlap;
	static std::string prefix;
	static std::string readsFile;
}

static const char* shortopts = "p:m:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "prefix",      required_argument, NULL, 'p' },
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
	ReadTable* pRT = new ReadTable(opt::readsFile);
	ReadTable* pRevRT = new ReadTable();
	pRevRT->initializeReverse(pRT);

	// Load the suffix arrays
	SuffixArray* pSA = loadSuffixArray(opt::prefix + ".sa");
	SuffixArray* pRSA = loadSuffixArray(opt::prefix + ".rsa");

	// Create the BWT
	BWT* pBWT = createBWT(pSA, pRT);
	BWT* pRBWT = createBWT(pRSA, pRevRT);

	// Open the writer
	std::string outFile = opt::prefix + ".ovr";
	std::ofstream outHandle(outFile.c_str());
	assert(outHandle.is_open());

	// Compute overlaps
	for(size_t i = 0; i < pRT->getCount(); ++i)
	{
		const SeqItem& read = pRT->getRead(i);

		// Align the read and its reverse complement
		Sequence seqs[2];
		seqs[0] = read.seq;//reverseComplement(read.seq);
		seqs[1] = reverseComplement(seqs[0]);

		HitData hits;
		// Get all the hits of this sequence to the forward and reverse BWT
		for(size_t sn = 0; sn <= 1; ++sn)
		{
			bool isRC = (sn == 1) ? true : false;
			const Sequence& currSeq = seqs[sn];
			
			pBWT->getHits(currSeq, opt::minOverlap, false, isRC, &hits);
			pRBWT->getHits(reverse(currSeq), opt::minOverlap, true, !isRC, &hits);
		}

		HitVector hitVec = hits.getHits();
		OverlapVector overlaps = processHits(i, hitVec, pRT, pRevRT);
		std::copy(overlaps.begin(), overlaps.end(), std::ostream_iterator<Overlap>(outHandle, "\n"));

		//	OverlapVector fwdOvr = alignRead(i, currSeq, pBWT, pRT, false, isRC);
		//	OverlapVector revOvr = alignRead(i, reverse(currSeq), pRBWT, pRevRT, true, !isRC);

		//	std::copy(fwdOvr.begin(), fwdOvr.end(), std::ostream_iterator<Overlap>(outHandle, "\n"));
		//	std::copy(revOvr.begin(), revOvr.end(), std::ostream_iterator<Overlap>(outHandle, "\n"));

	}
	outHandle.close();
	
	delete pRT;
	delete pRevRT;
	delete pSA;
	delete pRSA;
	delete pBWT;
	delete pRBWT;
}

// Process all the hits into overlaps
OverlapVector processHits(size_t seqIdx, const HitVector& hitVec, const ReadTable* pFwdRT, const ReadTable* pRevRT)
{
	OverlapVector overlaps;
	// Convert the hits to overlaps
	for(size_t j = 0; j < hitVec.size(); ++j)
	{
		// Convert the hit to an overlap
		const Hit& hit = hitVec[j];
		
		const ReadTable* pCurrRT = (hit.targetRev) ? pRevRT : pFwdRT;

		// Get the read names for the strings
		SeqItem query = pCurrRT->getRead(seqIdx);
		SeqItem target = pCurrRT->getRead(hit.said.getID());
		
		// Skip self alignments and non-canonical (where the query read has a lexo. higher name)
		if(query.id != target.id && query.id < target.id)
		{	
			std::cout << "Query: " << query.id << " align flags: " << hit.targetRev << "," << hit.queryRev << "\n";
			// Compute the endpoints of the overlap
			int s1 = hit.qstart;
			int e1 = s1 + hit.len - 1;

			int s2 = hit.said.getPos();
			int e2 = s2 + hit.len - 1;
			
			// The alignment was to the reverse index, flip the coordinates of r2
			if(hit.targetRev)
			{
				flipCoords(target.seq.length(), s2, e2);
			}

			// The alignment was the reverse of the input read (complemented or otherwise), flip coordinates
			if(hit.queryRev)
			{
				flipCoords(query.seq.length(), s1, e1);
			}
			
			// If both intervals were reversed, swap the start and end to
			// indicate that the reads are in the same orientation
			if(hit.targetRev && hit.queryRev)
			{
				swap(s1, e1);
				swap(s2, e2);
			}

			Overlap o(query.id, s1, e1, target.id, s2, e2);
			overlaps.push_back(o);
		}
	}
	return overlaps;
}

//
void flipCoords(const int len, int& s, int &e)
{
	s = len - (s + 1);
	e = len - (e + 1);	
}

//
void swap(int& s, int& e)
{
	int temp = e;
	e = s;
	s = temp;
}

// Create a bwt from a suffix array file and read table
BWT* createBWT(SuffixArray* pSA, const ReadTable* pRT)
{
	// Convert SA to a BWT
	BWT* pBWT = new BWT(pSA, pRT);
	//pBWT->print(pRT);
	return pBWT;
}

//
SuffixArray* loadSuffixArray(std::string filename)
{
	std::ifstream inSA(filename.c_str());
	checkFileHandle(inSA, filename);

	SuffixArray* pSA = new SuffixArray();
	inSA >> *pSA;
	return pSA;
}


// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
	opt::minOverlap = 25;
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'm': arg >> opt::minOverlap; break;
			case 'p': arg >> opt::prefix; break;
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
		std::cout << "\n" << OVERLAP_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::prefix.empty())
	{
		opt::prefix = stripFilename(opt::readsFile);
	}
}
