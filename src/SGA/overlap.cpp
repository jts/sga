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
#include "LCPArray.h"
#include "SGACommon.h"

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
"      -e, --exact                      use the exact extraction algorithm instead of BWT-based alignment\n"
"      -m, --min-overlap=OVERLAP_LEN    minimum overlap required between two reads [30]\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static unsigned int minOverlap;
	static std::string prefix;
	static std::string readsFile;
	static bool bExactAlgo;
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
int overlapMain(int argc, char** argv)
{
	parseOverlapOptions(argc, argv);
	std::string hitsFile = computeHitsBWT();
	parseHits(hitsFile);
	return 0;
}

std::string computeHitsBWT()
{
	// Create the BWT
	BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
	BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);

	// Open the writers
	std::string hitsFile = opt::prefix + HITS_EXT;
	std::ofstream hitsHandle(hitsFile.c_str());
	assert(hitsHandle.is_open());

	// Initially, reserve enough room for 100 hits
	// Some reads may have (many) more hits so the vector will automatically expand to
	// contain them. That will be the memory highwater mark so there is no point deallocating
	// the vector after each trip through the loop.
	HitVector* pHits = new HitVector;
	pHits->reserve(100);

	size_t count = 0;
	SeqReader reader(opt::readsFile);
	SeqItem read;
	while(reader.get(read))
	{
		// Align the read and its reverse complement
		Sequence seqs[2];
		seqs[0] = read.seq.toString();
		seqs[1] = reverseComplement(seqs[0]);

		// Get all the hits of this sequence to the forward and reverse BWT
		for(size_t sn = 0; sn <= 1; ++sn)
		{
			bool isRC = (sn == 1) ? true : false;
			const Sequence& currSeq = seqs[sn];
			pBWT->getPrefixHits(count, currSeq, opt::minOverlap, false, isRC, pHits);
			pRBWT->getPrefixHits(count, reverse(currSeq), opt::minOverlap, true, !isRC, pHits);
		}
		++count;
		
		// Write the hits to the file
		for(size_t i = 0; i < pHits->size(); ++i)
			hitsHandle << (*pHits)[i] << "\n";
		pHits->clear();
	}

	delete pHits;
	delete pBWT;
	delete pRBWT;

	hitsHandle.close();
	return hitsFile;
}

// Parse all the hits and convert them to overlaps
void parseHits(std::string hitsFile)
{
	// Open files
	std::string overlapFile = opt::prefix + OVR_EXT;
	std::ofstream overlapHandle(overlapFile.c_str());
	assert(overlapHandle.is_open());

	std::string containFile = opt::prefix + CTN_EXT;
	std::ofstream containHandle(containFile.c_str());
	assert(containHandle.is_open());

	std::ifstream hitsHandle(hitsFile.c_str());
	assert(hitsHandle.is_open());

	// Load the suffix array index and the reverse suffix array index
	// Note these are not the full suffix arrays
	SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
	SuffixArray* pRevSAI = new SuffixArray(opt::prefix + RSAI_EXT);

	// Load the read tables
	ReadTable* pFwdRT = new ReadTable(opt::readsFile);
	ReadTable* pRevRT = new ReadTable();
	pRevRT->initializeReverse(pFwdRT);
	
	// Read each hit sequentially, converting it to an overlaps
	Hit hit;
	while(hitsHandle >> hit)
	{
		// Parse the hit, converting it into an overlap
		const ReadTable* pCurrRT = (hit.targetRev) ? pRevRT : pFwdRT;
		const SuffixArray* pCurrSAI = (hit.targetRev) ? pRevSAI : pFwdSAI;

		// Get the read names for the strings
		const SeqItem& query = pCurrRT->getRead(hit.readIdx);

		// The index of the second read is given as the position in the SuffixArray index
		const SeqItem& target = pCurrRT->getRead(pCurrSAI->get(hit.saIdx).getID());

		// Skip self alignments and non-canonical (where the query read has a lexo. higher name)
		if(query.id != target.id && query.id < target.id)
		{	
			// Compute the endpoints of the overlap
			int s1 = hit.qstart;
			int e1 = s1 + hit.len - 1;

			int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
			int e2 = s2 + hit.len - 1;
			
			// The alignment was to the reverse index, flip the coordinates of r2
			if(hit.targetRev)
				flipCoords(target.seq.length(), s2, e2);

			// The alignment was the reverse of the input read (complemented or otherwise), flip coordinates
			if(hit.queryRev)
				flipCoords(query.seq.length(), s1, e1);
			
			// If both intervals were reversed, swap the start and end to
			// indicate that the reads are in the same orientation
			if(hit.targetRev && hit.queryRev)
			{
				swap(s1, e1);
				swap(s2, e2);
			}
			Overlap o(query.id, s1, e1, query.seq.length(), target.id, s2, e2, target.seq.length());
			writeOverlap(o, containHandle, overlapHandle);
		}
	}

	// Delete allocated data
	delete pFwdSAI;
	delete pRevSAI;
	delete pFwdRT;
	delete pRevRT;

	// Close files
	overlapHandle.close();
	containHandle.close();
	hitsHandle.close();
}

// Before sanity checks on the overlaps and write them out
void writeOverlap(Overlap& ovr, std::ofstream& containHandle, std::ofstream& overlapHandle)
{
	// Ensure that the overlap is not a containment
	if(ovr.read[0].isContained() && ovr.read[1].isContained())
	{
		// Reads are mutually contained, consider the read with lexo. higher
		// id to be contained within the one with lexo. lower id
		if(ovr.read[0].id < ovr.read[1].id)
			writeContainment(containHandle, ovr.read[1].id, ovr.read[0].id);
		else
			writeContainment(containHandle, ovr.read[0].id, ovr.read[1].id);
		return;
	}
	else if(ovr.read[0].isContained())
	{
		writeContainment(containHandle, ovr.read[0].id, ovr.read[1].id);
		return;
	}
	else if(ovr.read[1].isContained())
	{
		writeContainment(containHandle, ovr.read[1].id, ovr.read[0].id);
		return;
	}

	// Unless both overlaps are extreme, skip
	if(!ovr.read[0].isExtreme() || !ovr.read[1].isExtreme())
	{
		std::cerr << "Skipping non-extreme overlap: " << ovr << "\n";
		return;
	}

	// Ensure that the overlaps are the correct orientation
	// If the reads are from the same strand, one should be left extreme and one should be right extreme
	// If they are from opposite strands, they should both be left (or right)
	bool sameStrand = !(ovr.read[0].isReverse() || ovr.read[1].isReverse());
	bool proper = false;
	if(sameStrand)
	{
		proper = (ovr.read[0].isLeftExtreme() != ovr.read[1].isLeftExtreme() && 
				  ovr.read[0].isRightExtreme() != ovr.read[1].isRightExtreme());
	}
	else
	{
		proper = (ovr.read[0].isLeftExtreme() == ovr.read[1].isLeftExtreme() && 
				  ovr.read[0].isRightExtreme() == ovr.read[1].isRightExtreme());
	}

	if(!proper)
	{
		std::cerr << "Skipping improper overlap: " << ovr << "\n";
		return;
	}
	// All checks passed, output the overlap
	overlapHandle << ovr << "\n";
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

#if 0
void computeOverlapsLCP()
{
	ReadTable* pRT = new ReadTable(opt::readsFile);
	ReadTable* pRevRT = new ReadTable();
	pRevRT->initializeReverse(pRT);

	// Load the suffix arrays
	SuffixArray* pSA = new SuffixArray(opt::prefix + ".sa");
	SuffixArray* pRSA = new SuffixArray(opt::prefix + ".rsa");

	// Open the writers
	std::string overlapFile = opt::prefix + ".ovr";
	std::ofstream overlapHandle(overlapFile.c_str());
	assert(overlapHandle.is_open());
	std::string containFile = opt::prefix + ".ctn";
	std::ofstream containHandle(containFile.c_str());
	assert(containHandle.is_open());

	//
	// Remove contained reads
	//
	
	// Detect all the contained reads in the forward suffix array, then use that list to remove them both
	SAElemPairVec containedReads = pSA->detectRedundantStrings(pRT);

	// Print the contained reads and build the set of suffix array ids to remove
	NumericIDSet idSet;
	for(size_t idx = 0; idx < containedReads.size(); ++idx)
	{
		std::string contained = pRT->getRead(containedReads[idx].first.getID()).id;
		std::string container = pRT->getRead(containedReads[idx].second.getID()).id;
		writeContainment(containHandle, contained, container);
		idSet.insert(containedReads[idx].first.getID());
	}

	pSA->removeReads(idSet);
	pRSA->removeReads(idSet);
	std::cerr << "Warning: validation is turned on\n";
	pSA->validate(pRT);
	pRSA->validate(pRevRT);

	if(opt::verbose > 0)
		pSA->print();

	//
	// Extract overlaps
	//
	
	OverlapVector overlapVec = pSA->extractPrefixSuffixOverlaps(opt::minOverlap, pRT);
	processOverlaps(overlapVec, containHandle, overlapHandle);

	// Cleanup
	overlapHandle.close();
	containHandle.close();
	delete pRT;
	delete pRevRT;
	delete pSA;
	delete pRSA;
}
#endif

// Write out a containmend
void writeContainment(std::ofstream& containHandle, const std::string& contained, const std::string& within)
{
	containHandle << contained << "\t" << within << "\n";
}

// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
	// Set defaults
	opt::minOverlap = 25;
	opt::bExactAlgo = false;

	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'm': arg >> opt::minOverlap; break;
			case 'p': arg >> opt::prefix; break;
			case 'e': opt::bExactAlgo = true; break;
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
