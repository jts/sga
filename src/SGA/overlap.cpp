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
#include "Timer.h"
#include "BWTAlgorithms.h"

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
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences aligned\n"
"      -m, --min-overlap=OVERLAP_LEN    minimum overlap required between two reads\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -d, --max-diff=D                 report all prefix-suffix matches that have at most D differences\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static unsigned int minOverlap;
	static unsigned int maxDiff;
	static double errorRate;
	static std::string prefix;
	static std::string readsFile;
	static bool bExactAlgo;
}

static const char* shortopts = "p:m:d:e:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "max-diff",    required_argument, NULL, 'd' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "error-rate",  required_argument, NULL, 'e' },
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

int total_hits = 0;
int output_hits = 0;
int num_blocks = 0;
int num_seeds = 0;
size_t total_seed_size = 0;

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
	// The hits to the forward and reverse BWT must be kept seperately since they use different indices
	HitVector* pHits = new HitVector;
	HitVector* pRevHits = new HitVector;

	pHits->reserve(100);
	pRevHits->reserve(100);

	size_t count = 0;
	SeqReader reader(opt::readsFile);
	SeqItem read;
	Timer timer("BWT Alignment", true);
	int cost = 0;
	while(reader.get(read))
	{
		if(opt::verbose > 0 && count % 50000 == 0)
			printf("[overlap] Aligned %zu sequences\n", count);

		// Align the read and its reverse complement
		Sequence seqs[2];
		seqs[0] = read.seq.toString();
		seqs[1] = reverseComplement(seqs[0]);
		// Get all the hits of this sequence to the forward and reverse BWT
		for(size_t sn = 0; sn <= 1; ++sn)
		{

			bool isRC = (sn == 1) ? true : false;
			const Sequence& currSeq = seqs[sn];
			if(opt::bExactAlgo)
			{
				pBWT->getPrefixHits(count, currSeq, opt::minOverlap, false, isRC, pHits);
				pRBWT->getPrefixHits(count, reverse(currSeq), opt::minOverlap, true, !isRC, pHits);
			}
			else
			{
				//cost += pBWT->getInexactPrefixHits(currSeq, pRBWT, opt::maxDiff, opt::minOverlap, count, false, isRC, pHits);
				//cost += pRBWT->getInexactPrefixHits(reverse(currSeq), pBWT, opt::maxDiff, opt::minOverlap, count, true, !isRC, pHits);
				/*
				Hit hitTemplate(count, 0, 0, 0, false, isRC, 0); 
				cost += alignSuffixMaxDiff(currSeq, pBWT, pRBWT, opt::maxDiff, opt::minOverlap, hitTemplate, pHits);
				hitTemplate.targetRev = true;
				hitTemplate.queryRev = !isRC;
				cost += alignSuffixMaxDiff(reverse(currSeq), pRBWT, pBWT, opt::maxDiff, opt::minOverlap, hitTemplate, pRevHits);
				*/
			
				Hit hitTemplate(count, 0, 0, 0, false, isRC, 0); 
				cost += BWTAlgorithms::alignSuffixInexact(currSeq, pBWT, pRBWT, opt::errorRate, opt::minOverlap, hitTemplate, pHits);
				hitTemplate.targetRev = true;
				hitTemplate.queryRev = !isRC;
				cost += BWTAlgorithms::alignSuffixInexact(reverse(currSeq), pRBWT, pBWT, opt::errorRate, opt::minOverlap, hitTemplate, pRevHits);
				
			}
		}

		// Write the hits to the file
		outputHits(hitsHandle, pHits);
		outputHits(hitsHandle, pRevHits);
		pHits->clear();
		pRevHits->clear();
		++count;
	}
	double align_time_secs = timer.getElapsedTime();
	printf("[bwt] aligned %zu sequences in %lfs (%lf sequences/s)\n", count, align_time_secs, (double)count / align_time_secs);
	printf("[bwt] performed %d iterations in the inner loop (%lf cost/sequence)\n", cost, (double)cost / (double)count);
	printf("[bwt] total hits: %d output hits: %d redundancy: %lf cost/hit: %lf\n", total_hits, output_hits, (double)total_hits / output_hits, (double)cost/output_hits);
	printf("[bwt] processed %d blocks, %d seeds (avg sl: %lf)\n", num_blocks, num_seeds, (double)total_seed_size / num_seeds);
	printf("[bwt] cost per seed: %lf\n", (double)cost / num_seeds);

	delete pHits;
	delete pRevHits;
	delete pBWT;
	delete pRBWT;

	hitsHandle.close();
	return hitsFile;
}

// Remove duplicate hits and write them to the filehandle
void outputHits(std::ofstream& handle, HitVector* pHits)
{
	// Some hits may be duplicate so only output the longest hit to a particular saIdx
	size_t prevID = std::numeric_limits<size_t>::max();
	size_t prevLen = 0;
	std::sort(pHits->begin(), pHits->end());
	total_hits += pHits->size();
	for(size_t i = 0; i < pHits->size(); ++i)
	{
		const Hit& curr_hit = (*pHits)[i];
		if(curr_hit.saIdx != prevID)
		{
			handle << curr_hit << "\n";
			++output_hits;
		}
		else
		{
			assert(curr_hit.len <= prevLen); 
		}

		prevID = curr_hit.saIdx;
		prevLen = curr_hit.len;
	}
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
	
	// Read each hit sequentially, converting it to an overlap
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

		//std::cout << "HIT: " << hit << " ids: " << query.id << " " << target.id << "\n";
		// Skip self alignments and non-canonical (where the query read has a lexo. higher name)
		if(query.id != target.id && query.id < target.id)
		{	
			// Compute the endpoints of the overlap
			int s1 = hit.qstart;
			int e1 = s1 + hit.len - 1;
			SeqCoord sc1(s1, e1, query.seq.length());

			int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
			int e2 = s2 + hit.len - 1;
			SeqCoord sc2(s2, e2, target.seq.length());

			// The coordinates are always with respect to the read, so flip them if
			// we aligned to/from the reverse of the read
			if(hit.queryRev)
				sc1.flip();
			if(hit.targetRev)
				sc2.flip();

			bool isRC = hit.targetRev != hit.queryRev;

			Overlap o(query.id, sc1, target.id, sc2, isRC, hit.numDiff);
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
	if(ovr.match.coord[0].isContained() || ovr.match.coord[1].isContained())
	{
		containHandle << ovr << "\n";
		return;
	}

	// Unless both overlaps are extreme, skip
	if(!ovr.match.coord[0].isExtreme() || !ovr.match.coord[1].isExtreme())
	{
		std::cerr << "Skipping non-extreme overlap: " << ovr << "\n";
		return;
	}

	bool sameStrand = !ovr.match.isRC();
	bool proper = false;
	
	if(sameStrand)
	{
		proper = (ovr.match.coord[0].isLeftExtreme() != ovr.match.coord[1].isLeftExtreme() && 
				  ovr.match.coord[0].isRightExtreme() != ovr.match.coord[1].isRightExtreme());
	}
	else
	{
		proper = (ovr.match.coord[0].isLeftExtreme() == ovr.match.coord[1].isLeftExtreme() && 
				  ovr.match.coord[0].isRightExtreme() == ovr.match.coord[1].isRightExtreme());
	}
	
	if(!proper)
	{
		std::cerr << "Skipping improper overlap: " << ovr << "\n";
		return;
	}

	// All checks passed, output the overlap
	overlapHandle << ovr << "\n";
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

// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
	// Set defaults
	opt::minOverlap = DEFAULT_MIN_OVERLAP;
	opt::bExactAlgo = false;
	opt::maxDiff = 0;

	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'm': arg >> opt::minOverlap; break;
			case 'p': arg >> opt::prefix; break;
			case 'e': arg >> opt::errorRate; break;
			case 'd': arg >> opt::maxDiff; break;
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
