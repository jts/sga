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
#include "AssembleExact.h"

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
"      -m, --min-overlap=LEN            minimum overlap required between two reads\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -d, --max-diff=D                 report all prefix-suffix matches that have at most D differences\n"
"      -i, --irreducible                only output the irreducible edges for each node\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static unsigned int minOverlap;
	static unsigned int maxDiff;
	static double errorRate;
	static std::string prefix;
	static std::string readsFile;
	static bool bIrreducibleOnly;
}

static const char* shortopts = "p:m:d:e:vi";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "max-diff",    required_argument, NULL, 'd' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "error-rate",  required_argument, NULL, 'e' },
	{ "irreducible", no_argument,       NULL, 'i' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int overlapMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("sga overlap");
	parseOverlapOptions(argc, argv);
	std::string hitsFile = computeHitsBWT();
	parseHits(hitsFile);
	delete pTimer;
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

	// Open the writer
	std::string hitsFile = opt::prefix + HITS_EXT;
	std::ofstream hitsHandle(hitsFile.c_str());
	assert(hitsHandle.is_open());

	// Initially, reserve enough room for 100 hits
	// Some reads may have (many) more hits so the vector will automatically expand to
	// contain them. That will be the memory highwater mark so there is no point deallocating
	// the vector after each trip through the loop.
	// The hits to the forward and reverse BWT must be kept seperately since they use different indices
	// if they were kept in the same structure there could be sa index collisions
	HitVector* pHits = new HitVector;
	HitVector* pRevHits = new HitVector;
	OverlapBlockList obOutputList;

	size_t count = 0;
	SeqReader reader(opt::readsFile);
	SeqItem read;
	Timer timer("BWT Alignment", true);
	int cost = 0;

	while(reader.get(read))
	{
		if(opt::verbose > 0 && count % 50000 == 0)
			printf("[overlap] Aligned %zu sequences\n", count);

		// We wish to collect all the prefix/suffix matches for this read. We first collect all the prefixes that
		// overlap the suffix of this read (including reverse complement reads) then find all the suffixes that match a prefix of this read.
		std::string seq = read.seq.toString();
		
		if(!opt::bIrreducibleOnly)
		{
			//cost += overlapReadExhaustive(count, read, pBWT, pRBWT, pHits, pRevHits);
		}
		else
		{
			cost += overlapReadIrreducible(read, pBWT, pRBWT, &obOutputList);
		}

		// Write the hits to the file
		if(!obOutputList.empty())
		{
			// Write the header info
			size_t numBlocks = obOutputList.size();
			hitsHandle << count << " " << numBlocks << " ";
			//std::cout << "<Wrote> idx: " << count << " count: " << numBlocks << "\n";
			for(OBLIter iter = obOutputList.begin(); iter != obOutputList.end(); ++iter)
			{
				OverlapBlockRecord record(*iter);
				hitsHandle << record << " ";
				//std::cout << "\t" << record << "\n";
			}
			hitsHandle << "\n";
			obOutputList.clear();
		}
		++count;
	}

	double align_time_secs = timer.getElapsedTime();
	printf("[bwt] aligned %zu sequences in %lfs (%lf sequences/s)\n", count, align_time_secs, (double)count / align_time_secs);
	printf("[bwt] performed %d iterations in the inner loop (%lf cost/sequence)\n", cost, (double)cost / (double)count);
	printf("[bwt] processed %d blocks, %d seeds (avg sl: %lf)\n", num_blocks, num_seeds, (double)total_seed_size / num_seeds);
	printf("[bwt] cost per seed: %lf\n", (double)cost / num_seeds);

	delete pHits;
	delete pRevHits;
	delete pBWT;
	delete pRBWT;

	hitsHandle.close();
	return hitsFile;
}

size_t overlapReadExhaustive(size_t index, SeqItem& read, const BWT* pBWT, const BWT* pRBWT, HitVector* pHits, HitVector* pRevHits, OverlapBlockList* /*pOBOut*/)
{
	size_t cost = 0;
	std::string seq = read.seq.toString();

	/*
	// Match the suffix of seq to prefixes
	hitTemplate.setRev(false, false);
	cost += BWTAlgorithms::alignSuffixInexact(seq, pBWT, pRBWT, opt::errorRate, opt::minOverlap, hitTemplate, pHits);

	hitTemplate.setRev(true, false);
	cost += BWTAlgorithms::alignSuffixInexact(complement(seq), pRBWT, pBWT, opt::errorRate, opt::minOverlap, hitTemplate, pHits);

	// Match the prefix of seq to suffixes
	hitTemplate.setRev(false, true);
	cost += BWTAlgorithms::alignSuffixInexact(reverseComplement(seq), pBWT, pRBWT, opt::errorRate, opt::minOverlap, hitTemplate, pRevHits);

	hitTemplate.setRev(true, true);
	cost += BWTAlgorithms::alignSuffixInexact(reverse(seq), pRBWT, pBWT, opt::errorRate, opt::minOverlap, hitTemplate, pRevHits);
	*/

	Hit hitTemplate(index, 0, 0, 0, false, false, 0); 

	hitTemplate.setRev(false, false);
	cost += BWTAlgorithms::alignSuffixExact(seq, pBWT, pRBWT, opt::minOverlap, hitTemplate, pHits);

	hitTemplate.setRev(true, false);
	cost += BWTAlgorithms::alignSuffixExact(complement(seq), pRBWT, pBWT, opt::minOverlap, hitTemplate, pRevHits);

	// Match the prefix of seq to suffixes
	hitTemplate.setRev(false, true);
	cost += BWTAlgorithms::alignSuffixExact(reverseComplement(seq), pBWT, pRBWT, opt::minOverlap, hitTemplate, pHits);

	hitTemplate.setRev(true, true);
	cost += BWTAlgorithms::alignSuffixExact(reverse(seq), pRBWT, pBWT, opt::minOverlap, hitTemplate, pRevHits);
	return cost;
}

// Construct the set of blocks describing irreducible overlaps with READ
// and write the blocks to pOBOut
size_t overlapReadIrreducible(SeqItem& read, const BWT* pBWT, const BWT* pRBWT, OverlapBlockList* pOBOut)
{
	// The complete set of overlap blocks are collected in obTemp
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	static OverlapBlockList obTemp;

	static AlignFlags sufPreAF(false, false, false);
	static AlignFlags prePreAF(false, true, true);
	static AlignFlags sufSufAF(true, false, true);
	static AlignFlags preSufAF(true, true, false);

	std::string seq = read.seq.toString();

	// Irreducible overlaps only
	WARN_ONCE("Irreducible-only assumptions: All reads are the same length")

	// Match the suffix of seq to prefixes
	BWTAlgorithms::findOverlapBlocks(seq, pBWT, pRBWT, opt::minOverlap, sufPreAF, &obTemp, pOBOut);
	BWTAlgorithms::findOverlapBlocks(complement(seq), pRBWT, pBWT, opt::minOverlap, prePreAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	BWTAlgorithms::calculateIrreducibleHits(&obTemp, pOBOut);
	obTemp.clear();

	// Match the prefix of seq to suffixes
	BWTAlgorithms::findOverlapBlocks(reverseComplement(seq), pBWT, pRBWT, opt::minOverlap, sufSufAF, &obTemp, pOBOut);
	BWTAlgorithms::findOverlapBlocks(reverse(seq), pRBWT, pBWT, opt::minOverlap, preSufAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	BWTAlgorithms::calculateIrreducibleHits(&obTemp, pOBOut);
	obTemp.clear();
	return 0;
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
		// If we are in irreducible mode don't filter hits, just output everything.
		if(curr_hit.saIdx != prevID || opt::bIrreducibleOnly)
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
	std::string line;
	while(getline(hitsHandle, line))
	{
		std::istringstream convertor(line);

		// Read the overlap blocks for a read
		size_t readIdx;
		size_t numBlocks;
		convertor >> readIdx >> numBlocks;

		//std::cout << "<Read> idx: " << readIdx << " count: " << numBlocks << "\n";
		for(size_t i = 0; i < numBlocks; ++i)
		{
			// Read the block
			OverlapBlockRecord record;
			convertor >> record;
			//std::cout << "\t" << record << "\n";
			Hit hit;
			hit.readIdx = readIdx;

			// Iterate through the range and write the overlaps
			for(int64_t j = record.range.lower; j <= record.range.upper; ++j)
			{
				const ReadTable* pCurrRT = (record.flags.isTargetRev()) ? pRevRT : pFwdRT;
				const SuffixArray* pCurrSAI = (record.flags.isTargetRev()) ? pRevSAI : pFwdSAI;
				const SeqItem& query = pCurrRT->getRead(hit.readIdx);

				hit.saIdx = j;
				hit.qstart = query.seq.length() - record.overlapLen;
				hit.len = record.overlapLen;
				hit.numDiff = record.numDiff;
				hit.targetRev = record.flags.isTargetRev();
				hit.queryRev = record.flags.isQueryRev();

				// The index of the second read is given as the position in the SuffixArray index
				const SeqItem& target = pCurrRT->getRead(pCurrSAI->get(hit.saIdx).getID());

				// Skip self alignments and non-canonical (where the query read has a lexo. higher name)
				if(query.id != target.id)
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
				
					// The alignment logic above has the potential to produce duplicate alignments
					// To avoid this, we skip overlaps where the id of the first coord is lexo. lower than 
					// the second or the match is a containment and the query is reversed (containments can be 
					// output up to 4 times total).
					// If we are running in irreducible mode this is not necessary
					if(o.id[0] < o.id[1] || (o.match.isContainment() && hit.queryRev))
						continue;
					writeOverlap(o, containHandle, overlapHandle);
				}
			}
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

// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
	// Set defaults
	opt::minOverlap = DEFAULT_MIN_OVERLAP;
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
			case 'i': opt::bIrreducibleOnly = true; break;
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
