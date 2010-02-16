//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
#include "OverlapAlgorithm.h"
#include <math.h>

#define SWAP_LIST(x, y) pSwap = (x); (x) = (y); (y) = pSwap;

// Perform the overlap
void OverlapAlgorithm::overlapRead(const SeqItem& read, OverlapBlockList* pOutList) const
{
	if(!m_bIrreducible)
	{
		//overlapReadExhaustive(read, pOutList);
		overlapReadInexact(read, pOutList);
	}
	else
	{
		 overlapReadIrreducible(read, pOutList);
	}
}

//
void OverlapAlgorithm::overlapReadExhaustive(const SeqItem& read, OverlapBlockList* pOBOut) const
{
	// Collect the complete set of overlaps in pOBOut
	static const AlignFlags sufPreAF(false, false, false);
	static const AlignFlags prePreAF(false, true, true);
	static const AlignFlags sufSufAF(true, false, true);
	static const AlignFlags preSufAF(true, true, false);

	std::string seq = read.seq.toString();

	// Match the suffix of seq to prefixes
	findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, pOBOut, pOBOut);
	findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, pOBOut, pOBOut);

	// Match the prefix of seq to suffixes
	findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, pOBOut, pOBOut);
	findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, pOBOut, pOBOut);
}

void OverlapAlgorithm::overlapReadInexact(const SeqItem& read, OverlapBlockList* pOBOut) const
{
	// Collect the complete set of overlaps in pOBOut
	static const AlignFlags sufPreAF(false, false, false);
	static const AlignFlags prePreAF(false, true, true);
	static const AlignFlags sufSufAF(true, false, true);
	static const AlignFlags preSufAF(true, true, false);

	std::string seq = read.seq.toString();

	// Match the suffix of seq to prefixes
	overlapSuffixInexact(seq, m_pBWT, m_pRevBWT, m_errorRate, m_minOverlap, sufPreAF, pOBOut, pOBOut);
	overlapSuffixInexact(complement(seq), m_pRevBWT, m_pBWT, m_errorRate, m_minOverlap, prePreAF, pOBOut, pOBOut);

	// Match the prefix of seq to suffixes
	overlapSuffixInexact(reverseComplement(seq), m_pBWT, m_pRevBWT, m_errorRate, m_minOverlap, sufSufAF, pOBOut, pOBOut);
	overlapSuffixInexact(reverse(seq), m_pRevBWT, m_pBWT, m_errorRate, m_minOverlap, preSufAF, pOBOut, pOBOut);
}

// Construct the set of blocks describing irreducible overlaps with READ
// and write the blocks to pOBOut
void OverlapAlgorithm::overlapReadIrreducible(const SeqItem& read, OverlapBlockList* pOBOut) const
{
	// The complete set of overlap blocks are collected in obTemp
	// The filtered set (containing only irreducible overlaps) are placed into pOBOut
	// by calculateIrreducibleHits
	OverlapBlockList obTemp;

	static const AlignFlags sufPreAF(false, false, false);
	static const AlignFlags prePreAF(false, true, true);
	static const AlignFlags sufSufAF(true, false, true);
	static const AlignFlags preSufAF(true, true, false);

	std::string seq = read.seq.toString();

	// Irreducible overlaps only
	WARN_ONCE("Irreducible-only assumptions: All reads are the same length")

	// Match the suffix of seq to prefixes
	findOverlapBlocks(seq, m_pBWT, m_pRevBWT, m_minOverlap, sufPreAF, &obTemp, pOBOut);
	findOverlapBlocks(complement(seq), m_pRevBWT, m_pBWT, m_minOverlap, prePreAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);
	obTemp.clear();

	// Match the prefix of seq to suffixes
	findOverlapBlocks(reverseComplement(seq), m_pBWT, m_pRevBWT, m_minOverlap, sufSufAF, &obTemp, pOBOut);
	findOverlapBlocks(reverse(seq), m_pRevBWT, m_pBWT, m_minOverlap, preSufAF, &obTemp, pOBOut);

	// Process the first set of blocks and output the irreducible hits to pHits
	calculateIrreducibleHits(m_pBWT, m_pRevBWT, &obTemp, pOBOut);
}

// Write overlap blocks out to a file
void OverlapAlgorithm::writeOverlapBlocks(size_t readIdx, const OverlapBlockList* pList, std::ofstream& writer) const
{
	// Write the hits to the file
	if(!pList->empty())
	{
		// Write the header info
		size_t numBlocks = pList->size();
		writer << readIdx << " " << numBlocks << " ";
		//std::cout << "<Wrote> idx: " << count << " count: " << numBlocks << "\n";
		for(OverlapBlockList::const_iterator iter = pList->begin(); iter != pList->end(); ++iter)
		{
			writer << *iter << " ";
		}
		writer << "\n";
	}
}

// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
// overlaps with a suffix of w. The ranges are added to the pOBList
void OverlapAlgorithm::findOverlapBlocks(const std::string& w, const BWT* pBWT,
									  const BWT* pRevBWT, int minOverlap,
									  const AlignFlags& af, OverlapBlockList* pOBList,
                                      OverlapBlockList* pOBFinal)
{
	// All overlaps are added to this list and then sub-maximal overlaps are removed
	OverlapBlockList workingList;

	// The algorithm is as follows:
	// We perform a backwards search using the FM-index for the string w.
	// As we perform the search we collect the intervals 
	// of the significant prefixes (len >= minOverlap) that overlap w.
	BWTIntervalPair ranges;
	size_t l = w.length();
	int start = l - 1;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);
	
	// Collect the OverlapBlocks
	for(int i = start - 1; i >= 1; --i)
	{
		// Compute the range of the suffix w[i, l]
		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);
		int overlapLen = l - i;
		if(overlapLen >= minOverlap)
		{
			// Calculate which of the prefixes that match w[i, l] are terminal
			// These are the proper prefixes (they are the start of a read)
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			
			// The probe interval contains the range of proper prefixes
			if(probe.interval[1].isValid())
			{
				assert(probe.interval[1].lower > 0);
				workingList.push_back(OverlapBlock(probe, overlapLen, 0, af));
			}
		}
	}

	// Determine if this sequence is contained and should not be processed further
	BWTAlgorithms::updateBothL(ranges, w[0], pBWT);

	// Ranges now holds the interval for the full-length read
	// To handle containments, we output the overlapBlock to the final overlap block list
	// and it will be processed later
	// Two possible containment cases:
	// 1) This read is a substring of some other read
	// 2) This read is identical to some other read
	
	// Case 1 is indicated by the existance of a non-$ left or right hand extension
	// In this case we return no alignments for the string
	AlphaCount left_ext = BWTAlgorithms::getExtCount(ranges.interval[0], pBWT);
	AlphaCount right_ext = BWTAlgorithms::getExtCount(ranges.interval[1], pRevBWT);
	if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
	{
		// This case isn't handled yet
		pOBList->clear();
		return;
	}
	else
	{
		BWTAlgorithms::updateBothL(ranges, '$', pBWT);
		if(ranges.isValid())
			pOBFinal->push_back(OverlapBlock(ranges, w.length(), 0, af));
	}

	// Remove sub-maximal OverlapBlocks and move the remainder to the output list
	removeSubMaximalBlocks(&workingList);
	pOBList->splice(pOBList->end(), workingList);
	return;
}

// Set up the alignment blocks and call the alignment function on each block
int OverlapAlgorithm::overlapSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                           double error_rate, int minOverlap, const AlignFlags& af, 
										   OverlapBlockList* pOBList, OverlapBlockList* pOBFinal)
{
	WARN_ONCE("Using inexact overlap!");
	if((int)w.length() < minOverlap)
	{
		std::cerr << "Warning string " << w << " is shorter than minOverlap, it will not be aligned\n";
		return 0;
	}

	// Calculate the number of blocks and the size of each block
	int len = static_cast<int>(w.length());
	
	// Calculate the shortest overlap s.t. 1 difference in the overlap region is less than error_rate
	int max_differences = (int)floor(error_rate * (double)len);

	// Align the reads in segments of shortest_overlap
	// Each segment allows 1 more error than the previous
	int cost = 0;
	int covered = 0;
	for(int i = 0; i <= max_differences; ++i)
	{
		// Calculate the minimum amount of overlap s.t. i differences is less than error_rate
		int min_overlap_size = std::max((int)ceil(i / error_rate), minOverlap);

		// Calculate the endpoint of the block, it is the amount of overlap s.t. i + 1 is less than the error rate
		int max_overlap_size = std::min((int)ceil((i + 1) / error_rate) - 1, len);
		// If the max overlap in this block is less than the minOverlap parameter, skip the block
		if(max_overlap_size < minOverlap)
			continue;

		assert(double(i) / min_overlap_size <= error_rate);
		assert(double(i) / max_overlap_size <= error_rate);
		assert(double(i+1) / max_overlap_size > error_rate);

		// Convert the overlap sizes to block positions, which start from the end of w
		int block_start = len - max_overlap_size;
		int block_end = len - min_overlap_size;
		covered += (block_end - block_start + 1);
		cost += _alignSegmentSimple(w, block_start, block_end, pBWT, pRevBWT, i, af, pOBList, pOBFinal);
	}
	assert(covered == (len - minOverlap) + 1);
	return cost;
}

// Seeded blockwise BWT alignment of prefix-suffix for reads
// Each alignment is given a seed region and a block region
// The seed region is the terminal portion of w where maxDiff + 1 seeds are created
// at least 1 of these seeds must align exactly for there to be an alignment with 
// at most maxDiff differences between the prefix/suffix. Only alignments within the
// range [block_start, block_end] are output. The block_end coordinate is inclusive
int OverlapAlgorithm::_alignSegment(const std::string& w, int block_start, int block_end,
                                     const BWT* pBWT, const BWT* pRevBWT, int maxDiff, 
								     const AlignFlags& af, OverlapBlockList* pOBList, 
    								 OverlapBlockList* pOBFinal)
{
	OverlapSeedList* pLeftList = new OverlapSeedList;
	OverlapSeedList* pRightList = new OverlapSeedList;
	OverlapSeedList* pTempList = new OverlapSeedList;
	OverlapSeedList* pSwap;
	OverlapBlockList workingList;

	OverlapSeedList::iterator iter;
	int len = w.length();
	int seed_start = block_end;
	int seed_end = len;
	int num_seeds = maxDiff + 1;
	int seed_len = (seed_end - seed_start) / num_seeds;

	//printf("srs: %d sre: %d bs: %d be: %d d: %d\n", seed_start, seed_end, block_start, block_end, maxDiff);

	// Populate the initial seeds
	for(int i = 0; i < num_seeds; ++i)
	{
		int pos = seed_start + i*seed_len;
		OverlapSeed align;
		align.left_index = pos;
		align.right_index = pos;
		align.dir = ED_RIGHT;
		align.z = maxDiff;
		align.seed_len = std::min(seed_len, len - pos);
		assert(pos < len);

		//printf("Creating seed at %d to %d, str: %s\n", align.left_index, align.left_index + align.seed_len, w.substr(pos, align.seed_len).c_str());

		// Initialize the left and right suffix array intervals
		char b = w[pos];
		BWTAlgorithms::initIntervalPair(align.ranges, b, pBWT, pRevBWT);		
		pTempList->push_back(align);
	}

	// Extend all the seeds to the right without mismatches, then move them to the right extension list
	iter = pTempList->begin();
	while(iter != pTempList->end())
	{
		OverlapSeed& align = *iter;
		bool valid = true;
		while(align.isSeed())
		{
			++align.right_index;
			char b = w[align.right_index];
			BWTAlgorithms::updateBothR(align.ranges, b, pRevBWT);
			if(!align.isIntervalValid(RIGHT_INT_IDX))
			{
				valid = false;
				break;
			}
		}

		if(valid)
			pRightList->push_back(align);
		pTempList->erase(iter++);
	}
	assert(pTempList->empty());

	// Extend the right and left lists in lockstep, allowing mismatches
	while(!pLeftList->empty() || !pRightList->empty())
	{
		bool switched_list = false;

		// Process right extensions first
		iter = pRightList->begin();
		while(iter != pRightList->end())
		{
			OverlapSeed& align = *iter;
			//printf("ProcessingRIGHT: "); align.print(w);

			// If this alignment has run all the way to the end of the sequence
			// switch it to be a left extension sequence
			if(align.right_index == len - 1)
			{
				align.dir = ED_LEFT;
				pLeftList->push_back(align);
				switched_list = true;
			}
			else
			{
				++align.right_index;
				// If the length of the alignment is less than the seed length, do not allow mismatches
				if(align.z == 0)
				{
					char b = w[align.right_index];
					BWTAlgorithms::updateBothR(align.ranges, b, pRevBWT);
					if(align.isIntervalValid(RIGHT_INT_IDX))
						pTempList->push_back(align);
				}
				else
				{
					for(int i = 0; i < 4; ++i)
					{
						char b = ALPHABET[i];
						OverlapSeed branched = align;
						BWTAlgorithms::updateBothR(branched.ranges, b, pRevBWT);
						if(branched.isIntervalValid(RIGHT_INT_IDX))
						{
							if(b != w[align.right_index])
								--branched.z;
							pTempList->push_back(branched);
						}
					}
				}
			}
			pRightList->erase(iter++);
		}

		SWAP_LIST(pTempList, pRightList)
		pTempList->clear();

		// attempt to merge seeds
		if(switched_list)
		{
			pLeftList->sort(OverlapSeed::compareLeftRange);
			pLeftList->unique(OverlapSeed::equalLeftRange);

			/*
			for(iter = pLeftList->begin(); iter != pLeftList->end(); ++iter)
			{
				std::cout << "Align: "; iter->print();
			}
			*/
		}

		// Process left extensions
		iter = pLeftList->begin();
		while(iter != pLeftList->end())
		{
			OverlapSeed& align = *iter;
			//printf("ProcessingLEFT: "); align.print(w);

			// If the alignment is within the current block, attempt to output matching prefixes
			if(align.left_index <= block_end)
			{
				BWTIntervalPair probe = align.ranges;
				BWTAlgorithms::updateBothL(probe, '$', pBWT);
				
				// The probe interval contains the range of proper prefixes
				if(probe.interval[1].isValid())
				{
					assert(probe.interval[1].lower > 0);
					int overlapLen = len - align.left_index;
					OverlapBlock nBlock(OverlapBlock(probe, len - align.left_index, maxDiff - align.z, af));

					if(overlapLen == len)
						pOBFinal->push_back(nBlock);
					else
						workingList.push_back(nBlock);
				}

				/*
				int64_t t_lower = pBWT->getPC('$') + pBWT->getOcc('$', align.ranges.interval[LEFT_INT_IDX].lower - 1);
				int64_t t_upper = pBWT->getPC('$') + pBWT->getOcc('$', align.ranges.interval[LEFT_INT_IDX].upper) - 1;

				for(int64_t sa_idx = t_lower; sa_idx <= t_upper; ++sa_idx)
				{
					hitTemplate.saIdx = sa_idx;
					hitTemplate.qstart = align.left_index;
					hitTemplate.len = len - align.left_index;
					hitTemplate.numDiff = maxDiff - align.z;
					pHits->push_back(hitTemplate);
				}
				*/
			}

			// Extend hits
			--align.left_index;
			if(align.left_index >= block_start)
			{
				// If there cannot be a branch, only process the matching base
				if(align.z == 0)
				{
					char b = w[align.left_index];
					BWTAlgorithms::updateInterval(align.ranges.interval[LEFT_INT_IDX], b, pBWT);
					if(align.isIntervalValid(LEFT_INT_IDX))
						pTempList->push_back(align);
				}
				else
				{
					for(int i = 0; i < 4; ++i)
					{
						char b = ALPHABET[i];
						OverlapSeed branched = align;
						// Only update left interval
						BWTAlgorithms::updateInterval(branched.ranges.interval[LEFT_INT_IDX], b, pBWT);

						if(branched.isIntervalValid(LEFT_INT_IDX))
						{
							if(ALPHABET[i] != w[align.left_index])
								--branched.z;
							pTempList->push_back(branched);
						}
					}
				}
			}
			pLeftList->erase(iter++);
		}
		SWAP_LIST(pTempList, pLeftList)
		pTempList->clear();
	}

	// Remove sub-maximal OverlapBlocks and move the remainder to the output list
	removeSubMaximalBlocks(&workingList);
	pOBList->splice(pOBList->end(), workingList);

	delete pLeftList;
	delete pRightList;
	delete pTempList;
	return 0;
}

// Seeded blockwise BWT alignment of prefix-suffix for reads
// Each alignment is given a seed region and a block region
// The seed region is the terminal portion of w where maxDiff + 1 seeds are created
// at least 1 of these seeds must align exactly for there to be an alignment with 
// at most maxDiff differences between the prefix/suffix. Only alignments within the
// range [block_start, block_end] are output. The block_end coordinate is inclusive
int OverlapAlgorithm::_alignSegmentSimple(const std::string& w, int block_start, int block_end,
                                           const BWT* pBWT, const BWT* pRevBWT, int maxDiff, 
										   const AlignFlags& af, OverlapBlockList* pOBList, 
										   OverlapBlockList* pOBFinal)
{
	OverlapSeedList* pCurrList = new OverlapSeedList;
	OverlapSeedList* pNextList = new OverlapSeedList;
	OverlapBlockList blockWorkingList;
	OverlapSeedList::iterator iter;

	createOverlapSeeds(w, pBWT, pRevBWT, block_start, block_end, maxDiff, pCurrList);
	extendSeedsExactRight(w, pBWT, pRevBWT, ED_RIGHT, pCurrList, pNextList);

	pCurrList->clear();
	pCurrList->swap(*pNextList);
	assert(pNextList.empty());

	while(!pCurrList->empty())
	{
		iter = pCurrList->begin();
		while(iter != pCurrList->end())
		{
			OverlapSeed& align = *iter;
			
			if(align.dir == ED_RIGHT)
			{
				extendSeedInexactRight(align, w, pBWT, pRevBWT, ED_RIGHT, pNextList);
			}
			else
			{
                extendSeedInexactLeft(align, w, pBWT, pRevBWT, ED_LEFT, block_start, block_end, maxDiff,
				                      af, pNextList, &blockWorkingList, pOBFinal);

			}
			pCurrList->erase(iter++);
		}
		pCurrList->swap(*pNextList);
	}

	// Remove sub-maximal OverlapBlocks and move the remainder to the output list
	removeSubMaximalBlocks(&blockWorkingList);
	pOBList->splice(pOBList->end(), blockWorkingList);

	delete pCurrList;
	delete pNextList;
	return 0;
}

void OverlapAlgorithm::createOverlapSeeds(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                          int /*block_start*/, int block_end, int maxDiff, 
										  OverlapSeedList* pOutList)
{
	int len = w.length();
	int seed_start = block_end;
	int seed_end = len;
	int num_seeds = maxDiff + 1;
	int seed_len = (seed_end - seed_start) / num_seeds;

	//printf("srs: %d sre: %d bs: %d be: %d d: %d\n", seed_start, seed_end, block_start, block_end, maxDiff);

	// Populate the initial seeds
	for(int i = 0; i < num_seeds; ++i)
	{
		int pos = seed_start + i*seed_len;
		OverlapSeed align;
		align.left_index = pos;
		align.right_index = pos;
		align.dir = ED_RIGHT;
		align.z = maxDiff;
		align.seed_len = std::min(seed_len, len - pos);
		assert(pos < len);

		//printf("Creating seed at %d to %d, str: %s\n", align.left_index, align.left_index + align.seed_len, w.substr(pos, align.seed_len).c_str());

		// Initialize the left and right suffix array intervals
		char b = w[pos];
		BWTAlgorithms::initIntervalPair(align.ranges, b, pBWT, pRevBWT);		
		pOutList->push_back(align);
	}
}

void OverlapAlgorithm::extendSeedsExactRight(const std::string& w, const BWT* /*pBWT*/, const BWT* pRevBWT,
                                             ExtendDirection /*dir*/, const OverlapSeedList* pInList, 
											 OverlapSeedList* pOutList)
{

	for(OverlapSeedList::const_iterator iter = pInList->begin(); iter != pInList->end(); ++iter)
	{
		OverlapSeed align = *iter;
		bool valid = true;
		while(align.isSeed())
		{
			++align.right_index;
			char b = w[align.right_index];
			BWTAlgorithms::updateBothR(align.ranges, b, pRevBWT);
			if(!align.isIntervalValid(RIGHT_INT_IDX))
			{
				valid = false;
				break;
			}
		}
		
		if(valid)
			pOutList->push_back(align);
	}
}

void OverlapAlgorithm::extendSeedInexactRight(OverlapSeed& seed, const std::string& w, const BWT* /*pBWT*/, 
                                              const BWT* pRevBWT, ExtendDirection /*dir*/,
											  OverlapSeedList* pOutList)
{
	// If this alignment has run all the way to the end of the sequence
	// switch it to be a left extension sequence
	if(seed.right_index == int(w.length() - 1))
	{
		seed.dir = ED_LEFT;
		pOutList->push_back(seed);
	}
	else
	{
		++seed.right_index;
		
		if(seed.z == 0)
		{
			char b = w[seed.right_index];
			BWTAlgorithms::updateBothR(seed.ranges, b, pRevBWT);
			if(seed.isIntervalValid(RIGHT_INT_IDX))
				pOutList->push_back(seed);
		}
		else
		{
			for(int i = 0; i < 4; ++i)
			{
				char b = ALPHABET[i];
				OverlapSeed branched = seed;
				WARN_ONCE("UPDATEBOTHL????");
				BWTAlgorithms::updateBothR(branched.ranges, b, pRevBWT);
				if(branched.isIntervalValid(RIGHT_INT_IDX))
				{
					if(b != w[seed.right_index])
						--branched.z;
					pOutList->push_back(branched);
				}
			}
		}
	}
}

void OverlapAlgorithm::extendSeedInexactLeft(OverlapSeed& seed, const std::string& w, const BWT* pBWT, 
                                                 const BWT* /*pRevBWT*/, ExtendDirection /*dir*/,
                                                 int block_start, int block_end, int maxDiff, const AlignFlags& af,
												 OverlapSeedList* pOutList, OverlapBlockList* pOBPartialList, 
												 OverlapBlockList* pOBFullList)
{
	//printf("ProcessingLEFT: "); align.print(w);

	// If the alignment is within the current block, attempt to output matching prefixes
	// this implies the overlap is long enough to be valid
	if(seed.left_index <= block_end)
	{
		BWTIntervalPair probe = seed.ranges;
		BWTAlgorithms::updateBothL(probe, '$', pBWT);
		
		// The probe interval contains the range of proper prefixes
		if(probe.interval[1].isValid())
		{
			assert(probe.interval[1].lower > 0);
			int len = w.length();
			int overlapLen = len - seed.left_index;
			OverlapBlock nBlock(OverlapBlock(probe, overlapLen, maxDiff - seed.z, af));

			if(overlapLen == len)
				pOBFullList->push_back(nBlock);
			else
				pOBPartialList->push_back(nBlock);
		}
	}

	// Extend hits
	--seed.left_index;
	if(seed.left_index >= block_start)
	{
		if(seed.z == 0)
		{
			// Extend exact
			char b = w[seed.left_index];
			BWTAlgorithms::updateBothL(seed.ranges, b, pBWT);
			if(seed.isIntervalValid(LEFT_INT_IDX))
				pOutList->push_back(seed);
		}
		else
		{
			for(int i = 0; i < 4; ++i)
			{
				char b = ALPHABET[i];
				OverlapSeed branched = seed;
				// Only update left interval
				BWTAlgorithms::updateBothL(seed.ranges, b, pBWT);
				if(branched.isIntervalValid(LEFT_INT_IDX))
				{
					if(ALPHABET[i] != w[seed.left_index])
						--branched.z;
					pOutList->push_back(branched);
				}
			}
		}
	}
}


// Calculate the irreducible hits from the vector of OverlapBlocks
void OverlapAlgorithm::calculateIrreducibleHits(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList* pOBList, OverlapBlockList* pOBFinal)
{
	// processIrreducibleBlocks requires the pOBList to be sorted in descending order
	pOBList->sort(OverlapBlock::sortSizeDescending);
	processIrreducibleBlocks(pBWT, pRevBWT, *pOBList, pOBFinal);
}

// iterate through obList and determine the overlaps that are irreducible. This function is recursive.
// The final overlap blocks corresponding to irreducible overlaps are written to pOBFinal.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
// Invariant: each block corresponds to the same extension of the root sequence w.
void OverlapAlgorithm::processIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, OverlapBlockList* pOBFinal)
{
	if(obList.empty())
		return;
	
	// Count the extensions in the top level (longest) blocks first
	int topLen = obList.front().overlapLen;
	AlphaCount ext_count;
	OBLIter iter = obList.begin();
	while(iter != obList.end() && iter->overlapLen == topLen)
	{
		ext_count += iter->getCanonicalExtCount(pBWT, pRevBWT);
		++iter;
	}
	
	// Three cases:
	// 1) The top level block has ended, it contains the extension $. Output TLB and end.
	// 2) There is a singular unique extension base for all the blocks. Update all blocks and recurse.
	// 3) There are multiple extension bases, branch and recurse.
	// If some block other than the TLB ended, it must be contained within the TLB and it is not output
	// or considered further. Containments are handled elsewhere.
	// Likewise if multiple distinct strings in the TLB ended, we only output the top one. The rest
	// must have the same sequence as the top one and are hence considered to be contained with the top element.
	if(ext_count.get('$') > 0)
	{
		// An irreducible overlap has been found. It is possible that there are two top level blocks
		// (one in the forward and reverse direction). Since we can't decide which one
		// contains the other at this point, we output hits to both. Under a fixed 
		// length string assumption one will be contained within the other and removed later.
		OBLIter tlbIter = obList.begin();
		while(tlbIter != obList.end() && tlbIter->overlapLen == topLen)
		{
			pOBFinal->push_back(OverlapBlock(*tlbIter));
			++tlbIter;
		} 
		return;
	}
	else
	{
		// Count the rest of the blocks
		while(iter != obList.end())
		{
			ext_count += iter->getCanonicalExtCount(pBWT, pRevBWT);
			++iter;
		}

		if(ext_count.hasUniqueDNAChar())
		{
			// Update all the blocks using the unique extension character
			// This character is in the canonical representation wrt to the query
			char b = ext_count.getUniqueDNAChar();
			updateOverlapBlockRangesRight(pBWT, pRevBWT, obList, b);
			return processIrreducibleBlocks(pBWT, pRevBWT, obList, pOBFinal);
		}
		else
		{
			for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
			{
				char b = ALPHABET[idx];
				if(ext_count.get(b) > 0)
				{
					OverlapBlockList branched = obList;
					updateOverlapBlockRangesRight(pBWT, pRevBWT, branched, b);
					processIrreducibleBlocks(pBWT, pRevBWT, branched, pOBFinal);
				}
			}
		}
	}
}

// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
void OverlapAlgorithm::updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, char b)
{
	OverlapBlockList::iterator iter = obList.begin(); 
	while(iter != obList.end())
	{
		char cb = iter->flags.isQueryComp() ? complement(b) : b;
		BWTAlgorithms::updateBothR(iter->ranges, cb, iter->getExtensionBWT(pBWT, pRevBWT));
		// remove the block from the list if its no longer valid
		if(!iter->ranges.isValid())
			iter = obList.erase(iter);
		else
			++iter;
	}
}

