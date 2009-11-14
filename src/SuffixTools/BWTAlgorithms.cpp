//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// bwt_algorithms.cpp - Algorithms for aligning to a bwt structure
//
#include "BWTAlgorithms.h"
#include <math.h>

#define SWAP_LIST(x, y) pSwap = (x); (x) = (y); (y) = pSwap;

// Find the interval in pBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
BWTInterval BWTAlgorithms::findInterval(const BWT* pBWT, const std::string& w)
{
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	BWTInterval interval;
	initInterval(interval, curr, pBWT);
	--j;

	for(;j >= 0; --j)
	{
		curr = w[j];
		updateInterval(interval, curr, pBWT);
		if(!interval.isValid())
			return interval;
	}
	return interval;
}

// Find the intervals in pBWT/pRevBWT corresponding to w
// If w does not exist in the BWT, the interval 
// coordinates [l, u] will be such that l > u
BWTIntervalPair BWTAlgorithms::findIntervalPair(const BWT* pBWT, const BWT* pRevBWT, const std::string& w)
{
	BWTIntervalPair intervals;	
	int len = w.size();
	int j = len - 1;
	char curr = w[j];
	initIntervalPair(intervals, curr, pBWT, pRevBWT);
	--j;

	for(;j >= 0; --j)
	{
		curr = w[j];
		updateBothL(intervals, curr, pBWT);
		if(!intervals.isValid())
			return intervals;
	}
	return intervals;
}

// Set up the alignment blocks and call the alignment function on each block
int BWTAlgorithms::alignSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                       double error_rate, int minOverlap, Hit& hitTemplate, HitVector* pHits)
{
	if((int)w.length() < minOverlap)
	{
		std::cerr << "Warning string " << w << " is shorter than minOverlap, it will not be aligned\n";
		return 0;
	}

	// Calculate the number of blocks and the size of each block
	int len = static_cast<int>(w.length());
	
	// Calculate the shortest overlap s.t. 1 difference in the overlap region is less than error_rate
	int max_differences = (int)floor(error_rate * (double)len);

	// Align the reads in blocks of shortest_overlap
	// Each block allows 1 more error than the previous
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
		cost += _alignBlock(w, block_start, block_end, pBWT, pRevBWT, i, hitTemplate, pHits);
		extern int num_blocks;
		num_blocks++;
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
int BWTAlgorithms::_alignBlock(const std::string& w, int block_start, int block_end, const BWT* pBWT, const BWT* pRevBWT, 
                int maxDiff, Hit& hitTemplate, HitVector* pHits)
{
	BWTAlignList* pLeftList = new BWTAlignList;
	BWTAlignList* pRightList = new BWTAlignList;
	BWTAlignList* pTempList = new BWTAlignList;
	BWTAlignList* pSwap;

	BWTAlignList::iterator iter;
	int cost = 0;
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
		BWTAlign align;
		align.left_index = pos;
		align.right_index = pos;
		align.dir = ED_RIGHT;
		align.z = maxDiff;
		align.seed_len = std::min(seed_len, len - pos);
		assert(pos < len);

		//printf("Creating seed at %d to %d, str: %s\n", align.left_index, align.left_index + align.seed_len, w.substr(pos, align.seed_len).c_str());

		// Initialize the left and right suffix array intervals
		char b = w[pos];
		initIntervalPair(align.ranges, b, pBWT, pRevBWT);		
		pTempList->push_back(align);
		extern int num_seeds;
		extern size_t total_seed_size;
		++num_seeds;
		total_seed_size += align.seed_len;
	}


	// Extend all the seeds to the right, without mismatches, then move them to the right extension list
	iter = pTempList->begin();
	while(iter != pTempList->end())
	{
		BWTAlign& align = *iter;
		bool valid = true;
		while(align.isSeed())
		{
			++cost;
			++align.right_index;
			char b = w[align.right_index];
			updateBothR(align.ranges, b, pRevBWT);
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
			++cost;
			BWTAlign& align = *iter;
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
					updateBothR(align.ranges, b, pRevBWT);
					if(align.isIntervalValid(RIGHT_INT_IDX))
						pTempList->push_back(align);
				}
				else
				{
					for(int i = 0; i < 4; ++i)
					{
						char b = ALPHABET[i];
						BWTAlign branched = align;
						updateBothR(branched.ranges, b, pRevBWT);
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
			pLeftList->sort(BWTAlign::compareLeftRange);
			pLeftList->unique(BWTAlign::equalLeftRange);

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
			++cost;
			BWTAlign& align = *iter;
			//printf("ProcessingLEFT: "); align.print(w);

			// If the alignment is within the current block, attempt to output matching prefixes
			if(align.left_index <= block_end)
			{
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
			}

			// Extend hits
			--align.left_index;
			if(align.left_index >= block_start)
			{
				// If there cannot be a branch, only process the matching base
				if(align.z == 0)
				{
					char b = w[align.left_index];
					updateInterval(align.ranges.interval[LEFT_INT_IDX], b, pBWT);
					if(align.isIntervalValid(LEFT_INT_IDX))
						pTempList->push_back(align);
				}
				else
				{
					for(int i = 0; i < 4; ++i)
					{
						char b = ALPHABET[i];
						BWTAlign branched = align;
						// Only update left interval
						updateInterval(branched.ranges.interval[LEFT_INT_IDX], b, pBWT);

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

	delete pLeftList;
	delete pRightList;
	delete pTempList;
	return cost;
}

// Perform an exact overlapping alignment
int BWTAlgorithms::alignSuffixExact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                      int minOverlap, Hit& hitTemplate, HitVector* pHits)
{
	// The algorithm is as follows:
	// We perform a backwards search using the FM-index for the string w.
	// As we perform the search we collect the intervals matching proper prefixes and output them as hits
	BWTIntervalPair ranges;
	size_t l = w.length();
	int start = l - 1;
	int cost = 0;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);
	
	// Collect the overlaps
	for(int i = start - 1; i >= 0; --i)
	{
		++cost;
		// Compute the range of the suffix w[i, l]
		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);
		int overlapLen = l - i;
		if(overlapLen >= minOverlap)
		{
			// Calculate which of the prefixes that match w[i, l] are terminal
			// These are the proper prefixes (they are the start of a read)
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			
			for(int64_t sa_idx = probe.interval[0].lower; sa_idx <= probe.interval[0].upper; ++sa_idx)
			{
				hitTemplate.saIdx = sa_idx;
				hitTemplate.qstart = i;
				hitTemplate.len = overlapLen;
				hitTemplate.numDiff = 0;
				pHits->push_back(hitTemplate);
			}
		}
	}
	return cost;
}

void BWTAlgorithms::alignSuffixExactIrreducible(const std::string& w, const BWT* pBWT, const BWT* pRevBWT,
             	                           int minOverlap, Hit& hitTemplate, HitVector* pHits)
{
	// The algorithm is as follows:
	// We perform a backwards search using the FM-index for the string w.
	// As we perform the search we collect the intervals 
	// of the significant prefixes (len >= minOverlap) that overlap w.
	// Once all the intervals have been found, we extend them one base at a time and output the 
	// overlaps that are irreducible. An overlap between x,y is transitive if there is some other non-contained
	// read z, that overlaps x and the unmatched portion of z is a prefix of the unmatched portion of y. The 
	// read z is considered to be irreducible wrt x and is output. 
	OverlapBlockList obList;
	BWTIntervalPair ranges;
	size_t l = w.length();
	int start = l - 1;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);
	
	// Collect the blocks of overlaps
	// Do not use full-length overlaps (which must be containments) so 
	// stop at 1.
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
				obList.push_front(OverlapBlock(probe, overlapLen));
			}
		}
	}

	// Determine if this sequence is contained and should not be processed further
	BWTAlgorithms::updateBothL(ranges, w[0], pBWT);

	// Two cases in which this read is contained in another and should not be processed further
	// 1) It is a substring of a some other read/string
	// 2) It is identical to another read/string that has a lexographically lower ID.
	
	// Case 1 is indicated by the existance of a non-$ left or right hand extension
	// In this case we return no alignments for the string (could do better?)
	AlphaCount left_ext = getExtCount(ranges.interval[0], pBWT);
	AlphaCount right_ext = getExtCount(ranges.interval[1], pRevBWT);
	if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
	{
		std::cout << "found substring " << w << " " << hitTemplate.readIdx << "\n";
		return;
	}
	else
	{
		// Case 2 is handle in a slightly more tricky way. We output a hit to the head of the block of full-length strings
		// If the current read is the head of the block of full length strings, it will be considered a self-alignment later
		// and removed. All others will be considered proper containments since the head of the block will have the lexographically
		// lowest ID
		BWTAlgorithms::updateBothL(ranges, '$', pBWT);
		if(ranges.isValid())
		{
			size_t sa_idx = ranges.interval[0].lower; // index of the head
			hitTemplate.saIdx = sa_idx;
			hitTemplate.qstart = 0;
			hitTemplate.len = w.length();
			hitTemplate.numDiff = 0;
			pHits->push_back(hitTemplate);
		}
	}

	// pOBList contains a list of all the intervals in the FM-index
	// that are valid prefixs. Determine the irreducible edges.
	return processIrreducibleBlocks(obList, w.length(), pRevBWT, hitTemplate, pHits);
}

// iterate through obList and determine the overlaps that are irreducible. This function is recursive.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
// Invariant: each block corresponds to the same extension of the root sequence w.
void BWTAlgorithms::processIrreducibleBlocks(OverlapBlockList& obList, const size_t q_len, const BWT* pRevBWT, Hit& hitTemplate, HitVector* pHits)
{
	if(obList.empty())
		return;
	
	OverlapBlockList::iterator iter = obList.begin(); 
	AlphaCount top_count = getExtCount(iter->ranges.interval[1], pRevBWT);

	AlphaCount remain_count;
	++iter;
	for(; iter != obList.end(); ++iter)
	{
		remain_count += getExtCount(iter->ranges.interval[1], pRevBWT);
	}
	AlphaCount total_count = top_count + remain_count;

	// Three cases:
	// 1) The top level block has ended, it contains the extension $. Output TLB and end.
	// 2) There is a singular unique extension base for all the blocks. Update all blocks and recurse.
	// 3) There are multiple extension bases, branch and recurse.
	// If some block other than the TLB ended, it must be contained within the TLB and it is not output
	// or considered further. Containments are handled elsewhere.
	// Likewise if multiple distinct strings in the TLB ended, we only output the top one. The rest
	// must have the same sequence as the top one and are hence considered to be contained with the top element.
	if(top_count.get('$') > 0)
	{
		// Found the irreducible extension of the blocks, return it in pHits
		OverlapBlock& tlb = obList.front();
		size_t sa_idx = tlb.ranges.interval[0].lower;
		hitTemplate.saIdx = sa_idx;
		hitTemplate.qstart = q_len - tlb.overlapLen;
		hitTemplate.len = tlb.overlapLen;
		hitTemplate.numDiff = 0;
		pHits->push_back(hitTemplate);
		return;
	}
	else if(total_count.hasUniqueDNAChar())
	{
		char b = total_count.getUniqueDNAChar();
		updateOverlapBlockRangesRight(obList, b, pRevBWT);
		return processIrreducibleBlocks(obList, q_len, pRevBWT, hitTemplate, pHits);
	}
	else
	{
		for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
		{
			char b = ALPHABET[idx];
			if(total_count.get(b) > 0)
			{
				OverlapBlockList branched = obList;
				updateOverlapBlockRangesRight(branched, b, pRevBWT);
				processIrreducibleBlocks(branched, q_len, pRevBWT, hitTemplate, pHits);
			}
		}
	}
}

// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
void BWTAlgorithms::updateOverlapBlockRangesRight(OverlapBlockList& obList, char b, const BWT* pRevBWT)
{
	OverlapBlockList::iterator iter = obList.begin(); 
	while(iter != obList.end())
	{
		BWTAlgorithms::updateBothR(iter->ranges, b, pRevBWT);
		// remove the block from the list if its no longer valid
		if(!iter->ranges.isValid())
			iter = obList.erase(iter);
		else
			++iter;
	}
}



