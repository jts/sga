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

// Set up the alignment blocks and call the alignment function on each block
int BWTAlgorithms::alignSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                       double error_rate, int minOverlap, Hit& hitTemplate, HitVector* pHits)
{
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

