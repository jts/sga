//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// bwt_algorithms.cpp - Algorithms for aligning to a bwt structure
//
#include "bwt_algorithms.h"
#include <math.h>

#define LEFT_INT_IDX 0
#define RIGHT_INT_IDX 1

inline void updateLeft(BWTAlign& align, char b, const BWT* pB)
{
	size_t pb = pB->getC(b);
	align.r_lower[LEFT_INT_IDX] = pb + pB->getOcc(b, align.r_lower[LEFT_INT_IDX] - 1);
	align.r_upper[LEFT_INT_IDX] = pb + pB->getOcc(b, align.r_upper[LEFT_INT_IDX]) - 1;

}

inline void _initInterval(BWTAlign& align, int idx, char b, const BWT* pB)
{
	align.r_lower[idx] = pB->getC(b);
	align.r_upper[idx] = align.r_lower[idx] + pB->getOcc(b, pB->getBWLen() - 1) - 1;
}

// Update both the left and right intervals using the reverse BWT
inline void updateBoth(BWTAlign& align, char b, const BWT* pRevBWT)
{
	//_updateInterval(align, LEFT_INT_IDX, b, pBWT);
	//_updateInterval(align, RIGHT_INT_IDX, b, pRevBWT);

	// Update the left index using the difference between the AlphaCounts in the reverse table
	AlphaCount diff = pRevBWT->getOccDiff(align.r_lower[1] - 1, align.r_upper[1]);
	align.r_lower[0] = align.r_lower[0] + diff.getLessThan(b);
	align.r_upper[0] = align.r_lower[0] + diff.get(b) - 1;

	// Update the right index directly
	size_t pb = pRevBWT->getC(b);
	align.r_lower[1] = pb + pRevBWT->getOcc(b, align.r_lower[1] - 1);
	align.r_upper[1] = pb + pRevBWT->getOcc(b, align.r_upper[1]) - 1;

}

inline void initIntervals(BWTAlign& align, char b, const BWT* pBWT, const BWT* pRevBWT)
{
	_initInterval(align, LEFT_INT_IDX, b, pBWT);
	_initInterval(align, RIGHT_INT_IDX, b, pRevBWT);
}


// Set up the alignment blocks and call the alignment function on each block
int alignSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
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
	}
	assert(covered == (len - minOverlap) + 1);
	return cost;
}

// Set up the alignment blocks and call the alignment function on each block
int alignSuffixInexactExhaustive(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                       double error_rate, int minOverlap, Hit& hitTemplate, HitVector* pHits)
{
	// Calculate the number of blocks and the size of each block
	int len = static_cast<int>(w.length());
	int cost = 0;
	for(int i = minOverlap; i <= len; ++i)
	{
		// Calculate the maximum number of differences st diff / i < error_rate
		int diff = int(floor((i * error_rate)));
		cost += alignSuffixMaxDiff(w, pBWT, pRevBWT, diff, i, hitTemplate, pHits);
	}
	return cost;
}


// Seeded blockwise BWT alignment of prefix-suffix for reads
// Each alignment is given a seed region and a block region
// The seed region is the terminal portion of w where maxDiff + 1 seeds are created
// at least 1 of these seeds must align exactly for there to be an alignment with 
// at most maxDiff differences between the prefix/suffix. Only alignments within the
// range [block_start, block_end] are output. The block_end coordinate is inclusive
int _alignBlock(const std::string& w, int block_start, int block_end, const BWT* pBWT, const BWT* pRevBWT, 
                int maxDiff, Hit& hitTemplate, HitVector* pHits)
{
	BWTAlignQueue* pQueue = new BWTAlignQueue;
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
		initIntervals(align, b, pBWT, pRevBWT);		
		pQueue->push(align);
	}

	// Right extension phase
	while(!pQueue->empty())
	{
		++cost;
		BWTAlign align = pQueue->front();
		//std::cout << "Align: "; align.print(w);
		pQueue->pop();

		//printf("Processing: "); align.print(w);
		if(align.dir == ED_RIGHT)
		{
			// Update the interval using RevBWT
			++align.right_index;

			// Flip if we've reached the end of the right extension phase
			// This does not effect the subsequent updates
			if(align.right_index == len - 1)
				align.dir = ED_LEFT;

			// If the length of the alignment is less than the seed length, do not allow mismatches
			if(align.isSeed() || align.z == 0)
			{
				char b = w[align.right_index];
				updateBoth(align, b, pRevBWT);
				if(align.isIntervalValid(RIGHT_INT_IDX))
					pQueue->push(align);
			}
			else
			{
				for(int i = 0; i < 4; ++i)
				{
					char b = ALPHABET[i];
					BWTAlign branched = align;
					updateBoth(branched, b, pRevBWT);

					if(branched.isIntervalValid(RIGHT_INT_IDX))
					{
						if(b != w[align.right_index])
							--branched.z;
						pQueue->push(branched);
					}
				}
			}
		}
		else
		{
			// Left extension

			// If the alignment is within the current block, attempt to output matching prefixes
			if(align.left_index <= block_end)
			{
				int64_t t_lower = pBWT->getC('$') + pBWT->getOcc('$', align.r_lower[LEFT_INT_IDX] - 1);
				int64_t t_upper = pBWT->getC('$') + pBWT->getOcc('$', align.r_upper[LEFT_INT_IDX]) - 1;

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
			if(align.left_index < block_start)
				continue;

			// If there cannot be a branch, only process the matching base
			if(align.z == 0)
			{
				char b = w[align.left_index];
				updateLeft(align, b, pBWT);
				if(align.isIntervalValid(LEFT_INT_IDX))
					pQueue->push(align);
			}
			else
			{
				for(int i = 0; i < 4; ++i)
				{
					char b = ALPHABET[i];
					BWTAlign branched = align;
					// Only update left interval
					updateLeft(branched, b, pBWT);

					if(branched.isIntervalValid(LEFT_INT_IDX))
					{
						if(ALPHABET[i] != w[align.left_index])
							--branched.z;
						pQueue->push(branched);
					}
				}
			}
		}
	}
	delete pQueue;
	return cost;
}

// Seeded BWT alignment of suffix
// This algorithm proceeds in two phases:
// First, for the last minOverlap characters of w
// create maxDiff + 1 seeds. These seeds are extended
// to the right until they hit the end of w. No mismatches
// are allowed in the first ceil(minOverlap / (maxDiff + 1) bases
// Once all the seeds have hit the end of the string, the extension
// flips to going to the left, collecting all the matches prefixes
int alignSuffixMaxDiff(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, int maxDiff, int minOverlap, Hit& hitTemplate, HitVector* pHits)
{
	int cost = 0;
	BWTAlignQueue* pQueue = new BWTAlignQueue;
	// Calculate the initial seeds
	int len = w.length();
	int num_seeds = maxDiff + 1;
	int seed_len = minOverlap / num_seeds;
	int seed_start = len - minOverlap;

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
		//printf("Creating seed at %d to %d, str: %s\n", align.left_index, align.right_index, w.substr(pos, align.seed_len).c_str());

		// Initialize the left and right suffix array intervals
		char b = w[pos];
		initIntervals(align, b, pBWT, pRevBWT);		
		pQueue->push(align);
	}

	// Right extension phase
	while(!pQueue->empty())
	{
		++cost;
		BWTAlign align = pQueue->front();
		//std::cout << "Align: "; align.print(w);
		pQueue->pop();

		//printf("Processing: "); align.print(w);
		if(align.dir == ED_RIGHT)
		{
			// Update the interval using RevBWT
			++align.right_index;

			// Flip if we've reached the end of the right extension phase
			// This does not effect the subsequent updates
			if(align.right_index == len - 1)
				align.dir = ED_LEFT;

			// If the length of the alignment is less than the seed length, do not allow mismatches
			if(align.isSeed() || align.z == 0)
			{
				char b = w[align.right_index];
				updateBoth(align, b, pRevBWT);
				if(align.isIntervalValid(RIGHT_INT_IDX))
					pQueue->push(align);
			}
			else
			{
				for(int i = 0; i < 4; ++i)
				{
					char b = ALPHABET[i];
					BWTAlign branched = align;
					updateBoth(branched, b, pRevBWT);

					if(branched.isIntervalValid(RIGHT_INT_IDX))
					{
						if(b != w[align.right_index])
							--branched.z;
						pQueue->push(branched);
					}
				}
			}
		}
		else
		{
			// Left extension

			// If the overlap is large enough, output the hit
			int overlap_len = len - align.left_index;
			if(overlap_len >= minOverlap)
			{
				int64_t t_lower = pBWT->getC('$') + pBWT->getOcc('$', align.r_lower[LEFT_INT_IDX] - 1);
				int64_t t_upper = pBWT->getC('$') + pBWT->getOcc('$', align.r_upper[LEFT_INT_IDX]) - 1;

				for(int64_t sa_idx = t_lower; sa_idx <= t_upper; ++sa_idx)
				{
					hitTemplate.saIdx = sa_idx;
					hitTemplate.qstart = align.left_index;
					hitTemplate.len = overlap_len;
					hitTemplate.numDiff = maxDiff - align.z;
					//std::cout << "pushing hit of length " << overlap_len << " to saIdx " << sa_idx << "\n";
					pHits->push_back(hitTemplate);
				}
			}

			// Extend hits
			--align.left_index;
			if(align.left_index < 0)
				continue;

			// If there cannot be a branch, only process the matching base
			if(align.z == 0)
			{
				char b = w[align.left_index];
				updateLeft(align, b, pBWT);
				if(align.isIntervalValid(LEFT_INT_IDX))
					pQueue->push(align);
			}
			else
			{
				for(int i = 0; i < 4; ++i)
				{
					char b = ALPHABET[i];
					BWTAlign branched = align;
					// Only update left interval
					updateLeft(branched, b, pBWT);

					if(branched.isIntervalValid(LEFT_INT_IDX))
					{
						if(ALPHABET[i] != w[align.left_index])
							--branched.z;
						pQueue->push(branched);
					}
				}
			}
		}
	}
	delete pQueue;
	return cost;
}

