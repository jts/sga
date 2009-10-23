//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// bwt_algorithms.cpp - Algorithms for aligning to a bwt structure
//
#include "bwt_algorithms.h"

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


// Seeded BWT alignment of suffix
// This algorithm proceeds in two phases:
// First, for the last minOverlap characters of w
// create maxDiff + 1 seeds. These seeds are extended
// to the right until they hit the end of w. No mismatches
// are allowed in the first ceil(minOverlap / (maxDiff + 1) bases
// Once all the seeds have hit the end of the string, the extension
// flips to going to the left, collecting all the matches prefixes
int alignInexactSuffix(std::string w, const BWT* pBWT, const BWT* pRevBWT, int maxDiff, int minOverlap, Hit& hitTemplate, HitVector* pHits)
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
					std::cout << "pushing hit of length " << overlap_len << " to saIdx " << sa_idx << "\n";
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

