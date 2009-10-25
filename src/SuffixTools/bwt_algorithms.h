//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// bwt_algorithms.h - Algorithms for aligning to a bwt structure
//
#ifndef BWT_ALGORITHMS_H
#define BWT_ALGORITHMS_H

#include "STCommon.h"
#include "BWT.h"
#include <queue>

// types
enum ExtendDirection
{
	ED_LEFT,
	ED_RIGHT
};

struct BWTAlign
{
	inline int length() const { return right_index - left_index + 1; }
	inline bool isSeed() const { return length() <= seed_len; }
	inline bool isIntervalValid(int idx) { return r_lower[idx] <= r_upper[idx]; }

	int left_index; // inclusive
	int right_index; // inclusive
	int seed_len;
	ExtendDirection dir; // the direction that this alignment is being extended in
	int z;

	int64_t r_lower[2]; // idx 0  is the left interval, 1 is the right 
	int64_t r_upper[2]; // as above

	void print() const
	{
		printf("li: %d ri: %d sl: %d dir: %d z: %d lrl: %d lru: %d rlr: %d rlu: %d\n", 
		        left_index, right_index, seed_len, dir, z, (int)r_lower[0], (int)r_upper[0], (int)r_lower[1], (int)r_upper[1]);
	}


	void print(const std::string& w) const
	{
		printf("sub: %s li: %d ri: %d sl: %d dir: %d z: %d lrl: %d lru: %d rlr: %d rlu: %d\n", 
		        w.substr(left_index, length()).c_str(), left_index, right_index,
				seed_len, dir, z, (int)r_lower[0], (int)r_upper[0], (int)r_lower[1], (int)r_upper[1]);
	}
};

typedef std::queue<BWTAlign> BWTAlignQueue;

// functions
int alignSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                       double error_rate, int minOverlap, Hit& hitTemplate, HitVector* pHits);

int alignSuffixInexactExhaustive(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                       double error_rate, int minOverlap, Hit& hitTemplate, HitVector* pHits);

int alignSuffixMaxDiff(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, int maxDiff, int minOverlap, Hit& hitTemplate, HitVector* pHits);

int _alignBlock(const std::string& w, int block_start, int block_end,
                const BWT* pBWT, const BWT* pRevBWT, int maxDiff, Hit& hitTemplate, HitVector* pHits);




#endif
