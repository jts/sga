//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapSeed.h - Data structure holding a partial 
// alignment to a FM-index. 
//
#ifndef OVERLAPSEED_H
#define OVERLAPSEED_H

#include <queue>
#include <list>

// types
enum ExtendDirection
{
	ED_LEFT,
	ED_RIGHT
};

// Structure holding all the working variables for making an inexact alignment for a sequence
// to a BWT
struct OverlapSeed
{
	// Functions
	inline int length() const { return right_index - left_index + 1; }
	inline bool isSeed() const { return length() < seed_len; }
	inline bool isIntervalValid(int idx) { return ranges.interval[idx].isValid(); }
	inline bool allowMismatches() const { return z < maxDiff; }
	inline double calcErrorRate() const { return static_cast<double>(z) / static_cast<double>(length()); }
	 
	// Data
	
	// Index range is inclusive on both ends
	int left_index;
	int right_index;
	int seed_len;

	// The direction the alignment is currently being extended
	ExtendDirection dir; 
	
	// The current number of mismatches in the block
	int z;

	// The total number of allowed mismatches
	int maxDiff;

	// BWT interval coordinates, the first element of the pair is the left range
	BWTIntervalPair ranges; 
	
	// Sort the alignments based on their r_lower/r_upper
	static inline bool compareLeftRange(const OverlapSeed& a, const OverlapSeed& b)
	{
		return BWTInterval::compare(a.ranges.interval[0], b.ranges.interval[0]);
	}

	// Compare for equality based on the left range
	// If the length of the alignment is equal, then if the left ranges
	// are a match, the two alignment objects are redundant and one can be removed
	static inline bool equalLeftRange(const OverlapSeed& a, const OverlapSeed& b)
	{
#ifdef VALIDATE
		if(BWTInterval::equal(a.ranges.inteval[0], b.ranges.interval[0]))
			assert(a.length() == b.length() && a.z == b.z);
#endif
		return BWTInterval::equal(a.ranges.interval[0], b.ranges.interval[0]);
	}
	

	void print() const
	{
		printf("li: %d ri: %d sl: %d dir: %d z: %d lrl: %d lru: %d rlr: %d rlu: %d\n", 
		        left_index, right_index, seed_len, dir, z, 
				(int)ranges.interval[0].lower, (int)ranges.interval[0].upper, 
				(int)ranges.interval[1].lower, (int)ranges.interval[1].upper);
	}


	void print(const std::string& w) const
	{
		printf("sub: %s li: %d ri: %d sl: %d dir: %d z: %d lrl: %d lru: %d rlr: %d rlu: %d\n", 
		        w.substr(left_index, length()).c_str(), left_index, right_index,
				seed_len, dir, z, (int)ranges.interval[0].lower, (int)ranges.interval[0].upper, 
				(int)ranges.interval[1].lower, (int)ranges.interval[1].upper);
	}
};

// Collections
typedef std::queue<OverlapSeed> OverlapSeedQueue;
typedef std::list<OverlapSeed> OverlapSeedList;

#endif
