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
#include "BWTInterval.h"
#include "OverlapData.h"
#include <queue>
#include <list>

#define LEFT_INT_IDX 0
#define RIGHT_INT_IDX 1

// types
enum ExtendDirection
{
	ED_LEFT,
	ED_RIGHT
};

// Structure holding all the working variables for making an inexact alignment for a sequence
// to a BWT
struct BWTAlign
{
	inline int length() const { return right_index - left_index + 1; }
	inline bool isSeed() const { return length() < seed_len; }
	inline bool isIntervalValid(int idx) { return ranges.interval[idx].isValid(); }

	int left_index; // inclusive
	int right_index; // inclusive
	int seed_len;
	ExtendDirection dir; // the direction that this alignment is being extended in
	int z;

	BWTIntervalPair ranges; // ranges.interval[0] is the left interval, 1 is the right 
	
	// Sort the alignments based on their r_lower/r_upper
	static inline bool compareLeftRange(const BWTAlign& a, const BWTAlign& b)
	{
		return BWTInterval::compare(a.ranges.interval[0], b.ranges.interval[0]);
	}

	// Compare for equality based on the left range
	// If the length of the alignment is equal, then if the left ranges
	// are a match, the two alignment objects are redundant and one can be removed
	static inline bool equalLeftRange(const BWTAlign& a, const BWTAlign& b)
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

typedef std::queue<BWTAlign> BWTAlignQueue;
typedef std::list<BWTAlign> BWTAlignList;
typedef std::list<OverlapBlock> OverlapBlockList;
typedef OverlapBlockList::iterator OBLIter;

// functions
namespace BWTAlgorithms
{

// get the interval(s) in pBWT/pRevBWT that corresponds to the string w using a backward search algorithm
BWTInterval findInterval(const BWT* pBWT, const std::string& w);
BWTIntervalPair findIntervalPair(const BWT* pBWT, const BWT* pRevBWT, const std::string& w);

// Update the given interval using backwards search
// If the interval corrsponds to string S, it will be updated 
// for string bS
inline void updateInterval(BWTInterval& interval, char b, const BWT* pBWT)
{
	size_t pb = pBWT->getPC(b);
	interval.lower = pb + pBWT->getOcc(b, interval.lower - 1);
	interval.upper = pb + pBWT->getOcc(b, interval.upper) - 1;
}

//
// Update both the left and right intervals using pRevBWT
// This assumes that the left/right ranges in pair are for string S
// It returns the updated left/right ranges for string Sb (appending b)
// using the pRevBWT to update both
inline void updateBothR(BWTIntervalPair& pair, char b, const BWT* pRevBWT)
{
	// Update the left index using the difference between the AlphaCounts in the reverse table
	AlphaCount diff = pRevBWT->getOccDiff(pair.interval[1].lower - 1, pair.interval[1].upper);
	pair.interval[0].lower = pair.interval[0].lower + diff.getLessThan(b);
	pair.interval[0].upper = pair.interval[0].lower + diff.get(b) - 1;

	// Update the right index directly
	updateInterval(pair.interval[1], b, pRevBWT);
}

//
// Update both the left and right intervals using pBWT
// This assumes that the left/right ranges in pair are for string S
// It returns the updated left/right ranges for string bS (prepending b)
inline void updateBothL(BWTIntervalPair& pair, char b, const BWT* pBWT)
{
	// Update the left index using the difference between the AlphaCounts in the reverse table
	AlphaCount diff = pBWT->getOccDiff(pair.interval[0].lower - 1, pair.interval[0].upper);
	pair.interval[1].lower = pair.interval[1].lower + diff.getLessThan(b);
	pair.interval[1].upper = pair.interval[1].lower + diff.get(b) - 1;

	// Update the left index directly
	updateInterval(pair.interval[0], b, pBWT);

}

// Initialize the interval of index idx to be the range containining all the b suffixes
inline void initInterval(BWTInterval& interval, char b, const BWT* pB)
{
	interval.lower = pB->getPC(b);
	interval.upper = interval.lower + pB->getOcc(b, pB->getBWLen() - 1) - 1;
}

// Initialize the interval of index idx to be the range containining all the b suffixes
inline void initIntervalPair(BWTIntervalPair& pair, char b, const BWT* pBWT, const BWT* pRevBWT)
{
	initInterval(pair.interval[LEFT_INT_IDX], b, pBWT);
	initInterval(pair.interval[RIGHT_INT_IDX], b, pRevBWT);
}

// Return the counts of the bases between the lower and upper interval in pBWT
inline AlphaCount getExtCount(const BWTInterval& interval, const BWT* pBWT)
{
	return pBWT->getOccDiff(interval.lower - 1, interval.upper);
}


//
// Exact alignment algorithms
//


// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
// overlaps with a suffix of w. The ranges are added to the pOBList. If w has a full-length
// overlap with some other string, the containment overlap will be placed in pHits for 
// processing later
void findOverlapBlocks(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                       int minOverlap, const AlignFlags& af, OverlapBlockList* pOBTemp, OverlapBlockList* pOBFinal);


// Using the vector of OverlapBlocks, calculate the irreducible hits and output them to pHits
void calculateIrreducibleHits(OverlapBlockList* pOBList, OverlapBlockList* pOBFinal);

// Extend each block in obl until all the irreducible overlaps have been found. 
void processIrreducibleBlocks(OverlapBlockList& obList, OverlapBlockList* pOBFinal);

// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
void updateOverlapBlockRangesRight(OverlapBlockList& obList, char b);


// Return the count of all the possible one base extensions of the string w.
// This returns the number of times the suffix w[i, l]A, w[i, l]C, etc 
// appears in the FM-index for all i s.t. length(w[i, l]) >= minOverlap.
AlphaCount calculateExactExtensions(const unsigned int overlapLen, const std::string& w, const BWT* pBWT, const BWT* pRevBWT);

// Perform an exact suffix overlap
int alignSuffixExact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                     int minOverlap, Hit& hitTemplate, HitVector* pHits);


//
// Inexact alignment algorithms
//

// Perform an inexact suffix overlap, allowing at most error_rate errors
int alignSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                       double error_rate, int minOverlap, Hit& hitTemplate, HitVector* pHits);


// Align a subrange of a string against and fm-index while allowing maxDiff errors. Seeded.
int _alignBlock(const std::string& w, int block_start, int block_end,
                const BWT* pBWT, const BWT* pRevBWT, int maxDiff, Hit& hitTemplate, HitVector* pHits);


};

#endif
