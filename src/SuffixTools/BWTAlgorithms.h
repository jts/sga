//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// bwt_algorithms.h - Algorithms for aligning to a bwt structure
//
#ifndef BWT_ALGORITHMS_H
#define BWT_ALGORITHMS_H

#include "STCommon.h"
#include "BWT.h"
#include "BWTInterval.h"
#include <queue>
#include <list>

#define LEFT_INT_IDX 0
#define RIGHT_INT_IDX 1


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
    AlphaCount l = pRevBWT->getFullOcc(pair.interval[1].lower - 1);
    AlphaCount u = pRevBWT->getFullOcc(pair.interval[1].upper);
    AlphaCount diff = u - l;

    pair.interval[0].lower = pair.interval[0].lower + diff.getLessThan(b);
    pair.interval[0].upper = pair.interval[0].lower + diff.get(b) - 1;

    // Update the right index directly
    size_t pb = pRevBWT->getPC(b);
    pair.interval[1].lower = pb + l.get(b);
    pair.interval[1].upper = pb + u.get(b) - 1;
}

//
// Update both the left and right intervals using pBWT
// This assumes that the left/right ranges in pair are for string S
// It returns the updated left/right ranges for string bS (prepending b)
inline void updateBothL(BWTIntervalPair& pair, char b, const BWT* pBWT)
{
    // Update the left index using the difference between the AlphaCounts in the reverse table
    AlphaCount l = pBWT->getFullOcc(pair.interval[0].lower - 1);
    AlphaCount u = pBWT->getFullOcc(pair.interval[0].upper);
    AlphaCount diff = u - l;

    pair.interval[1].lower = pair.interval[1].lower + diff.getLessThan(b);
    pair.interval[1].upper = pair.interval[1].lower + diff.get(b) - 1;

    // Update the left index directly
    size_t pb = pBWT->getPC(b);
    pair.interval[0].lower = pb + l.get(b);
    pair.interval[0].upper = pb + u.get(b) - 1;
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

// Return the count of all the possible one base extensions of the string w.
// This returns the number of times the suffix w[i, l]A, w[i, l]C, etc 
// appears in the FM-index for all i s.t. length(w[i, l]) >= minOverlap.
AlphaCount calculateExactExtensions(const unsigned int overlapLen, const std::string& w, const BWT* pBWT, const BWT* pRevBWT);

};

#endif
