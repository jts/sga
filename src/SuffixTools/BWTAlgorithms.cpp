//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// bwt_algorithms.cpp - Algorithms for aligning to a bwt structure
//
#include "BWTAlgorithms.h"

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

// Count the number of occurrences of string w, including the reverse complement
size_t BWTAlgorithms::countSequenceOccurrences(const std::string& w, const BWT* pBWT)
{
    BWTInterval fwd_interval = findInterval(pBWT, w);
    BWTInterval rc_interval = findInterval(pBWT, reverseComplement(w));

    size_t count = 0;
    if(fwd_interval.isValid())
        count += fwd_interval.size();
    if(rc_interval.isValid())
        count += rc_interval.size();
    return count;
}

// Return the count of all the possible one base extensions of the string w.
// This returns the number of times the suffix w[i, l]A, w[i, l]C, etc 
// appears in the FM-index for all i s.t. length(w[i, l]) == overlapLen.
AlphaCount64 BWTAlgorithms::calculateExactExtensions(const unsigned int overlapLen, const std::string& w, const BWT* pBWT, const BWT* pRevBWT)
{
    // The algorithm is as follows:
    // We perform a backward search on the FM-index of w.
    // For each signficant suffix (length w[i,l] >= minOverlap)
    // we determine the proper prefixes that match w[i,l]. For each proper prefix matching, 
    // we compute the number of extensions of A,C,G,T for those prefix.
    AlphaCount64 ext_counts;
    BWTIntervalPair ranges;
    size_t l = w.length();
    int start = l - 1;
    BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);

    for(int i = start - 1; i >= 0; --i)
    {
        // Compute the range of the suffix w[i, l]
        BWTAlgorithms::updateBothL(ranges, w[i], pBWT);

        // Break if the suffix is no longer found
        if(!(ranges.interval[0].isValid() && ranges.interval[1].isValid())) 
            break;

        if((l - i) == overlapLen)
        {
            if(ranges.interval[1].isValid())
            {
                assert(ranges.interval[1].lower > 0);
                // The count for each extension is the difference between rank(B, upper) and rank(B, lower - 1)
                AlphaCount64 ac = pRevBWT->getOccDiff(ranges.interval[1].lower - 1, ranges.interval[1].upper);
                ext_counts += ac;
            }
        }
    }
    return ext_counts;
}

// 
std::string BWTAlgorithms::sampleRandomString(const BWT* pBWT)
{
    assert(RAND_MAX > 0x7FFF);
    size_t n = pBWT->getNumStrings();
    size_t idx = rand() % n;

    std::string out;

    // The range [0,n) in the BWT contains all the terminal
    // symbols for the reads. Search backwards from one of them
    // until the '$' is found gives a full string.

    BWTInterval interval(idx, idx);
    while(1)
    {
        assert(interval.isValid());
        char b = pBWT->getChar(interval.lower);
        if(b == '$')
            break;
        else
            out.push_back(b);
        updateInterval(interval, b, pBWT);
    } 
    return reverse(out);
}
