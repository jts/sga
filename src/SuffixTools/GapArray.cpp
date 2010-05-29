//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapArray - Data structure and functions used to count 
// the number of times a suffix of a given rank occurs in a data set
//
#include "GapArray.h"

// Increment the gap array for each suffix of seq
void updateGapArray(const DNAString& w, const BWT* pBWTInternal, GapArray& gap_array)
{
    size_t l = w.length();
    int i = l - 1;

    // Compute the rank of the last character of seq. We consider $ to be implicitly
    // terminated by a $ character. The first rank calculated is for this and it is given
    // by the C(a) array in BWTInternal
    int64_t rank = pBWTInternal->getPC('$'); // always zero
    incrementGapArray(rank, gap_array);

    // Compute the starting rank for the last symbol of w
    char c = w.get(i);
    rank = pBWTInternal->getPC(c);
    incrementGapArray(rank, gap_array);
    --i;

    // Iteratively compute the remaining ranks
    while(i >= 0)
    {
        char c = w.get(i);
        rank = pBWTInternal->getPC(c) + pBWTInternal->getOcc(c, rank - 1);
        //std::cout << "c: " << c << " rank: " << rank << "\n";
        incrementGapArray(rank, gap_array);
        --i;
    }
}

//
void incrementGapArray(int64_t rank, GapArray& gap_array)
{
    static size_t max_gap_count = std::numeric_limits<GAP_TYPE>::max();
    assert(gap_array[rank] < max_gap_count);
    assert(rank < (int64_t)gap_array.size());
    ++gap_array[rank];
}
