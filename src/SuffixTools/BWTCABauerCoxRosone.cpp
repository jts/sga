//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTCABauerCoxRosone - In-memory version of Illumina's
// BWT construction algorithm
#include "BWTCABauerCoxRosone.h"
#include <algorithm>

void BWTCA::runBauerCoxRosone(const ReadTable* pRT)
{
    size_t num_reads = pRT->getCount();
    size_t num_symbols = pRT->countSumLengths() + num_reads; // include 1 sentinel per read
    printf("Constructing bwt for %zu symbols, %zu reads\n", num_symbols, num_reads);

    // Allocate two working BWTs
    std::string out_bwt(num_symbols, '\0');
    std::string temp_bwt(num_symbols, '\0');

    // Allocate the N array, which stores the indices of the reads, in the order
    // in which they should be inserted for the next round

    // Allocate the P array, which stores the actual position of the inserting symbol
    // in the merged BWT
    NVector n_vector(num_reads);
    
    //
    // Algorithm
    // 

    // Iteration 0:
    // Output a BWT with the last symbol of every read, in the order they appear in the read table.
    // This is the ordering of all suffixes that start with the sentinel character
    outputInitialCycle(pRT, n_vector, out_bwt);

    // Sort the n vector by the symbol of the next cycle (which is the symbol 
    // the next suffix to insert starts with). By using a stable sort here,
    // we implicitly sort suffixes that begin with the same symbol by their
    // previously determined order
    std::stable_sort(n_vector.begin(), n_vector.end());

    // Iterations 1-k:
    // Output 

}

// Update N and the output BWT for the initial cycle, corresponding to the sentinel suffixes
void BWTCA::outputInitialCycle(const ReadTable* pRT, NVector& n_vector, std::string& bwt)
{
    size_t n = pRT->getCount();
    for(size_t i = 0; i < n; ++i)
    {
        size_t rl = pRT->getReadLength(i);
        char c = pRT->getChar(i, rl - 1);
        bwt[i] = c;

        assert(rl > 1);

        // Load the elements of the N vector with the next symbol
        n_vector[i].sym = c;
        n_vector[i].index = i;
    }
}

