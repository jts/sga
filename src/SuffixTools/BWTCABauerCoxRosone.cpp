//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTCABauerCoxRosone - In-memory version of Illumina's
// BWT construction algorithm
#include "BWTCABauerCoxRosone.h"
#include "Timer.h"
#include <algorithm>
#include <iterator>

void BWTCA::runBauerCoxRosone(const ReadTable* pRT)
{
    size_t num_reads = pRT->getCount();
    size_t num_symbols = pRT->countSumLengths() + num_reads; // include 1 sentinel per read
    printf("Constructing bwt for %zu symbols, %zu reads\n", num_symbols, num_reads);

    // Allocate two working BWTs
    BWTString out_bwt;
    BWTString temp_bwt;

    out_bwt.resize(num_symbols);
    temp_bwt.resize(num_symbols);

    std::cout << "nelem size: " << sizeof(BCRElem) << "\n";

    // Allocate the N array, which stores the indices of the reads, in the order
    // in which they should be inserted for the next round


    BCRVector bcrVector(num_reads);
    
    //
    // Algorithm
    // 

    // Track the number of suffixes that start with a given symbol
    AlphaCount64 suffixStartCounts;

    // Iteration 0:
    // Output a BWT with the last symbol of every read, in the order they appear in the read table.
    // This is the ordering of all suffixes that start with the sentinel character

    outputInitialCycle(pRT, bcrVector, out_bwt, suffixStartCounts);

    size_t maxCycles = pRT->getReadLength(0);
    size_t partial_bwt_length = pRT->getCount();
    Timer timer("cycles");
    for(size_t cycle = 2; cycle <= maxCycles; ++cycle)
    {
        std::cout << "Starting cycle " << cycle << " " << timer.getElapsedWallTime() << "\n";
        
        // Convert relative positions to absolute positions
        calculateAbsolutePositions(bcrVector, suffixStartCounts);
       
        // Sort nvector by the absolute position to insert the symbol
        std::sort(bcrVector.begin(), bcrVector.end());
        std::cout << "  done sorting..." << timer.getElapsedWallTime() << "\n";

        /*
        std::cout << "BCRVector absolute:\n";
        std::copy(bcrVector.begin(), bcrVector.end(), std::ostream_iterator<NElem>(std::cout, "\n"));
        */

        // Output the BWT for this cycle and update the vector and suffix start count
        partial_bwt_length = outputPartialCycle(cycle, pRT, bcrVector, out_bwt, partial_bwt_length, temp_bwt, suffixStartCounts);
        std::cout << "  done writing..." << timer.getElapsedWallTime() << "\n";

        // Swap the in/out bwt
        out_bwt.swap(temp_bwt);
    }
}

// Write out the next BWT for the next cycle. This updates BCRVector
// and suffixSymbolCounts. Returns the number of symbols written to writeBWT
size_t BWTCA::outputPartialCycle(int cycle,
                               const ReadTable* pRT,
                               BCRVector& bcrVector, 
                               const BWTString& readBWT, 
                               size_t total_read_symbols,
                               BWTString& writeBWT, 
                               AlphaCount64& suffixStartCounts)
{
    // We track the rank of each symbol as it is copied/inserted
    // into the new bwt
    AlphaCount64 rank;

    // Counters
    size_t num_copied = 0;
    size_t num_inserted = 0;
    size_t num_wrote = 0;

    for(size_t i = 0; i < bcrVector.size(); ++i)
    {
        BCRElem& ne = bcrVector[i];
        
        // Copy elements from the read bwt until we reach the target position
        while(num_copied + num_inserted < ne.position)
        {
            char c = readBWT.get(num_copied++);
            writeBWT.set(num_wrote++, c);
            rank.increment(c);
        }

        // Now insert the incoming symbol
        size_t rl = pRT->getReadLength(ne.index);
        char c = pRT->getChar(ne.index, rl - cycle);
        //std::cout << "Inserting " << c << " at position " << num_copied + num_inserted << "\n";
        writeBWT.set(num_wrote++, c);
        num_inserted += 1;

        // Update the nvector element
        ne.sym = c;

        // Record the rank of the inserted symbol
        ne.position = rank.get(c);

        // Update the rank and the number of suffixes that start with c
        rank.increment(c);
        suffixStartCounts.increment(c);
    }

    // Copy any remaining symbols in the bwt
    while(num_copied < total_read_symbols)
        writeBWT.set(num_wrote++, readBWT.get(num_copied++));

//    printf("Wrote: %zu Copied: %zu Inserted: %zu\n", num_wrote, num_copied, num_inserted);
//    std::cout << "new bwt: " << writeBWT.substr(0, num_wrote) << "\n";
    
    return num_wrote;
}

// Update N and the output BWT for the initial cycle, corresponding to the sentinel suffixes
// the symbolCounts vector is updated to hold the number of times each symbol has been inserted
// into the bwt
void BWTCA::outputInitialCycle(const ReadTable* pRT, BCRVector& bcrVector, BWTString& bwt, AlphaCount64& suffixSymbolCounts)
{
    AlphaCount64 incomingSymbolCounts;

    size_t n = pRT->getCount();
    for(size_t i = 0; i < n; ++i)
    {
        size_t rl = pRT->getReadLength(i);
        char c = pRT->getChar(i, rl - 1);
        bwt.set(i, c);

        assert(rl > 1);

        // Load the elements of the N vector with the next symbol
        bcrVector[i].sym = c;
        bcrVector[i].index = i;

        // Set the relative position of the symbol that is being inserted
        bcrVector[i].position = incomingSymbolCounts.get(c);

        // Increment the count of the first base of the suffix of the
        // incoming strings. This is $ for the initial cycle
        suffixSymbolCounts.increment('$');

        // Update the inserted symbols
        incomingSymbolCounts.increment(c);
    }

    suffixSymbolCounts += incomingSymbolCounts;
}

void BWTCA::calculateAbsolutePositions(BCRVector& bcrVector, const AlphaCount64& suffixSymbolCounts)
{
    // Calculate a predecessor array from the suffix symbol counts
    AlphaCount64 predCounts;
    for(int i = 0; i < BWT_ALPHABET::size; ++i)
    {
        char b = RANK_ALPHABET[i];
        int64_t pc = suffixSymbolCounts.getLessThan(b);
        predCounts.set(b, pc);
    }

    for(size_t i = 0; i < bcrVector.size(); ++i)
        bcrVector[i].position += predCounts.get(bcrVector[i].sym);
}

