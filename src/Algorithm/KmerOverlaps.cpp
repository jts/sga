///-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// KmerOverlaps - Overlap computation functions
// seeded by exact kmer matches
//
#include "KmerOverlaps.h"
#include "HashMap.h"
#include "BWTAlgorithms.h"
#include "Profiler.h"

//
MultipleAlignment KmerOverlaps::buildMultipleAlignment(const std::string& query, 
                                                       size_t k,
                                                       int min_overlap,
                                                       double min_identity,
                                                       int bandwidth,
                                                       const BWTIndexSet& indices)
{
    SequenceOverlapPairVector overlap_vector = retrieveMatches(query, k, min_overlap, min_identity, bandwidth, indices);
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", query, "");
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);
    return multiple_alignment;
}

// Struct to hold a partial match in the FM-index
// The position field is the location in the query sequence of this kmer.
// The index field is an index into the BWT. 
// The is_reverse flag indicates the strand of the partial match
struct KmerMatch
{
    int64_t position:16;
    int64_t index:47;
    int64_t is_reverse:1;

    friend bool operator<(const KmerMatch& a, const KmerMatch& b)
    {
        if(a.index == b.index)
            return a.is_reverse < b.is_reverse;
        else
            return a.index < b.index;
    }

    friend bool operator==(const KmerMatch& a, const KmerMatch& b)
    {
        return a.index == b.index && a.is_reverse == b.is_reverse;
    }
};

// Return a hash key for a KmerMatch
struct KmerMatchKey
{
    size_t operator()(const KmerMatch& a) const { return a.index; }
};

typedef std::set<KmerMatch> KmerMatchSet;
typedef HashMap<KmerMatch, bool, KmerMatchKey> KmerMatchMap;

//
SequenceOverlapPairVector KmerOverlaps::retrieveMatches(const std::string& query, size_t k, 
                                                        int min_overlap, double min_identity,
                                                        int bandwidth, const BWTIndexSet& indices)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::retrieveMatches")
    assert(indices.pBWT != NULL);
    assert(indices.pSSA != NULL);

    int64_t max_interval_size = 200;
    SequenceOverlapPairVector overlap_vector;

    // Use the FM-index to look up intervals for each kmer of the read. Each index
    // in the interval is stored individually in the KmerMatchMap. We then
    // backtrack to map these kmer indices to read IDs. As reads can share
    // multiple kmers, we use the map to avoid redundant lookups.
    // There is likely a faster algorithm which performs direct decompression
    // of the read sequences without having to expand the intervals to individual
    // indices. The current algorithm suffices for now.
    KmerMatchMap prematchMap;
    size_t num_kmers = query.size() - k + 1;
    for(size_t i = 0; i < num_kmers; ++i)
    {
        std::string kmer = query.substr(i, k);
        BWTInterval interval = BWTAlgorithms::findInterval(indices, kmer);
        if(interval.isValid() && interval.size() < max_interval_size) 
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
            {
                KmerMatch match = { i, j, false };
                prematchMap.insert(std::make_pair(match, false));
            }
        }

        kmer = reverseComplement(kmer);
        interval = BWTAlgorithms::findInterval(indices, kmer);
        if(interval.isValid() && interval.size() < max_interval_size) 
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
            {
                KmerMatch match = { i, j, true };
                prematchMap.insert(std::make_pair(match, false));
            }
        }
    }

    // Backtrack through the kmer indices to turn them into read indices.
    // This mirrors the calcSA function in SampledSuffixArray except we mark each entry
    // as visited once it is processed.
    KmerMatchSet matches;
    for(KmerMatchMap::iterator iter = prematchMap.begin(); iter != prematchMap.end(); ++iter)
    {
        // This index has been visited
        if(iter->second)
            continue;

        // Mark this as visited
        iter->second = true;

        // Backtrack the index until we hit the starting symbol
        KmerMatch out_match = iter->first;
        while(1) 
        {
            char b = indices.pBWT->getChar(out_match.index);
            out_match.index = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, out_match.index - 1);

            // Check if the hash indicates we have visited this index. If so, stop the backtrack
            KmerMatchMap::iterator find_iter = prematchMap.find(out_match);
            if(find_iter != prematchMap.end())
            {
                // We have processed this index already
                if(find_iter->second)
                    break;
                else
                    find_iter->second = true;
            }

            if(b == '$')
            {
                // We've found the lexicographic index for this read. Turn it into a proper ID
                out_match.index = indices.pSSA->lookupLexoRank(out_match.index);
                matches.insert(out_match);
                break;
            }
        }
    }

    // Refine the matches by computing proper overlaps between the sequences
    // Use the overlaps that meet the thresholds to build a multiple alignment
    for(KmerMatchSet::iterator iter = matches.begin(); iter != matches.end(); ++iter)
    {
        std::string match_sequence = BWTAlgorithms::extractString(indices.pBWT, iter->index);
        if(iter->is_reverse)
            match_sequence = reverseComplement(match_sequence);
        
        // Ignore identical matches
        if(match_sequence == query)
            continue;

        // Compute the overlap. If the kmer match occurs a single time in each sequence we use
        // the banded extension overlap strategy. Otherwise we use the slow O(M*N) overlapper.
        SequenceOverlap overlap;
        std::string match_kmer = query.substr(iter->position, k);
        size_t pos_0 = query.find(match_kmer);
        size_t pos_1 = match_sequence.find(match_kmer);
        assert(pos_0 != std::string::npos && pos_1 != std::string::npos);

        // Check for secondary occurrences
        if(query.find(match_kmer, pos_0 + 1) != std::string::npos || 
           match_sequence.find(match_kmer, pos_1 + 1) != std::string::npos) {
            // One of the reads has a second occurrence of the kmer. Use
            // the slow overlapper.
            overlap = Overlapper::computeOverlap(query, match_sequence);
        } else {
            overlap = Overlapper::extendMatch(query, match_sequence, pos_0, pos_1, bandwidth);
        }

        bool bPassedOverlap = overlap.getOverlapLength() >= min_overlap;
        bool bPassedIdentity = overlap.getPercentIdentity() / 100 >= min_identity;

        if(bPassedOverlap && bPassedIdentity)
        {
            SequenceOverlapPair op;
            op.sequence[0] = query;
            op.sequence[1] = match_sequence;
            op.overlap = overlap;
            op.is_reversed = iter->is_reverse;
            overlap_vector.push_back(op);
        }
    }
    return overlap_vector;
}
