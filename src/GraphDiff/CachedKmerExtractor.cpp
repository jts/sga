///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// CachedKmerExtractor - Extract reads from an 
// FM-index that share a kmer with a query sequence.
// Previously used kmers are cached to avoid duplicate
// lookups.
//
#include "CachedKmerExtractor.h"
#include "HapgenUtil.h"
#include "Profiler.h"

//
CachedKmerExtractor::CachedKmerExtractor(size_t k, const BWTIndexSet& index_set, 
                                         int minOverlap, double minIdentity) :m_k(k),
                                                                                 m_index_set(index_set),
                                                                                 m_minOverlap(minOverlap),
                                                                                 m_minIdentity(minIdentity)
{

}

// 
SequenceOverlapPairVector CachedKmerExtractor::queryOverlaps(const std::string& query)
{
    PROFILE_FUNC("CachedKmerExtractor::queryOverlaps")
    extractAndUpdate(query);
    extractAndUpdate(reverseComplement(query));

    SequenceOverlapPairVector out;
    // The cache now contains all reads that share a kmer with the query
    // Compute the overlaps
    populateOutput(query, false, &out);
    populateOutput(query, true, &out);
    return out;
}

void CachedKmerExtractor::populateOutput(const std::string& query, bool is_reverse, SequenceOverlapPairVector* out_vector)
{
    PROFILE_FUNC("CachedKmerExtractor::populateOutput")
    size_t nk = query.size() - m_k + 1;
    std::set<size_t> used_indices;

    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = query.substr(i, m_k);
        std::string q_kmer = is_reverse ? reverseComplement(kmer) : kmer;
        const std::vector<size_t>& indices = m_cache_map.find(q_kmer)->second;

        for(size_t j = 0; j < indices.size(); ++j)
        {
            size_t index = indices[j];
            if(used_indices.find(index) != used_indices.end())
                continue;
            used_indices.insert(index);

            std::string match_sequence = m_strings[index];
            if(is_reverse)
                match_sequence = reverseComplement(match_sequence);
        
            // Ignore identical matches
            if(match_sequence.empty() || match_sequence == query)
                continue;

            // Compute the overlap. If the kmer match occurs a single time in each sequence we use
            // the banded extension overlap strategy. Otherwise we use the slow O(M*N) overlapper.
            SequenceOverlap overlap;
            size_t pos_0 = query.find(kmer);
            size_t pos_1 = match_sequence.find(kmer);

            // If there is a single occurrence of the kmer in each read,
            // use that position to seed the overlap calculation
            if(pos_0 == std::string::npos ||
               pos_1 == std::string::npos ||
               query.find(kmer, pos_0 + 1) != std::string::npos || 
               match_sequence.find(kmer, pos_1 + 1) != std::string::npos) 
            {
                // One of the reads has a second occurrence of the kmer. Use
                // the slow overlapper.
                overlap = Overlapper::computeOverlap(query, match_sequence);
            } else {
                overlap = Overlapper::extendMatch(query, match_sequence, pos_0, pos_1, 2);
            }
            bool bPassedOverlap = overlap.getOverlapLength() >= m_minOverlap;
            bool bPassedIdentity = overlap.getPercentIdentity() / 100 >= m_minIdentity;
            if(bPassedOverlap && bPassedIdentity)
            {
                SequenceOverlapPair op;
                op.sequence[0] = query;
                op.sequence[1] = match_sequence;
                op.overlap = overlap;
                op.is_reversed = is_reverse;
                out_vector->push_back(op);
            }
        }
    }

}

// Update the cache such that it contains all the reads containing any k-mer in the query
void CachedKmerExtractor::extractAndUpdate(const std::string& query)
{
    PROFILE_FUNC("CachedKmerExtractor::extractAndUpdate")

    // Check which kmers have not been cached yet
    std::vector<std::string> not_cached_kmers;
    if(query.size() < m_k)
        return;
            
    size_t nk = query.size() - m_k + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = query.substr(i, m_k);
        if(m_cache_map.find(kmer) == m_cache_map.end())
            not_cached_kmers.push_back(kmer);
    }

    printf("Kmers not in cache: %zu of %zu\n", not_cached_kmers.size(), nk);

    // For the kmers not in the cache, lookup the reads using the FM-index
    SeqRecordVector left_vector;
    HapgenUtil::extractHaplotypeReads(not_cached_kmers, m_index_set,
                                      m_k, false, 1000000, &left_vector, NULL);
    
    // Insert the new kmers into the cache. The updateStrings() function below will update
    // the strings vector and the kmer -> index cache
    for(size_t i = 0; i < not_cached_kmers.size(); ++i)
        m_cache_map.insert(std::make_pair(not_cached_kmers[i], std::vector<size_t>()));

    printf("Found %zu new reads\n", left_vector.size());
    // Insert the new strings to the collection.
    for(size_t i = 0; i < left_vector.size(); ++i)
        updateStrings(left_vector[i].seq.toString());
}

//
size_t CachedKmerExtractor::updateStrings(const std::string& incoming)
{
    m_strings.push_back(incoming);
    size_t new_index = m_strings.size() - 1;

    updateCache(incoming, new_index);
    return new_index;
}

// Add the index of a read pair into the cache map.
void CachedKmerExtractor::updateCache(const std::string& str, size_t index)
{
    assert(str.size() >= m_k);
    size_t nk = str.size() - m_k + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = str.substr(i, m_k);
        std::map<std::string, std::vector<size_t> >::iterator iter = m_cache_map.find(kmer);
        if(iter != m_cache_map.end())
            iter->second.push_back(index);
    }
}
