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
#ifndef CACHED_KMER_EXTRACTOR_H
#define CACHED_KMER_EXTRACTOR_H

#include <string>
#include <map>
#include "overlapper.h"
#include "BWTIndexSet.h"
#include "KmerOverlaps.h"

class CachedKmerExtractor
{
    public:
        CachedKmerExtractor(size_t k, const BWTIndexSet& index_set, int minOverlap, double minIdentity);
        
        // Find read pairs that share a k-mer with one of the query strings
        SequenceOverlapPairVector queryOverlaps(const std::string& query);

    private:

        //
        void populateOutput(const std::string& query, bool is_reverse, SequenceOverlapPairVector* out_vector);

        // Extract reads with k-mers matching the query from the FM-index
        // and update the cache
        void extractAndUpdate(const std::string& query);
        
        // Update the collection of cached sequences
        size_t updateStrings(const std::string& incoming);

        // Update the cache map
        void updateCache(const std::string& str, size_t index);
        
        // Data
        size_t m_k;
        BWTIndexSet m_index_set;
        int m_minOverlap;
        double m_minIdentity;

        // The first and second reads in a pair
        std::vector<std::string> m_strings;
        std::map<std::string, std::vector<size_t> > m_cache_map;
};

#endif

