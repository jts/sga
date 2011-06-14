//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRCorrectionAlgorithmAlgorithm - Collection of algorithms for 
// correcting long reads with a set of short reads
//
#include "LRCorrectionAlgorithm.h"
#include "StringThreader.h"

namespace LRCorrectionAlgorithm
{

std::string correct(const std::string& query,
                    const BWT* pTargetBWT, 
                    const BWT* pRevTargetBWT, 
                    const SampledSuffixArray* pTargetSSA,
                    const LRAlignment::LRParams& params)
{
    if(0)
        return correctAlignment(query, pTargetBWT, pTargetSSA, params);
    else
        return correctGraphThread(query, pTargetBWT, pRevTargetBWT, pTargetSSA, params);
}

// Correct the read using an alignment-based approach
std::string correctAlignment(const std::string& query,
                             const BWT* pTargetBWT, 
                             const SampledSuffixArray* pTargetSSA,
                             const LRAlignment::LRParams& params)
{
   LRAlignment::LRHitVector hits;
   LRAlignment::bwaswAlignment(query, pTargetBWT, pTargetSSA, params, hits);
   MultiAlignment ma = LRAlignment::convertHitsToMultiAlignment(query, pTargetBWT, pTargetSSA, params, hits);
   std::string consensus = ma.generateConsensus();
   ma.print();
   if(consensus.size() > query.size() * 0.8f)
       return consensus;
    else
        return "";
}


// Correct the read using an alignment-based approach
std::string correctAlignmentPartial(const std::string& query,
                                    const BWT* pTargetBWT, 
                                    const SampledSuffixArray* pTargetSSA,
                                    const LRAlignment::LRParams& params)
{

    size_t ss_length = 150;
    size_t stride = 50;
    LRAlignment::LRHitVector allHits;

    for(size_t i = 0; i < query.size() - ss_length; i += stride)
    {
        if(query.size() < ss_length)
            break;
        std::string sub = query.substr(i, ss_length);
    
        LRAlignment::LRHitVector hits;
        LRAlignment::bwaswAlignment(sub, pTargetBWT, pTargetSSA, params, hits);

        for(size_t j = 0; j < hits.size(); ++j)
        {
            LRAlignment::LRHit shiftHit = hits[j];
            shiftHit.q_start += i;
            shiftHit.q_end += i;
            allHits.push_back(shiftHit);
        }
    }

    MultiAlignment ma = LRAlignment::convertHitsToMultiAlignment(query, pTargetBWT, pTargetSSA, params, allHits);
    ma.print();
    std::string consensus = ma.generateConsensus();
    if(consensus.size() > query.size() * 0.9f)
       return consensus;
    else
        return "";
}

// Attempt to correct the sequence via threading through a graph
std::string correctGraphThread(const std::string& query,
                               const BWT* pTargetBWT, 
                               const BWT* pRevTargetBWT, 
                               const SampledSuffixArray* pTargetSSA,
                               const LRAlignment::LRParams& params)
{
    // Calculate hits using bwa-sw alignment over the first portion
    // of the long read. The hits will be used to seed the threading
    // process.
    size_t extension_kmer = 51;
    size_t ss_length = 150;
    std::string sub = query.substr(0, ss_length);
    LRAlignment::LRHitVector hits;
    LRAlignment::bwaswAlignment(sub, pTargetBWT, pTargetSSA, params, hits);

    // No hits found
    if(hits.empty())
        return "";

    int bestScore = std::numeric_limits<int>::min();
    std::string bestString = "";

    for(size_t i = 0; i < hits.size(); ++i)
    {
        LRAlignment::LRHit seedHit = hits[i];
        std::string seed = BWTAlgorithms::extractSubstring(pTargetBWT, seedHit.targetID, seedHit.t_start, seedHit.length);
        if(seed.size() < extension_kmer)
            continue;

        // Perform the threading for this query seed
        std::string querySeed = query.substr(seedHit.q_start);
        StringThreader threader(seed, &querySeed, seedHit.q_end, extension_kmer, pTargetBWT, pRevTargetBWT);

        StringThreaderResultVector results;
        threader.run(results);
        
        // Process the results and save the best-scoring thread
        for(size_t j = 0; j < results.size(); ++j)
        {
            StringThreaderResult& r = results[j];
            int score = StdAlnTools::globalAlignment(querySeed.substr(0, r.query_align_length), r.thread);
            
            if(score > bestScore)
            {
                bestScore = score;
                bestString = r.thread;
            }
        }
    }
    return bestString;
}

} // namespace
