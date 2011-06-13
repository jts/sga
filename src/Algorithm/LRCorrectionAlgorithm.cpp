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

std::string LRCorrectionAlgorithm::correct(const std::string& query,
                                           const BWT* pTargetBWT, 
                                           const BWT* pRevTargetBWT, 
                                           const SampledSuffixArray* pTargetSSA,
                                           const LRAlignment::LRParams& params)
{
    if(1)
        return correctAlignment(query, pTargetBWT, pTargetSSA, params);
    else
        return correctGraphThread(query, pTargetBWT, pRevTargetBWT, pTargetSSA, params);
}

// Correct the read using an alignment-based approach
std::string LRCorrectionAlgorithm::correctAlignment(const std::string& query,
                                           const BWT* pTargetBWT, 
                                           const SampledSuffixArray* pTargetSSA,
                                           const LRAlignment::LRParams& params)
{
   LRAlignment::LRHitVector hits;
   LRAlignment::bwaswAlignment(query, pTargetBWT, pTargetSSA, params, hits);
   MultiAlignment ma = LRAlignment::convertHitsToMultiAlignment(query, pTargetBWT, pTargetSSA, params, hits);
   std::string consensus = ma.generateConsensus();
   if(consensus.size() > query.size() * 0.9f)
       return consensus;
    else
        return "";
}


// Correct the read using an alignment-based approach
std::string LRCorrectionAlgorithm::correctAlignmentPartial(const std::string& query,
                                                  const BWT* pTargetBWT, 
                                                  const SampledSuffixArray* pTargetSSA,
                                                  const LRAlignment::LRParams& params)
{

    size_t ss_length = 125;
    size_t stride = 50;
    for(size_t i = 0; i < query.size() - ss_length; i += stride)
    {
        if(query.size() < ss_length)
            break;
        std::string sub = query.substr(i, ss_length);
    
 //       std::string sub = query;
        LRAlignment::LRHitVector hits;
        LRAlignment::bwaswAlignment(sub, pTargetBWT, pTargetSSA, params, hits);
        std::cout << "hits: " << hits.size() << "\n";
        MultiAlignment ma = LRAlignment::convertHitsToMultiAlignment(sub, pTargetBWT, pTargetSSA, params, hits);
        std::string consensus = ma.generateConsensus();
    }

    return query;
}

// Attempt to correct the sequence via threading through a graph
std::string LRCorrectionAlgorithm::correctGraphThread(const std::string& query,
                                             const BWT* pTargetBWT, 
                                             const BWT* pRevTargetBWT, 
                                             const SampledSuffixArray* pTargetSSA,
                                             const LRAlignment::LRParams& params)
{
    // Calculate a seed sequence
    size_t ss_length = 200;
    std::string sub = query.substr(0, ss_length);
    LRAlignment::LRHitVector hits;
    LRAlignment::bwaswAlignment(sub, pTargetBWT, pTargetSSA, params, hits);

    if(hits.empty())
        return "";

    LRAlignment::LRHit seedHit = hits[0];
    std::string seed = BWTAlgorithms::extractSubstring(pTargetBWT, seedHit.targetID, seedHit.t_start, seedHit.length);
    if(seed.size() < 51)
        return "";

    std::string querySeed = query.substr(seedHit.q_start);
    StringThreader threader(seed, &querySeed, seedHit.q_end, 51, pTargetBWT, pRevTargetBWT);

    StringVector matchedStrings;
    threader.run(matchedStrings);

    if(matchedStrings.empty())
        return "";

    for(size_t i = 0; i < matchedStrings.size(); ++i)
    {
        StdAlnTools::printGlobalAlignment(query, matchedStrings[i]);
    }
    return matchedStrings.back();
}
