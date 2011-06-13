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
    if(0)
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
   ma.print();
   if(consensus.size() > query.size() * 0.8f)
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
std::string LRCorrectionAlgorithm::correctGraphThread(const std::string& query,
                                             const BWT* pTargetBWT, 
                                             const BWT* pRevTargetBWT, 
                                             const SampledSuffixArray* pTargetSSA,
                                             const LRAlignment::LRParams& params)
{
    // Calculate a seed sequence
    size_t ss_length = 150;
    std::string sub = query.substr(0, ss_length);
    LRAlignment::LRHitVector hits;
    LRAlignment::bwaswAlignment(sub, pTargetBWT, pTargetSSA, params, hits);

    if(hits.empty())
        return "";

    size_t bestScore = 0;
    std::string bestString = "";
    for(size_t i = 0; i < hits.size(); ++i)
    {
        LRAlignment::LRHit seedHit = hits[i];
        std::string seed = BWTAlgorithms::extractSubstring(pTargetBWT, seedHit.targetID, seedHit.t_start, seedHit.length);
        if(seed.size() < 51)
            continue;

        std::string querySeed = query.substr(seedHit.q_start);
        StringThreader threader(seed, &querySeed, seedHit.q_end, 51, pTargetBWT, pRevTargetBWT);

        StringVector matchedStrings;
        threader.run(matchedStrings);

        for(size_t j = 0; j < matchedStrings.size(); ++j)
        {
            //StdAlnTools::printGlobalAlignment(query, matchedStrings[j]);
            if(matchedStrings[j].size() > bestScore)
            {
                bestScore = matchedStrings[j].size();
                bestString = matchedStrings[j];
            }
        }
    }
    return bestString;
}
