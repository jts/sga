//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRCorrection - Collection of algorithms for 
// correcting long reads with a set of short reads
//
#include "LRCorrection.h"
#include "StringBuilder.h"

std::string LRCorrection::correctLongRead(const std::string& query,
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
std::string LRCorrection::correctAlignment(const std::string& query,
                                           const BWT* pTargetBWT, 
                                           const SampledSuffixArray* pTargetSSA,
                                           const LRAlignment::LRParams& params)
{
    LRAlignment::LRHitVector hits;
    LRAlignment::bwaswAlignment(query, pTargetBWT, pTargetSSA, params, hits);
    MultiAlignment ma = LRAlignment::convertHitsToMultiAlignment(query, pTargetBWT, pTargetSSA, params, hits);
    std::string consensus = ma.generateConsensus();
    
    addOverlappingHits(query, pTargetBWT, pTargetSSA, params, hits);

    return query;
}

// Attempt to correct the sequence via threading through a graph
std::string LRCorrection::correctGraphThread(const std::string& query,
                                             const BWT* pTargetBWT, 
                                             const BWT* pRevTargetBWT, 
                                             const SampledSuffixArray* /*pTargetSSA*/,
                                             const LRAlignment::LRParams& /*params*/)
{
    int kmer = 61;
    int seedPosition = -1;
    std::string seed = findSeedStringNaive(query, pTargetBWT, kmer, seedPosition);
    std::string querySub = query.substr(seedPosition);

    if(!seed.empty())
    {
        StringBuilder stringBuilder(seed, kmer, pTargetBWT, pRevTargetBWT);
        int extendCount = query.size() - kmer;

        int cullDistance = 20;
        int distanceFromLastCull = 0;
        while(extendCount--)
        {
            stringBuilder.extendOnce();
            distanceFromLastCull += 1;
            if(distanceFromLastCull == cullDistance)
            {
                stringBuilder.cull(querySub, 20);
                distanceFromLastCull = 0;
            }
        }
        stringBuilder.cull(querySub, 2);
        //stringBuilder.print(querySub);
    }
    else
    {
        std::cout << "No seed string found\n";
    }
    return query;
}

//
std::string LRCorrection::findSeedStringNaive(const std::string& query,
                                              const BWT* pTargetBWT, 
                                              int k,
                                              int& position)
{
    for(size_t i = 0; i < query.size() - k + 1; ++i)
    {
        std::string kmer = query.substr(i, k);
        size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, pTargetBWT);
        if(count > 0)
        {
            position = i;
            std::cout << "SEED: " << kmer << " c: " << count << " p: " << position << "\n";
            return kmer;
        }
    }
    std::string empty;
    return empty;
}

// Search pTargetBWT for sequences that overlap the current set of hits
void LRCorrection::addOverlappingHits(const std::string& query,
                                      const BWT* pTargetBWT, 
                                      const SampledSuffixArray* pTargetSSA,
                                      const LRAlignment::LRParams& params,
                                      LRAlignment::LRHitVector& hits)
{
    (void)query;
    (void)pTargetBWT;
    (void)pTargetSSA;
    (void)params;

    // Map from a read ID to a hit
    typedef std::map<uint64_t, LRAlignment::LRHit> LRHitsMap;
    LRHitsMap hitsMap;
    std::queue<LRAlignment::LRHit*> hitsQueue;

    static const int minOverlap = 60;

    for(size_t i = 0; i < hits.size(); ++i)
    {
        LRHitsMap::iterator iter = hitsMap.insert(std::make_pair(hits[i].targetID, hits[i])).first;
        hitsQueue.push(&iter->second);
    }

    while(!hitsQueue.empty())
    {
        LRAlignment::LRHit* pCurrHit = hitsQueue.front();
        hitsQueue.pop();

        std::string w[2];
        w[0] = BWTAlgorithms::extractString(pTargetBWT, pCurrHit->targetID);
        w[1] = reverseComplement(w[0]);

        // Loop over the sequence and its reverse complement and find overlaps
        std::vector<uint64_t> newTargetIDs;

        for(size_t i = 0; i < 2; ++i)
        {
            std::string& x = w[i];
            int64_t l = x.size() - 1;
            BWTInterval interval;
            BWTAlgorithms::initInterval(interval, x[l], pTargetBWT);

            while(interval.isValid() && l >= 0)
            {
                int64_t o = x.size() - l;
                if(o >= minOverlap)
                {
                    BWTInterval prefix_interval = interval;
                    BWTAlgorithms::updateInterval(prefix_interval, '$', pTargetBWT);
                    
                    // Check if a valid overlap has been found. If so, add
                    // the reads to the map
                    if(prefix_interval.isValid())
                    {
                        for(int64_t k = prefix_interval.lower; k <= prefix_interval.upper; ++k)
                        {
                            size_t overlappedTargetID = pTargetSSA->lookupLexoRank(k);
                            newTargetIDs.push_back(overlappedTargetID);
                        }
                    }
                }

                // Perform a backwards extension
                BWTAlgorithms::updateInterval(interval, x[l], pTargetBWT);
                l -= 1;
            }
        }

        // Process the new targets found
        for(size_t i = 0; i < newTargetIDs.size(); ++i)
        {
            if(hitsMap.count(newTargetIDs[i]) == 0)
                std::cout << "Found new target: " << newTargetIDs[i] << "\n";
        }
    }
}
