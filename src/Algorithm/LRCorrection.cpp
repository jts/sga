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
#include "StringThreader.h"

std::string LRCorrection::correctLongRead(const std::string& query,
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
std::string LRCorrection::correctAlignment(const std::string& query,
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
        return query;
}


// Correct the read using an alignment-based approach
std::string LRCorrection::correctAlignmentPartial(const std::string& query,
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
        //addOverlappingHits(query, pTargetBWT, pTargetSSA, params, hits);        
    }

    /*
    MultiAlignment ma = LRAlignment::convertHitsToMultiAlignment(query, pTargetBWT, pTargetSSA, params, hits);
    std::string consensus = ma.generateConsensus();
    addOverlappingHits(query, pTargetBWT, pTargetSSA, params, hits);
    */
    return query;
}

// Attempt to correct the sequence via threading through a graph
std::string LRCorrection::correctGraphThread(const std::string& query,
                                             const BWT* pTargetBWT, 
                                             const BWT* pRevTargetBWT, 
                                             const SampledSuffixArray* /*pTargetSSA*/,
                                             const LRAlignment::LRParams& /*params*/)
{
    std::string seed_ecoli_test = "GATTTCCAGCGCGCCATCGCCACAGGCAATCAGCAGTGGCGCAACAGAAATCACGCTCCCCGGCTGTGCTTTGCTGGCATG";
    std::string seed_yeast_test = "TTTTATTTTTCCAAGAGTGCAATCAGCGGGTTTCCTCCTTATTTGCGTTTGGAGCAATCTTCCCTTTTTCGAACACAAAAT";
    std::string& seed = seed_yeast_test;
    StringThreader threader(seed, &query, seed.size(), 51, pTargetBWT, pRevTargetBWT);
    threader.run();
    threader.printAll();

    return query;
}
