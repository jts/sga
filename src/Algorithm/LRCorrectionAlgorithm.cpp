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
#include "HaplotypeBuilder.h"

namespace LRCorrectionAlgorithm
{

std::string correct(const std::string& query,
                    const BWT* pTargetBWT, 
                    const BWT* pRevTargetBWT, 
                    const SampledSuffixArray* pTargetSSA,
                    const LRAlignment::LRParams& params)
{
    (void)pTargetSSA;
    (void)params;
    //return correctAlignment(query, pTargetBWT, pTargetSSA, params);
    //return correctGraphThread(query, pTargetBWT, pRevTargetBWT, pTargetSSA, params);
    return correctDeBruijnPath(query, pTargetBWT, pRevTargetBWT);
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
                               const BWT* /*pRevTargetBWT*/, 
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
        StringThreader threader(seed, &querySeed, seedHit.q_end, extension_kmer, pTargetBWT);

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

// Correct the read using a de Bruijn graph threading approach
std::string correctDeBruijnPath(const std::string& query,
                                const BWT* pTargetBWT, 
                                const BWT* pRevTargetBWT,
                                bool verbose)
{
    size_t kmer = 41;
    size_t kmerThreshold = 5;

    if(verbose)
        std::cout << "Processing query: " << query << "\n";

    // Compute the set of anchor sequences. These are kmers with frequency
    // above the given threshold. For non-adjacent anchors, we construct haplotypes
    // from the graph and use these sequences to correct the read.
    AnchorVector anchors;

    // Create a table of kmer frequences over the query sequence
    for(size_t i = 0; i < query.size() - kmer + 1; ++i)
    {
        std::string ks = query.substr(i, kmer);
        size_t count = BWTAlgorithms::countSequenceOccurrences(ks, pTargetBWT);

        AnchorSequence anchor = { ks, static_cast<int>(i), static_cast<int>(count) };
        if(verbose)
            std::cout << anchor << "\n";
        if(count >= kmerThreshold)
            anchors.push_back(anchor);
    }

    // Build haplotypes between non-adjacent anchors

    // If there is zero or one anchor, do nothing
    if(anchors.size() <= 1)
        return "";

    AnchorSequence& previousAnchor = anchors.front();
    std::string correctedSequence = previousAnchor.sequence;

    for(size_t i = 1; i < anchors.size(); ++i)
    {
        AnchorSequence& currAnchor = anchors[i];
        
        // Check if the anchors are non-adjacent
        if(currAnchor.position != previousAnchor.position + 1)
        {
            // Build haplotypes between anchors
            if(verbose)
            {
                std::cout << "Anchor1: " << previousAnchor << "\n";
                std::cout << "Anchor2: " << currAnchor << "\n";
            }
            
            // Skip degenerate anchors, leave the read uncorrected
            if(previousAnchor.sequence == currAnchor.sequence)
                return "";

            HaplotypeBuilder builder;
            builder.setTerminals(previousAnchor, currAnchor);
            builder.setIndex(pTargetBWT, pRevTargetBWT);
            builder.setKmerParameters(kmer, kmerThreshold);

            HaplotypeBuilderResult result;
            builder.run();
            builder.parseWalks(result);

            if(result.haplotypes.empty())
                return ""; // no correction found
            
            if(verbose)
                std::cout << "Built " << result.haplotypes.size() << " haplotype\n";

            // Take a subsequence of the query that these haplotypes replace
            size_t queryStart = previousAnchor.position;
            size_t queryStop = currAnchor.position + kmer - 1;
            std::string querySub = query.substr(queryStart, queryStop - queryStart + 1);
            
            int64_t bestScore = 0;
            size_t bestIndex = 0;
            for(size_t j = 0; j < result.haplotypes.size(); ++j)
            {
                int score = StdAlnTools::globalAlignment(querySub, result.haplotypes[j], verbose);
                if(score > bestScore)
                {
                    bestScore = score;
                    bestIndex = j;
                }
            }

            // Append in the sequence of the best haplotype
            // We skip the first k bases as these are part
            // of the previous anchor, which has already been appended
            std::string toAppend = result.haplotypes[bestIndex].substr(kmer);
            correctedSequence.append(toAppend);
        }
        else
        {
            // Adjancent anchors, append in the last base of the anchor
            char lastBase = currAnchor.sequence[currAnchor.sequence.size() - 1];
            correctedSequence.append(1, lastBase);
        }
        previousAnchor = currAnchor;
    }

    if(verbose)
        StdAlnTools::globalAlignment(query, correctedSequence, true);
    return correctedSequence;
}

} // namespace
