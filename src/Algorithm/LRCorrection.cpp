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

std::string LRCorrection::correctLongRead(const std::string& query,
                                          const BWT* pTargetBWT, 
                                          const SampledSuffixArray* pTargetSSA,
                                          const LRAlignment::LRParams& params)
{
    LRAlignment::LRHitVector hits;
    LRAlignment::bwaswAlignment(query, pTargetBWT, pTargetSSA, params, hits);
    MultiAlignment ma = LRAlignment::convertHitsToMultiAlignment(query, pTargetBWT, pTargetSSA, params, hits);
    std::string consensus = ma.generateConsensus();
    
    //ma.print(&consensus);

    return query;
}

