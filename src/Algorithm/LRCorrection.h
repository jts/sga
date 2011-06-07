//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRCorrection - Collection of algorithms for 
// correcting long reads with a set of short reads
//
#ifndef LRCORRECTION_H
#define LRCORRECTION_H

#include "LRAlignment.h"

namespace LRCorrection
{
    // Main correction function
    std::string correctLongRead(const std::string& query,
                                const BWT* pTargetBWT, 
                                const BWT* pRevTargetBWT, 
                                const SampledSuffixArray* pTargetSSA,
                                const LRAlignment::LRParams& params);

    // 
    std::string correctAlignment(const std::string& query,
                                 const BWT* pTargetBWT, 
                                 const SampledSuffixArray* /*pTargetSSA*/,
                                 const LRAlignment::LRParams& /*params*/);


   //
   std::string correctGraphThread(const std::string& query,
                                  const BWT* pTargetBWT, 
                                  const BWT* pRevTargetBWT, 
                                  const SampledSuffixArray* /*pTargetSSA*/,
                                  const LRAlignment::LRParams& /*params*/);

    
    // Find the first portion of query that is found in the BWT
    std::string findSeedStringNaive(const std::string& query,
                                    const BWT* pTargetBWT, 
                                    int k,
                                    int& position);


    // Extend the hit vector by adding in reads that overlap the given set of reads
    void addOverlappingHits(const std::string& query,
                            const BWT* pTargetBWT, 
                            const SampledSuffixArray* pTargetSSA,
                            const LRAlignment::LRParams& params,
                            LRAlignment::LRHitVector& hits);

};

#endif
