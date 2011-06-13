//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRCorrectionAlgorithm - Collection of algorithms for 
// correcting long reads with a set of short reads
//
#ifndef LRCORRECTION_H
#define LRCORRECTION_H

#include "LRAlignment.h"

namespace LRCorrectionAlgorithm
{
    // Main correction function
    std::string correct(const std::string& query,
                        const BWT* pTargetBWT, 
                        const BWT* pRevTargetBWT, 
                        const SampledSuffixArray* pTargetSSA,
                        const LRAlignment::LRParams& params);

    // Correct the query read by aligning it to all the reads in pTargetBWT
    std::string correctAlignment(const std::string& query,
                                 const BWT* pTargetBWT, 
                                 const SampledSuffixArray* pTargetSSA,
                                 const LRAlignment::LRParams& params);
    
    // Correct the query read by aligning it in overlapping pieces to pTargetBWT
    std::string correctAlignmentPartial(const std::string& query,
                                        const BWT* pTargetBWT, 
                                        const SampledSuffixArray* pTargetSSA,
                                        const LRAlignment::LRParams& params);


   // Experimental function for correcting a read by threading it through a graph
   std::string correctGraphThread(const std::string& query,
                                  const BWT* pTargetBWT, 
                                  const BWT* pRevTargetBWT, 
                                  const SampledSuffixArray* /*pTargetSSA*/,
                                  const LRAlignment::LRParams& /*params*/);


};

#endif
