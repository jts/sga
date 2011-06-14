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
    // Main correction function, dispatching the call to one 
    // of the actual correction algorithms
    std::string correct(const std::string& query,
                        const BWT* pTargetBWT, 
                        const BWT* pRevTargetBWT, 
                        const SampledSuffixArray* pTargetSSA,
                        const LRAlignment::LRParams& params);

    // Correct the query read by aligning it to all the reads in pTargetBWT
    // and calling a consensus sequence
    std::string correctAlignment(const std::string& query,
                                 const BWT* pTargetBWT, 
                                 const SampledSuffixArray* pTargetSSA,
                                 const LRAlignment::LRParams& params);
    
    // Correct the query read by breaking it into pieces, aligning those
    // pieces to pTargetBWT then constructing a multialignment from the sub-alignments.
    std::string correctAlignmentPartial(const std::string& query,
                                        const BWT* pTargetBWT, 
                                        const SampledSuffixArray* pTargetSSA,
                                        const LRAlignment::LRParams& params);

    // Correct a read by threading it through an implicit assembly graph
    // represented by pTargetBWT and pRevTargetBWT.
    std::string correctGraphThread(const std::string& query,
                                   const BWT* pTargetBWT, 
                                   const BWT* pRevTargetBWT, 
                                   const SampledSuffixArray* /*pTargetSSA*/,
                                   const LRAlignment::LRParams& /*params*/);


};

#endif
