///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DeBruijnHaplotypeBuilder - Build haplotypes
// using a de Bruijn graph
//
#ifndef DEBRUIJN_HAPLOTYPE_BUILDER_H
#define DEBRUIJN_HAPLOTYPE_BUILDER_H

#include "BWT.h"
#include "SampledSuffixArray.h"
#include "BWTInterval.h"
#include "BWTIntervalCache.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "VariationBuilderCommon.h"
#include "HaplotypeBuilder.h"
#include "multiple_alignment.h"
#include "GraphCompare.h"
#include "ErrorCorrectProcess.h"
#include "SGWalk.h"
#include <queue>

// Build haplotypes starting from a given sequence.
class DeBruijnHaplotypeBuilder
{
    public:

        DeBruijnHaplotypeBuilder(const GraphCompareParameters& parameters);
        ~DeBruijnHaplotypeBuilder();
        
        // Set the string to start from
        void setInitialHaplotype(const std::string& str);

        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run(StringVector& out_haplotypes);
        
    private:

        //
        // Data
        //
        GraphCompareParameters m_parameters;
        std::string m_startingKmer;
};

#endif
