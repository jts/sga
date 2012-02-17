///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadCoherentHaplotypeBuilder - Build read coherent haplotypes
// using walks through a de Bruijn graph.
//
#ifndef READ_COHERENT_HAPLOTYPE_BUILDER_H
#define READ_COHERENT_HAPLOTYPE_BUILDER_H
#include "BWT.h"
#include "SampledSuffixArray.h"
#include "BWTInterval.h"
#include "BWTIntervalCache.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "VariationBubbleBuilder.h"
#include "HaplotypeBuilder.h"
#include <queue>

// Build haplotypes starting from a given sequence.
class ReadCoherentHaplotypeBuilder
{
    public:

        ReadCoherentHaplotypeBuilder();
        ~ReadCoherentHaplotypeBuilder();

        void setInitialHaplotype(const std::string& str);
        void setIndex(const BWT* pBWT, const BWTIntervalCache* pCache, const SampledSuffixArray* pSSA);
        void setKmer(size_t k);

        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run();
        
    private:
        
        // Get all reads that share a kmer with a haplotype.
        void getReads();

        // Align the reads onto the haplotypes
        void alignReadsToHaplotypes();

        //
        // Data
        //
        std::list<std::string> m_haplotypes;
        std::vector<std::string> m_reads;

        const BWT* m_pBWT;
        const BWTIntervalCache* m_pIntervalCache;
        const SampledSuffixArray* m_pSSA;

        size_t m_kmer;
        size_t m_kmerThreshold;
        size_t m_maxEditDistance;
};

#endif
