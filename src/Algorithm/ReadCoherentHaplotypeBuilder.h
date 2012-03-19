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
#include "multiple_alignment.h"
#include <queue>

// Structure holding a read and its inferred position on the haplotype
struct HaplotypePositionedRead
{
    // functions
    friend bool operator<(const HaplotypePositionedRead& a, const HaplotypePositionedRead& b) { return a.position < b.position; }

    // data
    std::string sequence;
    int position; // relative to the initial kmer, which is position 0.
};

typedef std::vector<HaplotypePositionedRead> HaplotypeReadVector;

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
        HaplotypeBuilderReturnCode run(StringVector& out_haplotypes);
        
    private:
        
        // Get all reads that share a kmer with a haplotype.
        void getReadsWithKmer(const std::string& kmer, std::vector<std::string>* out_reads);
    
        // Build a multiple alignment of all the reads that are potentially part of this haplotype
        MultipleAlignment buildMultipleAlignment(HaplotypeReadVector& positioned_reads) const;

        // Compute a consensus sequence for the multiple alignment
        std::string getConsensus(MultipleAlignment* multiple_alignment, int max_differences) const;

        //
        // Data
        //
        const BWT* m_pBWT;
        const BWTIntervalCache* m_pIntervalCache;
        const SampledSuffixArray* m_pSSA;
        
        std::string m_initial_kmer_string;

        size_t m_kmer;
        size_t m_kmerThreshold;
        size_t m_maxEditDistance;
};

#endif
