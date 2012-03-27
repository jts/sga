///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapHaplotypeBuilder - Build read coherent haplotypes
// using read overlaps.
//
#ifndef OVERLAP_HAPLOTYPE_BUILDER_H
#define OVERLAP_HAPLOTYPE_BUILDER_H
#include "BWT.h"
#include "SampledSuffixArray.h"
#include "BWTInterval.h"
#include "BWTIntervalCache.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "VariationBubbleBuilder.h"
#include "HaplotypeBuilder.h"
#include "multiple_alignment.h"
#include "GraphCompare.h"
#include <queue>

// Build haplotypes starting from a given sequence.
class OverlapHaplotypeBuilder
{
    public:

        OverlapHaplotypeBuilder();
        ~OverlapHaplotypeBuilder();

        void setInitialHaplotype(const std::string& str);
        void setParameters(const GraphCompareParameters& parameters); 

        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run(StringVector& out_haplotypes);
        
    private:
    
        // Get all the reads that contain at least one of the kmers in the vector
        void getReadsForKmers(const StringVector& kmers, StringVector* reads);

        // Coerece the set of reads into an ordered overlapping sequences
        // This function relies on the fact that all reads share the same kmer
        void orderReadsInitial(const std::string& initial_kmer, const StringVector& reads, StringVector* ordered_reads);

        // Order reads by inserting them into an already ordered list
        void orderReadsExtended(const StringVector& incoming_reads, StringVector* ordered_vector);
        
        // Remove duplicated reads from the ordered list
        void removeDuplicates(StringVector* ordered_vector);

        // Build a multiple alignment from an ordered set of reads
        MultipleAlignment buildMultipleAlignment(const StringVector& ordered_vector);

        // Extract new unique kmers from the reads in the set
        StringVector getExtensionKmers(const StringVector& reads);
        
        // Compute a consensus sequence for the multiple alignment
        std::string getConsensus(MultipleAlignment* multiple_alignment, int min_call_coverage, int max_differences) const;
    
        //
        // Data
        //
        GraphCompareParameters m_parameters;
        std::string m_initial_kmer_string;
        std::set<std::string> m_used_kmers;
};

#endif
