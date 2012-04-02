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
#include "ErrorCorrectProcess.h"
#include "SGWalk.h"
#include <queue>

// Build haplotypes starting from a given sequence.
class OverlapHaplotypeBuilder
{
    public:

        OverlapHaplotypeBuilder(const GraphCompareParameters& parameters);
        ~OverlapHaplotypeBuilder();

        void setInitialHaplotype(const std::string& str);

        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run(StringVector& out_haplotypes);
        
    private:
    
        // Check for walks through the graph that cover all seed vertices
        // If there is a valid covering assembly of the reads, it will be
        // written to the vector
        void checkWalks(StringVector* out_strings);

        // Get all the reads that contain at least one of the kmers in the vector
        void getReadsForKmers(const StringVector& kmers, size_t k, StringVector* reads);

        // Correct sequencing errors in the set of reads
        void correctReads(StringVector* reads);

        // Initialize the string graph using the corrected reads
        void buildInitialGraph(const StringVector& reads);

        // Extend the graph by finding overlaps for the tip vertices
        void extendGraph();

        // Insert a new vertex into the graph with the specified sequence
        void insertVertexIntoGraph(const std::string& prefix, const std::string& sequence);

        // Find corrected reads that share a perfect overlap to the input sequence
        StringVector getCorrectedOverlaps(const std::string& sequence);

        // Find tip vertices in the graph
        VertexPtrVec findTips() const;

        // Returns true if the sequence represents a junction in the variation graph
        bool isJoinSequence(const std::string& sequence);

        // Returns the number of vertices in the vector that are contained in any walk
        size_t countCoveredVertices(const VertexPtrVec& vertices, const SGWalkVector& walks);

        // Coerece the set of reads into an ordered overlapping sequences
        // This function relies on the fact that all reads share the same kmer
        void orderReadsInitial(const std::string& initial_kmer, const StringVector& reads, StringVector* ordered_reads);

        // Order reads by inserting them into an already ordered list
        void orderReadsExtended(const StringVector& incoming_reads, StringVector* ordered_vector);
        
        // Remove duplicated reads from the ordered list
        void removeDuplicates(StringVector* ordered_vector);

        // Build a multiple alignment from an ordered set of reads
        MultipleAlignment buildMultipleAlignment(const StringVector& ordered_vector);

        //
        // Data
        //
        GraphCompareParameters m_parameters;
        ErrorCorrectProcess* m_corrector;
        StringGraph* m_graph;
        size_t m_numReads;
        static const int MIN_OVERLAP = 61;

        std::string m_initial_kmer_string;
        std::set<std::string> m_used_reads;

        std::map<std::string, std::string> m_correction_cache;
};

#endif
