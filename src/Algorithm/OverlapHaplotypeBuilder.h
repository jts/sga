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
        void addInitialReadsToGraph(const StringVector& reads);

        // Add a single vertex to the graph
        Vertex* addVertexToGraph(const std::string& id, const std::string& sequence);
        
        // Insert a new vertex into the graph, which is an extension of source
        Vertex* addExtensionVertexToGraph(Vertex* source, const std::string& sequence);

        // Returns true if this vertex is a node that contains kmer(s) also present
        // in the base reads/reference. If so, we can merge the sequence back into the
        // graph at this point
        bool isVertexJoinNode(const Vertex* vertex) const;

        // Recruit new edges into the graph
        VertexPtrVec extendGraph();

        // Check if there are paths that form candidate haplotypes that cover the base vertices
        std::string findHaplotypes();

        // Remove duplicate reads and containments from the graph
        void cleanDuplicates();

        // Remove sequences from the vector that are already
        // present in the graph
        void removeUsedSequences(StringVector* sequences);

        // Build a multiple alignment from a walk through a graph
        MultipleAlignment buildMultipleAlignmentFromWalk(const SGWalk& walk);
        
        // Compute a consensus sequence for the multiple alignment
        std::string getConsensus(MultipleAlignment* multiple_alignment, int min_call_coverage, int max_differences) const;
    
        // Get overlapping reads using the FM-index
        StringVector getOverlappingReads(const std::string& str) const;

        // Write the graph to disk as a dot file
        void writeGraph(const std::string& filename) const;

        //
        // Data
        //
        StringGraph* m_graph;
        VertexPtrVec m_seedVertices;
        VertexPtrVec m_joinVertices[2];
        
        std::set<std::string> m_used_sequences;

        double m_minIdentity;
        int m_minOverlap;

        GraphCompareParameters m_parameters;
        std::string m_initial_kmer_string;
        size_t m_extend_vertices;
};

#endif
