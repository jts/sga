///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HaplotypeBuilder - Construct candidate
// haplotypes from a pair of k-mer seeds.
//
#ifndef HAPLOTYPE_BUILDER_H
#define HAPLOTYPE_BUILDER_H
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "VariationBubbleBuilder.h"
#include <queue>

// The actual result structure
struct HaplotypeBuilderResult
{
    StringVector haplotypes;
};

//
// Class to build a variant bubble starting at a particular sequence
//
class HaplotypeBuilder
{
    public:

        HaplotypeBuilder();
        ~HaplotypeBuilder();

        void setTerminals(const std::string& leftTerminal, size_t leftCoverage, const std::string& rightTerminal, size_t rightCoverage);
        void setIndex(const BWT* pBWT, const BWT* pRBWT);

        // Set the threshold of kmer occurrences to use it as an edge
        void setKmerParameters(size_t k, size_t t);
    
        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        bool run();
        
        // Parse walks from the constructed graph
        void parseWalks(HaplotypeBuilderResult& results) const;

    private:
        
        // Add a vertex to the graph and record the sequence coverage value
        void addVertex(Vertex* pVertex, int coverage);

        // Make the sequence of a new deBruijn vertex using the edge details
        std::string makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction);
        void addDeBruijnEdges(const Vertex* pX, const Vertex* pY, EdgeDir direction);

        // Count the number of extensions of a de Bruijn node that are above
        // the required k-mer coverage
        size_t countValidExtensions(const AlphaCount64& ac) const;
        
        //
        // Data
        //
        const BWT* m_pBWT;
        const BWT* m_pRevBWT;

        StringGraph* m_pGraph;
        StrIntMap m_vertexCoverageMap;

        BubbleExtensionQueue m_queue;
        Vertex* m_pStartVertex;
        Vertex* m_pJoinVertex;
        
        //
        size_t m_kmerThreshold;
        size_t m_kmerSize;
};


#endif
