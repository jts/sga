///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VariationBubbleBuilder - Construct a variation
// bubble from an initial seed k-mer which only
// appears in one out of a pair of abstract
// deBruijn graphs.
//
#ifndef VARIATION_BUBBLE_BUILDER_H
#define VARIATION_BUBBLE_BUILDER_H
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include <queue>

// Result of the bubble construction
// Used to track statistics
enum BubbleResultCode
{
    BRC_UNKNOWN,
    BRC_OK,
    BRC_SOURCE_BROKEN,
    BRC_SOURCE_BRANCH,
    BRC_TARGET_BROKEN,
    BRC_TARGET_BRANCH,
    BRC_WALK_FAILED,
    BRC_NO_SOLUTION
};

// The actual result structure
struct BubbleResult
{
    std::string targetString;
    std::string sourceString;
    double targetCoverage;
    double sourceCoverage;
    BubbleResultCode returnCode;
};

// A directed node in the de Bruijn graph that
// has not had it's neighbors visited
struct BuilderExtensionNode
{
    BuilderExtensionNode(Vertex* pX, EdgeDir d) : pVertex(pX), direction(d) {}

    Vertex* pVertex; // the vertex to extend
    EdgeDir direction; // the direction to extend to
};
typedef std::queue<BuilderExtensionNode> BuilderExtensionQueue;
typedef std::map<std::string, int> StrIntMap;

//
// Class to build a variant bubble starting at a particular sequence
//
class VariationBubbleBuilder
{
    public:
        VariationBubbleBuilder();
        ~VariationBubbleBuilder();

        // The source string is the string the bubble starts from
        void setSourceString(const std::string& str, int coverage);

        // The source index is the index that the contains the source string
        void setSourceIndex(const BWT* pBWT, const BWT* pRBWT);

        // The target index is the index that we try to build the bubble onto
        void setTargetIndex(const BWT* pBWT, const BWT* pRBWT);

        // Set the threshold of kmer occurrences to use it as an edge
        void setKmerThreshold(size_t t);
    
        // Set the maximum number of allowed branches when searching
        // for the completion of the target half of the bubble
        void setAllowedBranches(size_t b);

        // Run the bubble construction process
        // The found strings are placed in the StringVector
        // If this vector is empty, a bubble could not be found
        BubbleResult run();
        
        // Get all the kmers on the source branch
        StringVector getSourceKmers() const;

    private:
        
        // Build the source portion of the graph
        BubbleResultCode buildSourceBubble();

        // Build the target portion of the graph
        BubbleResultCode buildTargetBubble();

        // After the bubble has been built into the graph, this function
        // finds and compares the two sequences
        void parseBubble(BubbleResult& result);

        // Returns true if the walk is the part of the target sequence
        // Also calculates the kmer coverage of the walk
        bool classifyWalk(const SGWalk& walk, int& outCoverage) const;
    
        // Add a vertex to the graph and record the sequence coverage value
        void addVertex(Vertex* pVertex, int coverage);

        // Make the sequence of a new deBruijn vertex using the edge details
        void addDeBruijnEdges(const Vertex* pX, const Vertex* pY, EdgeDir direction);

        //
        // Data
        //
        const BWT* m_pSourceBWT;
        const BWT* m_pSourceRevBWT;

        const BWT* m_pTargetBWT;
        const BWT* m_pTargetRevBWT;

        StringGraph* m_pGraph;
        StrIntMap m_vertexCoverageMap;

        BuilderExtensionQueue m_queue;

        VertexPtrVec m_senseJoins;
        VertexPtrVec m_antisenseJoins;
        
        //
        size_t m_kmerThreshold;
        size_t m_allowedTargetBranches;

        static const GraphColor SOURCE_COLOR = GC_BLUE;
        static const GraphColor TARGET_COLOR = GC_RED;
        static const GraphColor JOIN_COLOR = GC_BLACK;
};


#endif
