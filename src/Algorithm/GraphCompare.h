///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphCompare - Compare two (abstractly represented)
// graphs against each other to find strings
// only present in one.
//
// The graphs are abstractly represented as
// an FM-index.
//
#ifndef GRAPH_COMPARE_H
#define GRAPH_COMPARE_H

#include <list>
#include <stack>
#include <queue>
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "BitVector.h"

// Typedefs
typedef std::vector<const BWT*> BWTVector;

struct GraphCompareStackNode
{
    // data

    static const size_t NUM_GRAPHS = 2;

    BWTIntervalPair intervalPairs[NUM_GRAPHS];
    AlphaCount64 lowerCounts[NUM_GRAPHS];
    AlphaCount64 upperCounts[NUM_GRAPHS];
    std::string str;
    uint8_t alphaIndex;
    int length;

    // functions

    // Initialize the intervals to be the range of all strings starting with b
    void initialize(char b, const BWTVector& bwts, const BWTVector& rbwts);

    // Update the intervals to extend to symbol b
    void update(char b, const BWTVector& bwts, const BWTVector& rbwts);
    
    //
    AlphaCount64 getAggregateExtCount() const;

    //
    void print() const;

};
typedef std::stack<GraphCompareStackNode> GraphCompareStack;

class GraphCompare
{
    public:

        //
        // Functions
        //
        GraphCompare(const BWT* pBaseBWT, const BWT* pBaseRBWT,
                     const BWT* pVariantBWT, const BWT* pVariantRBWT,
                     int kmer);

        ~GraphCompare();
        
        // Perform the comparison
        void run();

    private:
        
        bool processVariantKmer(const std::string& str, const BWTVector& bwts, const BWTVector& rbwts, int varIndex);

        // Mark all the kmers in seq (and their reverse complements as seen)
        void markVariantSequenceKmers(const std::string& str);

        // returns true if the given sequence is marked in the bitvector
        bool isKmerMarked(const std::string& str) const;

        //
        // Functions
        //

        //
        // Data
        //
        const BWT* m_pBaseBWT; 
        const BWT* m_pBaseRevBWT;
        const BWT* m_pVariantBWT; 
        const BWT* m_pVariantRevBWT;
        BitVector* m_pUsedVariantKmers;
        size_t m_kmer;
};

//
struct BubbleResult
{
    std::string targetString;
    std::string sourceString;
    bool success;
};

//
struct BubbleExtensionNode
{
    BubbleExtensionNode(Vertex* pX, EdgeDir d) : pVertex(pX), direction(d) {}

    Vertex* pVertex; // the vertex to extend
    EdgeDir direction; // the direction to extend to
};
typedef std::queue<BubbleExtensionNode> BubbleExtensionQueue;

//
// Class to build a variant bubble starting at a particular sequence
//
class BubbleBuilder
{
    public:
        BubbleBuilder();
        ~BubbleBuilder();

        // The source string is the string the bubble starts from
        void setSourceString(const std::string& str);

        // The source index is the index that the contains the source string
        void setSourceIndex(const BWT* pBWT, const BWT* pRBWT);

        // The target index is the index that we try to build the bubble onto
        void setTargetIndex(const BWT* pBWT, const BWT* pRBWT);

        // Run the bubble construction process
        // The found strings are placed in the StringVector
        // If this vector is empty, a bubble could not be found
        BubbleResult run();

    private:
        
        // Build the source portion of the graph
        bool buildSourceBubble();

        // Build the target portion of the graph
        bool buildTargetBubble();

        // After the bubble has been built into the graph, this function
        // finds and compares the two sequences
        BubbleResult parseBubble();

        // Returns true if the walk is the part of the target sequence
        bool classifyWalk(const SGWalk& walk) const;
    
        // Make the sequence of a new deBruijn vertex using the edge details
        std::string makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction);
        void addDeBruijnEdges(const Vertex* pX, const Vertex* pY, EdgeDir direction);

        const BWT* m_pSourceBWT;
        const BWT* m_pSourceRevBWT;

        const BWT* m_pTargetBWT;
        const BWT* m_pTargetRevBWT;

        StringGraph* m_pGraph;
        BubbleExtensionQueue m_queue;

        VertexPtrVec m_senseJoins;
        VertexPtrVec m_antisenseJoins;

        static const GraphColor SOURCE_COLOR = GC_BLUE;
        static const GraphColor TARGET_COLOR = GC_RED;
        static const GraphColor JOIN_COLOR = GC_BLACK;
};

#endif
