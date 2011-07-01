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
#include "VariationBubbleBuilder.h"

// Structures and typedefs
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
    size_t length;

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

// Parameters structure
struct GraphCompareParameters
{
    const BWT* pBaseBWT; 
    const BWT* pBaseRevBWT;
    const BWT* pVariantBWT;
    const BWT* pVariantRevBWT;
    size_t kmer;
};

//
//
//
class GraphCompare
{
    public:

        //
        // Functions
        //
        GraphCompare(const GraphCompareParameters& params);

        ~GraphCompare();
        
        // Perform the comparison
        void run();

    private:
        
        // Returns true if the kmer represented by the node is a variant
        bool isVariantKmer(GraphCompareStackNode* pNode) const;

        // Update the node stack by extending the given node
        bool updateNodeAndStack(GraphCompareStackNode* pNode, GraphCompareStack& stack, const BWTVector& bwts, const BWTVector& rbwts);

        // When a kmer that is found in only one index, this function is called to attempt to build the full variation
        // string
        bool processVariantKmer(const std::string& str, const BWTVector& bwts, const BWTVector& rbwts, int varIndex);

        // Mark all the kmers in seq (and their reverse complements as seen)
        void markVariantSequenceKmers(const std::string& str);

        // returns true if the given sequence is marked in the bitvector
        bool isKmerMarked(const std::string& str) const;

        // Update statistics 
        void updateVariationCount(const BubbleResult& result);
        
        //
        // Functions
        //

        //
        // Data
        //
        GraphCompareParameters m_parameters;
        BitVector* m_pUsedVariantKmers;
        std::ostream* m_pWriter;

        // Results stats
        int m_numBubbles;
        int m_numAttempted;
        int m_numTargetBranched;
        int m_numSourceBranched;
        int m_numTargetBroken;
        int m_numSourceBroken;
        int m_numWalkFailed;
        int m_numNoSolution;

        int m_numInsertions;
        int m_numDeletions;
        int m_numSubs;
};

#endif
