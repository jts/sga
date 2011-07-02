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
    std::string reverseStr;
    uint8_t alphaIndex;
    size_t length;

    // functions

    // Initialize the intervals for the string str
    void initialize(const std::string& str, const BWTVector& bwts, const BWTVector& rbwts);

    // Update the intervals to extend to symbol b
    void update(char b, const BWTVector& bwts, const BWTVector& rbwts);
    
    // Return the count of each extension of str in the bwt collection
    AlphaCount64 getAggregateExtCount() const;

    //
    void print() const;

};
typedef std::stack<GraphCompareStackNode> GraphCompareStack;

// A range of kmers that each GraphCompare
// instance in multi-threaded mode will operate over
struct KmerRange
{
    KmerRange() : start(-1), end(-1) {}

    void initializeBatch(int batchIdx, int totalBatches);
    std::string decode(int64_t i) const;

    const static int keyLength = 4;
    int64_t start;
    int64_t end; // exclusive
};

// Parameters structure
class GraphCompareAggregateResults;
struct GraphCompareParameters
{
    const BWT* pBaseBWT; 
    const BWT* pBaseRevBWT;
    const BWT* pVariantBWT;
    const BWT* pVariantRevBWT;
    size_t kmer;
    BitVector* pBitVector;
    GraphCompareAggregateResults* pResults;

    // Parameters to control the range of k-mers to process
    KmerRange range;
};

//
// Statistics tracking object
//
struct GraphCompareStats
{
    // functions
    GraphCompareStats() { clear(); }
    void clear();
    void add(const GraphCompareStats& other);
    void print() const;

    // data
    int numBubbles;
    int numAttempted;
    int numTargetBranched;
    int numSourceBranched;
    int numTargetBroken;
    int numSourceBroken;
    int numWalkFailed;
    int numNoSolution;

    int numInsertions;
    int numDeletions;
    int numSubs;
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

        // 
        static void* runThreaded(void* obj);

    private:
        
        //
        // Functions
        //

        // Returns true if the kmer represented by the node is a variant
        bool isVariantKmer(GraphCompareStackNode* pNode) const;

        // Initialize the stack with a range of kmers
        void initializeStack(KmerRange range, GraphCompareStack& stack, const BWTVector& bwts, const BWTVector& rbwts);

        // Update the node stack by extending the given node
        bool updateNodeAndStack(GraphCompareStackNode* pNode, GraphCompareStack& stack, 
                                const BWTVector& bwts, const BWTVector& rbwts);

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
        // Data
        //
        GraphCompareParameters m_parameters;
        std::ostream* m_pWriter;

        // Results stats
        GraphCompareStats m_stats;
};

// Shared result object that the threaded
// GraphCompare instances write to. Protected
// by a mutex
class GraphCompareAggregateResults
{

    public:
        GraphCompareAggregateResults();
        ~GraphCompareAggregateResults();

        void updateShared(const GraphCompareStats stats);
        void printStats() const;

    private:
        pthread_mutex_t m_mutex;
        GraphCompareStats m_stats;
};

#endif
