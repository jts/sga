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
#include "SequenceProcessFramework.h"
#include "BWTIntervalCache.h"
#include "SampledSuffixArray.h"
#include "DindelRealignWindow.h"

// Structures and typedefs
typedef std::vector<const BWT*> BWTVector;

// Parameters structure
class GraphCompareAggregateResults;
struct GraphCompareParameters
{
    // Indices for the base reads
    const BWT* pBaseBWT; 
    const BWTIntervalCache* pBaseBWTCache;
    const SampledSuffixArray* pBaseSSA;
    
    // Indices for the variant reads
    const BWT* pVariantBWT;
    const BWTIntervalCache* pVariantBWTCache;
    const SampledSuffixArray* pVariantSSA;
    
    // Reference genome
    const BWT* pReferenceBWT;
    const BWT* pReferenceRevBWT;
    const SampledSuffixArray* pReferenceSSA;
    const ReadTable* pRefTable;

    // Bitvector to mark used kmers
    BitVector* pBitVector;
 
    // Parameters
    size_t kmer;
    size_t kmerThreshold;
    size_t maxBranches;
    bool bReferenceMode;
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
    int numHBFailed;
    int numNoSolution;

    int numInsertions;
    int numDeletions;
    int numSubs;
};

struct GraphCompareResult
{
    StringVector varStrings;
    StringVector baseStrings;
    DoubleVector varCoverages;
    DoubleVector baseCoverages;

    StringVector baseVCFStrings;
    StringVector variantVCFStrings;
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
        
        // Process a read and all its kmers
        GraphCompareResult process(const SequenceWorkItem& item);
        
        //
        void updateSharedStats(GraphCompareAggregateResults* pSharedStats);

    private:
        
        //
        // Functions
        //

        // When a kmer that is found in only one index, this function is called to attempt to build the full variation
        // string
        BubbleResult processVariantKmer(const std::string& str, int count, const BWTVector& bwts, int varIndex);
        BubbleResult processVariantKmerAggressive(const std::string& str, int count);
        
        // Mark all the kmers in str as being visited
        void markVariantSequenceKmers(const std::string& str);

        // Update statistics 
        void updateVariationCount(const BubbleResult& result);

        // Debug/testing functions
        bool buildVariantStringGreedy(const std::string& startingKmer, std::string& outString, size_t& flanking_k_length);
        bool buildVariantStringConservative(const std::string& startingKmer, std::string& outString, size_t& flanking_k_length);
        bool transformVariantString(const std::string& inStr, std::string& outStr);
        IntVector makeCountProfile(const std::string& str, size_t k, const BWT* pBWT, int max);


        //
        // Data
        //
        GraphCompareParameters m_parameters;

        // Results stats
        GraphCompareStats m_stats;

};

// Shared result object that the threaded
// GraphCompare instances write to. Protected
// by a mutex
class GraphCompareAggregateResults
{

    public:
        GraphCompareAggregateResults(const std::string& fileprefix);
        ~GraphCompareAggregateResults();

        void process(const SequenceWorkItem& item, const GraphCompareResult& result);

        void updateShared(const GraphCompareStats stats);
        void printStats() const;

    private:
        pthread_mutex_t m_mutex;
        GraphCompareStats m_stats;
        
        // Fasta output file
        std::ostream* m_pWriter;
        
        // VCF output files
        VCFFile m_baseVCFFile;
        VCFFile m_variantVCFFile;

        size_t m_numVariants;
};

#endif
