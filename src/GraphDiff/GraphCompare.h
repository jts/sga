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
#include "BWTIndexSet.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "BitVector.h"
#include "VariationBuilderCommon.h"
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
    // Base read index
    BWTIndexSet baseIndex;
    
    // Indices for the variant reads
    BWTIndexSet variantIndex;
    
    // Reference genome
    BWTIndexSet referenceIndex;
    const ReadTable* pRefTable;

    // Bitvector to mark used kmers
    BitVector* pBitVector;
 
    // Parameters
    size_t kmer;
    size_t kmerThreshold;
    size_t maxKmerThreshold; // skip kmers seen this many times or more
    size_t maxBranches;
    int maxSingletons; // the maximum number of times we can use k-mers with count = 1 in the extension
    size_t minKmerThreshold;
    int minOverlap;
    bool bReferenceMode;
    int verbose;

    DindelRealignParameters dindelRealignParameters;
};

struct GraphBuildResult
{
    StringVector variant_haplotypes;
    StringVector base_haplotypes;
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
        void debug(const std::string& debugFilename);
        void testKmersFromFile(const std::string& kmerFilename);
        void testKmer(const std::string& kmer);

        //
        void updateSharedStats(GraphCompareAggregateResults* pSharedStats);

    private:
        
        //
        // Functions
        //

        // Attempt to assemble a variant kmer into haplotypes
        GraphBuildResult processVariantKmer(const std::string& str, int count);
        
        // Generate a bitmask of kmers that have already been used 
        std::vector<bool> generateKmerMask(const std::string& str) const;

        // Mark all the kmers in str as being visited
        void markVariantSequenceKmers(const std::string& str);
        
        // Calculate the largest k such that every k-mer in the sequence is present at least min_depth times in the BWT
        size_t calculateMaxCoveringK(const std::string& sequence, int min_depth, const BWTIndexSet& indices);

        // Calculate the number of high coverage branches off a haplotype path through the de Bruijn graph
        size_t calculateHaplotypeBranches(const std::string& sequence, size_t k, size_t min_branch_depth, const BWTIndexSet& indices);

        // Update statistics 
        void updateVariationCount(const BubbleResult& result);

        // Debug/testing functions
        bool buildVariantStringGraph(const std::string& startingKmer, StringVector& haplotypes);

        bool transformVariantString(const std::string& inStr, std::string& outStr);
        IntVector makeCountProfile(const std::string& str, size_t k, const BWT* pBWT, int max);
        void showMappingLocations(const std::string& str);

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
