//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRAlignment - Collection of algorithms for performing
// long read matches against an FM-index
//
#ifndef LRALIGNMENT_H
#define LRALIGNMENT_H

#include "BWT.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "HashMap.h"
#include "MultiAlignment.h"
#include "stdaln.h"
#include <stack>

namespace LRAlignment
{

// Enum of identifiers for cell-filtering heuristics
enum CutAlgorithm
{
    LRCA_DEFAULT,
    LRCA_Z_BEST,
    LRCA_Z_BEST_STRATA,
    LRCA_SCORE_FRAC,
    LRCA_NONE
};

// Parameters object
struct LRParams
{
    LRParams() { setDefaults(); }

    void setDefaults()
    {
        gap_open = 5;
        gap_ext = 2;
        mismatch = 3;
        match = 1;
        gap_open_extend = gap_open + gap_ext;
        threshold = 30;
        bandwidth = 50;
        zBest = 20;
        percentCutoff = 0.90f;
        cutTailAlgorithm = LRCA_Z_BEST_STRATA;
    }

    //
    void setPacBio()
    {
        setDefaults();

        gap_open = 2;
        gap_ext = 1;
        mismatch = 5;
        gap_open_extend = gap_open + gap_ext;
    }

    int gap_open;
    int gap_ext;
    int mismatch;
    int match;
    int gap_open_extend;
    int threshold;
    int bandwidth;
    
    // Cell filtering heuristics
    int zBest;
    double percentCutoff;
    CutAlgorithm cutTailAlgorithm;
};

// Structure holding an alignment between a query sequence and a BWT of a collection of sequences
struct LRHit
{
    LRHit() : interval(0,-1), flag(0), num_seeds(0), targetID(-1), t_start(-1), length(0), G(0), q_start(0), q_end(0) {}
    BWTInterval interval;
    std::string targetString;

    uint32_t flag;
    uint32_t num_seeds;

    uint64_t targetID; // ID of target sequence
    int t_start; // alignment start position on target
    int length; // length of the target alignment
    int G;
    int q_start; // alignment start on query
    int q_end; // alignment end (exclusive) on query

    static bool compareG(const LRHit& a, const LRHit& b) { return a.G > b.G; }
    static bool compareIDandG(const LRHit& a, const LRHit& b) 
    { 
        if(a.targetID == b.targetID)
            return a.G > b.G;
        else
            return a.targetID < b.targetID;
    }
};
typedef std::vector<LRHit> LRHitVector;

// Structure holding the score between a node of the prefix DAWG of the query
// and the prefix trie of the target
struct LRCell
{
    // Functions
    void initializeDefault();
    void clearChildren();
    bool hasUninitializedChild() const;

    // Data Members
    
    // interval on the target bwt
    BWTInterval interval;

    // scores
    int I;
    int D;
    int G;

    // Character code from parent
    uint8_t parent_cidx;

    // Query and target lengths
    int q_len;
    int t_len;

    // Index of the parent in cell array of the stack entry
    int parent_idx;

    // Index of this cell in the cell array
    int u_idx;

    // Indices of the children in the cell array
    int children_idx[4];

    // 
    std::string revTargetString;
};
typedef std::vector<LRCell> LRCellVector;

//
typedef HashMap<uint64_t, uint64_t> LRHash;

// A stack entry holds the an interval into the query BWT
// and an array of cells with the scores for that interval
struct LRStackEntry
{
    BWTInterval interval;
    LRCellVector cells;
};

//
typedef std::stack<LRStackEntry*> LRStack;
typedef std::vector<LRStackEntry*> LRPendingVector;

// Core alignment function - align the sequence query 
// against all sequences in pTargetBWT
void bwaswAlignment(const std::string& query, 
                    const BWT* pTargetBWT, 
                    const SampledSuffixArray* pTargetSSA,
                    const LRParams& params,
                    LRHitVector& outHits);

MultiAlignment convertHitsToMultiAlignment(const std::string& query, 
                                           const BWT* pTargetBWT, 
                                           const SampledSuffixArray* pTargetSSA,
                                           const LRParams& params,
                                           const LRHitVector& hits);


//
// Helper functions
//

// Initialize a hash table of BWT intervals representing
// the nodes in a DAWG
void initializeDAWGHash(BWT* pQueryBWT, LRHash& dawgHash);

// Merge the cells of the two stack entries
void mergeStackEntries(LRStackEntry* u, LRStackEntry* v);

// Update the given LRStack to contain the new StackEntry after
// performing any necessary merges with pending Stacks
int updateStack(LRStack* pStack, 
                LRStackEntry* u, 
                LRPendingVector* pPendingVector, 
                LRHash* pDawgHash, 
                const LRParams& params);

// Cull duplicated cells in the given stack entry
void removeDuplicateCells(LRStackEntry* u, LRHash& hash);

// Filter duplicated hits from hitsVector using their position on the query sequence
int resolveDuplicateHits(const BWT* pTargetBWT, 
                         const SampledSuffixArray* pTargetSSA, 
                         LRHitVector& hits, 
                         int IS);

// Filter hits out of the hits vector using their ID
// At most 1 hit per target sequence is kept
int resolveDuplicateHitsByID(const BWT* pTargetBWT, 
                             const SampledSuffixArray* pTargetSSA, 
                             LRHitVector& hits, 
                             int IS);

// add hits to the vector for cells that score above threshold
void saveHits(const SuffixArray* pQuerySA, 
              const SampledSuffixArray* pTargetSSA, 
              const BWT* pTargetBWT, 
              LRStackEntry* u, 
              int threshold, 
              LRHitVector& hits);

// save hits to cells that contain intervals that represent prefixes of reads
void saveTerminalHits(const SuffixArray* pQuerySA, 
                      const SampledSuffixArray* pTargetSSA, 
                      const BWT* pTargetBWT, 
                      LRStackEntry* u, 
                      int threshold, 
                      LRHitVector& hits);

// Using the cell points in c, calculate scores for cell c[0]
int fillCells(const LRParams& params, int match_score, LRCell* c[4]);

// Functions to heuristically remove low-scoring cells
void cutTail(LRStackEntry* u, const LRParams& params);
void cutTailByScorePercent(LRStackEntry* u, const LRParams& params);
void cutTailByZBest(LRStackEntry* u, const LRParams& params);
void cutTailByStratifiedZBest(LRStackEntry* u, const LRParams& params);

// Attempt to extend a hit coordinate to a full-length match
void extendHitFullLength(LRHit& hit, const std::string& query, 
                         const std::string& target, const LRParams& params, AlnParam* pStdAlnPar);

// Convert a dynamic programming path to a pair of padded strings representing the alignment
void path2padded(const std::string& s1, 
                 const std::string& s2, 
                 std::string& out1, 
                 std::string& out2, 
                 std::string& outm, 
                 path_t* path, 
                 int path_len);

// Convert a std::string into the stdAln required packed format.
// This function allocates memory which the caller must free
uint8_t* createStdAlnPacked(const std::string& s, size_t start = 0, size_t length = std::string::npos);

// Calculate the maximum target length for a query of length ql
size_t calculateMaxTargetLength(int ql, const LRParams& params);

};

#endif
