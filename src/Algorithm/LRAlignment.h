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
#include "stdaln.h"
#include <stack>

namespace LRAlignment
{

// 
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
    }

    int gap_open;
    int gap_ext;
    int mismatch;
    int match;
    int gap_open_extend;
    int threshold;
    int bandwidth;
    int zBest;
};

//
struct LRHit
{
    LRHit() : interval(0,-1), flag(0), num_seeds(0), targetID(-1), position(-1), length(0), G(0), G2(0), beg(0), end(0) {}
    BWTInterval interval;
    std::string targetString;

    uint32_t flag;
    uint32_t num_seeds;

    int targetID; // ID of target sequence
    int position; // alignment start position on target
    int length;
    int G;
    int G2;
    int beg; // alignment start on query
    int end; // alignment end (exclusive) on query

    static bool compareG(const LRHit& a, const LRHit& b) { return a.G > b.G; }
};
typedef std::vector<LRHit> LRHitVector;

//
struct LRCell
{
    // Functions
    void initializeDefault();
    void clearChildren();
    bool hasUninitializedChild() const;

    // Data Members
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
typedef HashMap<uint64_t, uint64_t> LRHash;

//
struct LRStackEntry
{
    BWTInterval interval;
    LRCellVector cells;
};

typedef std::stack<LRStackEntry*> LRStack;
typedef std::vector<LRStackEntry*> LRPendingVector;

//
void initializeDAWGHash(BWT* pQueryBWT, LRHash& dawgHash);

//
void bwaswAlignment(const std::string& query, 
                    const BWT* pTargetBWT, 
                    const SampledSuffixArray* pTargetSSA,
                    const LRParams& params);

//
void mergeStackEntries(LRStackEntry* u, LRStackEntry* v);

//
void removeDuplicateCells(LRStackEntry* u, LRHash& hash);

//
int resolveDuplicateHits(const BWT* pTargetBWT, 
                         const SampledSuffixArray* pTargetSSA, 
                         LRHitVector& hits, 
                         int IS);

//
void saveHits(const SuffixArray* pQuerySA, LRStackEntry* u, int threshold, LRHitVector& hits);

//
int fillCells(const LRParams& params, int match_score, LRCell* c[4]);

void cutTail(LRStackEntry* u, int T);

// Generate a cigar string for all hits in the vector
void generateCIGAR(const std::string& query, const LRParams& params, LRHitVector& hits);

void path2padded(const std::string& s1, const std::string& s2, std::string& out1, std::string& out2, std::string& outm, path_t* path, int path_len);

};

#endif
