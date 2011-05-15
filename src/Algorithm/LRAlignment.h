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
#include "BWTAlgorithms.h"
#include "HashMap.h"
#include <stack>

namespace LRAlignment
{

// 
struct LRParams
{
    static const int gap_open = 5;
    static const int gap_ext = 2;
    static const int mismatch = 3;
    static const int match = 1;
    static const int gap_open_extend = gap_open + gap_ext;
    static const int threshold = 30;
};

//
struct LRHit
{
    LRHit() : interval(0,0), flag(0), num_seeds(0), length(0), G(0), G2(0), beg(0), end(0) {}
    BWTInterval interval;
    uint32_t flag;
    uint32_t num_seeds;
    int length;
    int G;
    int G2;
    int beg;
    int end;
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
void bwaswAlignment(const std::string& query, BWT* pTargetBWT);

//
void mergeStackEntries(LRStackEntry* u, LRStackEntry* v);

//
void removeDuplicateCells(LRStackEntry* u, LRHash& hash);

//
void saveHits(const SuffixArray* pQuerySA, LRStackEntry* u, int threshold, LRHitVector& hits);

//
int fillCells(const LRParams& params, int match_score, LRCell* c[4]);

};

#endif
