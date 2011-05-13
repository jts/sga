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
#include <stack>

namespace LRAlignment
{


struct LRParams
{
    static const int gap_open = 5;
    static const int gap_ext = 3;
    static const int mismatch = 3;
    static const int match = 1;
    static const int gap_open_extend = gap_open + gap_ext;
};

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

//
struct LRStackEntry
{
    BWTInterval interval;
    LRCellVector cells;
};

typedef std::stack<LRStackEntry> LRStack;

void bwaswAlignment(const std::string& query, BWT* pTargetBWT);
int fill_cells(const LRParams& params, int match_score, LRCell* c[4]);

};

#endif
