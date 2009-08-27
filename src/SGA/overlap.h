//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// overlap - Overlap reads using a bwt
//
#ifndef OVERLAP_H
#define OVERLAP_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"

// typedefs
typedef std::vector<Overlap> OverlapVector;

// functions

//
int overlapMain(int argc, char** argv);

// overlap computation
void computeOverlaps();
OverlapVector processHits(size_t seqIdx, const HitVector& hitVec, const ReadTable* pFwdRT, const ReadTable* pRevRT);
OverlapVector alignRead(size_t seqIdx, const Sequence& seq, const BWT* pBWT, const ReadTable* pRT, bool isRevIdx, bool isReversed);
void flipCoords(const int len, int& s, int &e);
void swap(int& s, int& e);

// data structure creation
BWT* createBWT(SuffixArray* pSA, const ReadTable* pRT);
SuffixArray* loadSuffixArray(std::string filename);

// options
void parseOverlapOptions(int argc, char** argv);

#endif
