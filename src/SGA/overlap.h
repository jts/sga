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

// functions

//
int overlapMain(int argc, char** argv);

// overlap computation
void computeOverlapsLCP();
void computeOverlapsBWT();

//
OverlapVector processHits(size_t seqIdx, const HitVector& hitVec, const ReadTable* pFwdRT, const ReadTable* pRevRT);
OverlapVector alignRead(size_t seqIdx, const Sequence& seq, const BWT* pBWT, const ReadTable* pRT, bool isRevIdx, bool isReversed);
void processOverlaps(const OverlapVector& overlapVec, std::ofstream& containHandle, std::ofstream& overlapHandle);

// utility functions
void flipCoords(const int len, int& s, int &e);
void swap(int& s, int& e);
void writeContainment(std::ofstream& containHandle, const std::string& contained, const std::string& within);

// data structure creation
BWT* createBWT(SuffixArray* pSA, const ReadTable* pRT);
SuffixArray* loadSuffixArray(std::string filename);

// options
void parseOverlapOptions(int argc, char** argv);

#endif
