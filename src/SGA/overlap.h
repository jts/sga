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
#include "Match.h"
#include "BWTAlgorithms.h"


// functions

//
int overlapMain(int argc, char** argv);

// overlap computation
std::string computeHitsBWT();

// Overlap functions
size_t overlapReadExhaustive(size_t index, SeqItem& read, const BWT* pBWT, const BWT* pRBWT, HitVector* pHits, HitVector* pRevHits, OverlapBlockList* pOBOut);
size_t overlapReadIrreducible(SeqItem& read, const BWT* pBWT, const BWT* pRBWT, OverlapBlockList* pOBOut);

// Output processing
void outputHits(std::ofstream& handle, HitVector* pHits);
void parseHits(std::string hitsFile);
void writeOverlap(Overlap& overlap, std::ofstream& containHandle, std::ofstream& overlapHandle);

// utility functions
void flipCoords(const int len, int& s, int &e);
void swap(int& s, int& e);

// options
void parseOverlapOptions(int argc, char** argv);

#endif
