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
#include "OverlapAlgorithm.h"

// functions

//
int overlapMain(int argc, char** argv);

// overlap computation
std::string computeHitsBWT();

size_t computeHitsSerial(SeqReader& reader, std::ofstream& writer, const OverlapAlgorithm* pOverlapper);
size_t computeHitsParallel(SeqReader& reader, const OverlapAlgorithm* pOverlapper);

//
void writeOverlapBlockList(std::ofstream& writer, size_t idx, const OverlapBlockList* pList);

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
