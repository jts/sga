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
void computeHitsBWT(StringVector& filenames);
size_t computeHitsSerial(SeqReader& reader, const OverlapAlgorithm* pOverlapper, StringVector& filenames);
size_t computeHitsParallel(SeqReader& reader, const OverlapAlgorithm* pOverlapper, StringVector& filenames);
size_t computeHitsParallelBatch(SeqReader& reader, const OverlapAlgorithm* pOverlapper, StringVector& filenames);


// Output processing
void convertHitsToASQG(const StringVector& hitsFilenames);
void convertHitsToOverlaps(const StringVector& hitsFilenames);

void writeOverlap(Overlap& overlap, std::ofstream& containHandle, std::ofstream& overlapHandle);

// Convert a line from a hits file into a vector of overlaps
OverlapVector hitStringToOverlaps(const std::string& hitString, 
                                  const ReadTable* pFwdRT, const ReadTable* pRevRT, 
								  const SuffixArray* pFwdSAI, const SuffixArray* pRevSAI);

// utility functions
void flipCoords(const int len, int& s, int &e);
void swap(int& s, int& e);

// options
void parseOverlapOptions(int argc, char** argv);

#endif
