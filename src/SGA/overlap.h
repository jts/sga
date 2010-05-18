//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
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
void computeHitsBWT(OverlapMode mode, std::ostream* pASQGWriter, StringVector& filenameVec);

size_t computeHitsSerial(SeqReader& reader, OverlapMode mode, 
                         const OverlapAlgorithm* pOverlapper, 
                         StringVector& filenameVec,
                         std::ostream* pASQGWriter);

size_t computeHitsParallel(SeqReader& reader, OverlapMode mode, 
                           const OverlapAlgorithm* pOverlapper, 
                           StringVector& filenameVec,
                           std::ostream* pASQGWriter);

// Output processing
void convertHitsToASQG(const StringVector& hitsFilenames, std::ostream* pASQGWriter);

void convertHitsToOverlaps(const StringVector& hitsFilenames);
void writeOverlap(Overlap& overlap, std::ofstream& containHandle, std::ofstream& overlapHandle);

// Convert a line from a hits file into a vector of overlaps
OverlapVector hitStringToOverlaps(const std::string& hitString, 
                                  const ReadTable* pFwdRT, const ReadTable* pRevRT, 
                                  const SuffixArray* pFwdSAI, const SuffixArray* pRevSAI);

// options
void parseOverlapOptions(int argc, char** argv);

#endif
