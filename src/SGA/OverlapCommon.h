//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapCommon - Common wrapper used for finding overlaps
// for a set of reads
//
#ifndef OVERLAPCOMMON_H
#define OVERLAPCOMMON_H

#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "LCPArray.h"
#include "SGACommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "AssembleExact.h"
#include "OverlapThread.h"
#include "ASQG.h"

namespace OverlapCommon
{

// overlap computation
size_t computeHitsSerial(const std::string& prefix, const std::string& readsFile, 
                         const OverlapAlgorithm* pOverlapper, OverlapMode mode,
                         int minOverlap, StringVector& filenameVec, std::ostream* pASQGWriter);

size_t computeHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                           const OverlapAlgorithm* pOverlapper, OverlapMode mode,
                           int minOverlap, StringVector& filenameVec, std::ostream* pASQGWriter);

};

#endif
