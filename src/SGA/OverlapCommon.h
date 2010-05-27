//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapCommon - Common functions used in overlap methods
//
#ifndef OVERLAPCOMMON_H
#define OVERLAPCOMMON_H

#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"
#include "SGACommon.h"
#include "Timer.h"

namespace OverlapCommon
{

void parseHitsString(const std::string& hitString, 
                     const ReadTable* pFwdRT, 
                     const SuffixArray* pFwdSAI, const SuffixArray* pRevSAI, 
                     size_t& readIdx, OverlapVector& outVector, bool& isSubstring);
};

#endif
