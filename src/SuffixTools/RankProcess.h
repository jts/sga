///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RankProcess - Compute a vector of BWT ranks for
// SequenceWorkItems
//
#ifndef RANKPROCESS_H
#define RANKPROCESS_H

#include "Util.h"
#include "BWT.h"
#include "SequenceWorkItem.h"
#include "GapArray.h"

typedef std::vector<int64_t> RankVector;

// Compute the overlap blocks for reads
class RankProcess
{
    public:
        RankProcess(const BWT* pBWT, bool doReverse);
        ~RankProcess();

        RankVector process(const SequenceWorkItem& item);
    
    private:
        const BWT* m_pBWT;
        bool m_doReverse;
};

// Update the gap array with 
class RankPostProcess
{
    public:
        RankPostProcess(GapArray* pGapArray);
        void process(const SequenceWorkItem& item, const RankVector& rankVec);

    private:
        GapArray* m_pGapArray;
};

#endif
