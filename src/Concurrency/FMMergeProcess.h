///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// FMMergeProcess - Merge unambiguously overlapping sequences
//
#ifndef FMMERGEPROCESS_H
#define FMMERGEPROCESS_H

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"

struct FMMergeResult
{

    int numReads;
};

// Compute the overlap blocks for reads
class FMMergeProcess
{
    public:
        FMMergeProcess(const OverlapAlgorithm* pOverlapper, 
                       int minOverlap);

        ~FMMergeProcess();

        FMMergeResult process(const SequenceWorkItem& item);
    
    private:
        const OverlapAlgorithm* m_pOverlapper;
        const int m_minOverlap;
};

// Write the results from the overlap step to an ASQG file
class FMMergePostProcess
{
    public:
        FMMergePostProcess(std::ostream* pWriter);
        void process(const SequenceWorkItem& item, const FMMergeResult& result);

    private:
        std::ostream* m_pWriter;
};

#endif
