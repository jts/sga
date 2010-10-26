///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StatsProcess - Compute statistics about the reads
//
#ifndef STATSPROCESS_H
#define STATSPROCESS_H

#include "Util.h"
#include "BWT.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"

class StatsResult
{
    public:
        StatsResult() {}
        std::vector<int> kmerCoverage;
};

//
class StatsProcess
{
    public:
        StatsProcess(const BWT* pBWT, const BWT* pRBWT, int kmerLength);
        ~StatsProcess();
        StatsResult process(const SequenceWorkItem& item);

    private:
        
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        const int m_kmerLength;
};

// Write the results from the overlap step to an ASQG file
class StatsPostProcess
{
    public:
        StatsPostProcess();
        ~StatsPostProcess();

        void process(const SequenceWorkItem& item, const StatsResult& result);

    private:

        std::map<int, int> kmerCovHist;
};

#endif
