///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StatsProcess - Compute statistics about the reads
//
#include "StatsProcess.h"
#include "BWTAlgorithms.h"

//
//
//
StatsProcess::StatsProcess(const BWT* pBWT, const BWT* pRBWT, int kmerLength) :
                            m_pBWT(pBWT),
                            m_pRBWT(pRBWT),
                            m_kmerLength(kmerLength)
{

}

//
StatsProcess::~StatsProcess()
{

}

//
StatsResult StatsProcess::process(const SequenceWorkItem& workItem)
{
    // Perform a kmer-based qc check on the read
    StatsResult result;

    std::string readSequence = workItem.read.seq.toString();
    int k = m_kmerLength;
    int n = readSequence.size();
    int nk = n - m_kmerLength + 1;

    for(int i = 0; i < nk; ++i)
    {
        std::string kmer = readSequence.substr(i, k);
        int count = BWTAlgorithms::countSequenceOccurrences(kmer, m_pBWT, m_pRBWT);
        result.kmerCoverage.push_back(count);
    }

    return result;
}

//
//
//
StatsPostProcess::StatsPostProcess()
{
}

//
StatsPostProcess::~StatsPostProcess()
{
    int max = 200;
    int maxCount = 0;
    printf("Kmer coverage histogram\n");
    printf("cov\tcount\n");
    for(std::map<int,int>::iterator iter = kmerCovHist.begin(); iter != kmerCovHist.end(); ++iter)
    {
        if(iter->first <= max)
            printf("%d\t%d\n", iter->first, iter->second);
        else
            maxCount += iter->second;
    }
    printf(">%d\t%d\n", max, maxCount);
}

//
void StatsPostProcess::process(const SequenceWorkItem& /*item*/, const StatsResult& result)
{
    for(size_t i = 0; i < result.kmerCoverage.size(); ++i)
    {
        kmerCovHist[result.kmerCoverage[i]]++;
    }
}
