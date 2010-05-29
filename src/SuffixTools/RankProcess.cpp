///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RankProcess - Compute a vector of BWT ranks for
// SequenceWorkItems
//
#include "RankProcess.h"

//
//
//
RankProcess::RankProcess(const BWT* pBWT, bool doReverse) : m_pBWT(pBWT), m_doReverse(doReverse)
{

}

//
RankProcess::~RankProcess()
{

}

// Calculate the vector of ranks for the given sequence
RankVector RankProcess::process(const SequenceWorkItem& workItem)
{
    RankVector out;
    DNAString w = workItem.read.seq;
    if(m_doReverse)
        w.reverse();

    size_t l = w.length();
    int i = l - 1;

    // Compute the rank of the last character of seq. We consider $ to be implicitly
    // terminated by a $ character. The first rank calculated is for this and it is given
    // by the C(a) array in BWTInternal
    int64_t rank = m_pBWT->getPC('$'); // always zero
    out.push_back(rank);

    // Compute the starting rank for the last symbol of w
    char c = w.get(i);
    rank = m_pBWT->getPC(c);
    out.push_back(rank);
    --i;

    // Iteratively compute the remaining ranks
    while(i >= 0)
    {
        char c = w.get(i);
        rank = m_pBWT->getPC(c) + m_pBWT->getOcc(c, rank - 1);
        //std::cout << "c: " << c << " rank: " << rank << "\n";
        out.push_back(rank);
        --i;
    }
    return out;
}

//
//
//
RankPostProcess::RankPostProcess(GapArray* pGapArray) : m_pGapArray(pGapArray)
{

}

//
void RankPostProcess::process(const SequenceWorkItem& /*item*/, const RankVector& result)
{
    for(RankVector::const_iterator iter = result.begin(); iter != result.end(); ++iter)
        incrementGapArray(*iter, *m_pGapArray);
}
