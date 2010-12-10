///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// FMMergeProcess - Merge unambiguously overlapping sequences
//
#include "FMMergeProcess.h"

//
FMMergeProcess::FMMergeProcess(const OverlapAlgorithm* pOverlapper, int minOverlap, const BitVector* pMarkedReads) : 
                                     m_pOverlapper(pOverlapper), 
                                     m_minOverlap(minOverlap), 
                                     m_pMarkedReads(pMarkedReads)
{

}

//
FMMergeProcess::~FMMergeProcess()
{


}

//
FMMergeResult FMMergeProcess::process(const SequenceWorkItem& item)
{
    (void)item;

    // Calculate the intervals in the forward FM-index for this read
    const BWT* pBWT = m_pOverlapper->getBWT();

    // Find the interval in the fm-index containing the read
    BWTInterval readInterval = BWTAlgorithms::findInterval(pBWT, item.read.seq.toString());

    // Update the interval by looking for the $ to map the interval into indices in the bit array
    BWTAlgorithms::updateInterval(readInterval, '$', pBWT);

    // Check if this read has been used yet
    bool used = false;
    for(int64_t i = readInterval.lower; i <= readInterval.upper; ++i)
    {
        used = m_pMarkedReads->test(i) || used;
    }

    assert(!used);

    FMMergeResult result;
    result.usedReads = readInterval;
    result.numReads = 0;
    return result;
}

//
FMMergePostProcess::FMMergePostProcess(std::ostream* pWriter, BitVector* pMarkedReads) : m_pWriter(pWriter), m_pMarkedReads(pMarkedReads)
{

}

//
void FMMergePostProcess::process(const SequenceWorkItem& item, const FMMergeResult& result)
{
    (void)item;
    (void)result;

    // Set a bit mask for the indicated values
    for(int64_t i = result.usedReads.lower; i <= result.usedReads.upper; ++i)
    {
        m_pMarkedReads->set(i, true);
    }

}
