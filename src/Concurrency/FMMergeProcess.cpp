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
FMMergeProcess::FMMergeProcess(const OverlapAlgorithm* pOverlapper, int minOverlap) : m_pOverlapper(pOverlapper), m_minOverlap(minOverlap)
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

    FMMergeResult result;
    result.numReads = 0;
    return result;
}

//
FMMergePostProcess::FMMergePostProcess(std::ostream* pWriter) : m_pWriter(pWriter)
{

}

//
void FMMergePostProcess::process(const SequenceWorkItem& item, const FMMergeResult& result)
{
    (void)item;
    (void)result;
}
