///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ErrorCorrectProcess - Wrapper to perform error correction
// for a sequence work item
//
#include "ErrorCorrectProcess.h"

//
//
//
ErrorCorrectProcess::ErrorCorrectProcess(const OverlapAlgorithm* pOverlapper, 
                                         int minOverlap) : m_pOverlapper(pOverlapper), 
                                                           m_minOverlap(minOverlap)
{

}

//
ErrorCorrectProcess::~ErrorCorrectProcess()
{

}

//
ErrorCorrectResult ErrorCorrectProcess::process(const SequenceWorkItem& workItem)
{
    OverlapResult overlap_result = m_pOverlapper->overlapRead(workItem.read, m_minOverlap, &m_blockList);
    m_blockList.clear();
    ErrorCorrectResult result;
    return result;
}

//
//
//
ErrorCorrectPostProcess::ErrorCorrectPostProcess(std::ostream* pWriter) : m_pWriter(pWriter)
{

}

//
void ErrorCorrectPostProcess::process(const SequenceWorkItem& /*item*/, const ErrorCorrectResult& /*result*/)
{
    //m_pOverlapper->writeResultASQG(*m_pASQGWriter, item.read, result);
}
