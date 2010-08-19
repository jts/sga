///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ConnectProcess - Algorithm to resolve the sequence
// connecting two ends of a paired end read
//
#include "ConnectProcess.h"

//
//
//
ConnectProcess::ConnectProcess(const OverlapAlgorithm* pOverlapper, 
                               int minOverlap) : 
                                            m_pOverlapper(pOverlapper), 
                                            m_minOverlap(minOverlap)
{

}

//
ConnectProcess::~ConnectProcess()
{

}

//
ConnectResult ConnectProcess::process(const SequenceWorkItemPair& workItemPair)
{
    (void)workItemPair;
    std::cout << "Work pair: " << workItemPair.first.read.id << "," << workItemPair.second.read.id << "\n";
    ConnectResult result;

    return result;
}

//
//
//
ConnectPostProcess::ConnectPostProcess(std::ostream* pWriter) : 
                                                      m_pWriter(pWriter)
{

}

//
void ConnectPostProcess::process(const SequenceWorkItemPair& workItemPair, const ConnectResult& result)
{
    (void)workItemPair;
    (void)result;
}
