//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRCorrectionProcess - Wrapper class implementing
// the long-read correction process
//
#include "LRCorrectionProcess.h"
#include "LRCorrectionAlgorithm.h"

// Constructor
LRCorrectionProcess::LRCorrectionProcess(const LRCorrectionParameters params) : m_params(params)
{

}

// Destructor
LRCorrectionProcess::~LRCorrectionProcess()
{

}

// Perform the actual correction
LRCorrectionResult LRCorrectionProcess::process(const SequenceWorkItem& item)
{
    LRCorrectionResult result;
    if((int)item.read.seq.length() > m_params.minLength)
    {
        result.correctedSequence = LRCorrectionAlgorithm::correct(item.read.seq.toString(), 
                                                                  m_params.pBWT, 
                                                                  m_params.pRBWT,
                                                                  m_params.pSSA, 
                                                                  m_params.alignParams);
    }
    return result;
}

//
// Post processor
//
LRCorrectionPostProcess::LRCorrectionPostProcess(std::ostream* pCorrectedWriter) : m_pCorrectedWriter(pCorrectedWriter),
                                                                                   m_readsKept(0),
                                                                                   m_readsDiscarded(0)
{

}

//
LRCorrectionPostProcess::~LRCorrectionPostProcess()
{
    printf("Long read corrector kept %zu reads, discarded %zu\n", m_readsKept, m_readsDiscarded);
}

//
void LRCorrectionPostProcess::process(const SequenceWorkItem& item, const LRCorrectionResult& result)
{
    SeqRecord record = item.read;
    record.seq = result.correctedSequence;
    record.qual = "";

    if(!record.seq.empty())
    {
        record.write(*m_pCorrectedWriter);
        m_readsKept += 1;
    }
    else
    {
        m_readsDiscarded += 1;
    }
}
