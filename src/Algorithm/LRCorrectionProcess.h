//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRCorrectionProcess - Wrapper class implementing
// the long-read correction process
//
#ifndef LRCORRECTION_PROCESS_H
#define LRCORRECTION_PROCESS_H

#include <string>
#include "SequenceProcessFramework.h"
#include "BWT.h"
#include "SampledSuffixArray.h"
#include "LRAlignment.h"

// Parameters to the corrector
class LRCorrectionParameters
{
    public:
        const BWT* pBWT;
        const BWT* pRBWT;
        const SampledSuffixArray* pSSA;
        LRAlignment::LRParams alignParams;
};

// Object holding the result of correcting one read
class LRCorrectionResult
{
    public:
        LRCorrectionResult() {}

        // data
        std::string correctedSequence;

};

//
class LRCorrectionProcess
{
    public:
        LRCorrectionProcess(const LRCorrectionParameters params); 
        ~LRCorrectionProcess();

        LRCorrectionResult process(const SequenceWorkItem& item);

    private:
        
        LRCorrectionParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class LRCorrectionPostProcess
{
    public:
        LRCorrectionPostProcess(std::ostream* pCorrectedWriter);
        ~LRCorrectionPostProcess();

        void process(const SequenceWorkItem& item, const LRCorrectionResult& result);

    private:

        std::ostream* m_pCorrectedWriter;

        size_t m_readsKept;
        size_t m_readsDiscarded;
};



#endif
