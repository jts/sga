///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QCProcess - Process to perform quality checks 
// for a sequence work item
//
#ifndef QCPROCESS_H
#define QCPROCESS_H

#include "Util.h"
#include "BWT.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"

class QCResult
{
    public:
        QCResult() : qcPassed(false) {}

        bool qcPassed;
};

//
class QCProcess
{
    public:
        QCProcess(const BWT* pBWT, const BWT* pRBWT, int kmerLength, int kmerThreshold);
        ~QCProcess();
        QCResult process(const SequenceWorkItem& item);

    private:
        
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        const int m_kmerLength;
        const int m_kmerThreshold;
};

// Write the results from the overlap step to an ASQG file
class QCPostProcess
{
    public:
        QCPostProcess(std::ostream* pCorrectedWriter, std::ostream* pDiscardWriter);
        ~QCPostProcess();

        void process(const SequenceWorkItem& item, const QCResult& result);

    private:

        std::ostream* m_pCorrectedWriter;
        std::ostream* m_pDiscardWriter;

        size_t m_readsKept;
        size_t m_readsDiscarded;
};

#endif
