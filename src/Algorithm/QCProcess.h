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
#include "BitVector.h"

// Parameters
struct QCParameters
{
    const BWT* pBWT;
    const BWT* pRevBWT;
    BitVector* pSharedBV;

    int kmerLength;
    int kmerThreshold;
    bool checkDuplicates;
    bool checkKmer;
};

// Results object
class QCResult
{
    public:
        QCResult() : kmerPassed(false), dupPassed(false) {}

        bool kmerPassed;
        bool dupPassed;
};

// Perform quality checks on the input stream of reads
class QCProcess
{
    public:
        QCProcess(QCParameters params); 
        ~QCProcess();
        QCResult process(const SequenceWorkItem& item);
        
        // Discard reads with low-frequency kmers
        bool performKmerCheck(const SequenceWorkItem& item);

        // Discard reads that are identical to, or a substring of, some other read
        bool performDuplicateCheck(const SequenceWorkItem& item);

    private:
        
        const QCParameters m_params;
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
        size_t m_readsFailedKmer;
        size_t m_readsFailedDup;
};

#endif
