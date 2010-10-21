///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QCProcess - Process to perform quality checks 
// for a sequence work item
//
#include "QCProcess.h"
#include "BWTAlgorithms.h"

//
//
//
QCProcess::QCProcess(const BWT* pBWT, const BWT* pRBWT, int kmerLength, int kmerThreshold) :
                     m_pBWT(pBWT),
                     m_pRBWT(pRBWT),
                     m_kmerLength(kmerLength),
                     m_kmerThreshold(kmerThreshold)
{

}

//
QCProcess::~QCProcess()
{

}

//
QCResult QCProcess::process(const SequenceWorkItem& workItem)
{
    // Perform a kmer-based qc check on the read
    QCResult result;

    std::string readSequence = workItem.read.seq.toString();
    int k = m_kmerLength;
    int n = readSequence.size();
    int nk = n - m_kmerLength + 1;
    int threshold = m_kmerThreshold;

    // Are all kmers in the read well-represented?
    bool allSolid = true;

    for(int i = 0; i < nk; ++i)
    {
        std::string kmer = readSequence.substr(i, k);
        int count = BWTAlgorithms::countSequenceOccurrences(kmer, m_pBWT, m_pRBWT);
        if(count <= threshold)
        {
            allSolid = false;
            break;
        }
    }

    if(allSolid)
        result.qcPassed = true;
    else
        result.qcPassed = false;
    return result;
}

//
//
//
QCPostProcess::QCPostProcess(std::ostream* pCorrectedWriter,
                             std::ostream* pDiscardWriter) :
                                m_pCorrectedWriter(pCorrectedWriter),
                                m_pDiscardWriter(pDiscardWriter),
                                m_readsKept(0), m_readsDiscarded(0)
{

}

//
QCPostProcess::~QCPostProcess()
{
    std::cout << "Reads kept: " << m_readsKept << "\n";
    std::cout << "Reads discarded: " << m_readsDiscarded << "\n";
}

//
void QCPostProcess::process(const SequenceWorkItem& item, const QCResult& result)
{
    SeqRecord record = item.read;
    if(result.qcPassed)
    {
        record.write(*m_pCorrectedWriter);
        ++m_readsKept;
    }
    else
    {
        // To be able to rebuild the index after discarding the read, we need to write
        // the rank of the string (its position in the original read file into the read name)
        std::stringstream newID;
        newID << item.read.id << ",seqrank=" << item.idx;
        record.id = newID.str();

        record.write(*m_pDiscardWriter);
        ++m_readsDiscarded;
    }
}
