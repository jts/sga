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
struct KmerWindow
{
    int64_t getCount() const
    {
        int64_t count = 0;
        if(fwdIntervals.interval[0].isValid())
        {
            count += fwdIntervals.interval[0].size();
        }
        
        if(rcIntervals.interval[0].isValid())
        {
            count += rcIntervals.interval[0].size();
        }
        return count;
    }

    int start;
    int end;
    BWTIntervalPair fwdIntervals;
    BWTIntervalPair rcIntervals;

    bool isInitialized;
};

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
    // Perform a k-mer filter on the read
    // Each k-mer must be seen at least m times for the read to be kept
    // Naively, this would required k*(l - k + 1) calls to the BWT
    // if each k-mer was computed independently. We improve this by performing
    // a two-stage test. If the kmer k_i is seen at least m times, we extend
    // it by one base. If this k+1-mer is also seen at least m times, it implies
    // that the kmer k_(i+1) is also seen at least m times so it does not need
    // to be recomputed. If the (k+1)-mer is seen less than m times, we recompute
    // kmer k_(i+1) as a final check. Using this optimistic algorithm (checking
    // longer kmers than is required) is significantly faster.
    QCResult result;

    std::string w = workItem.read.seq.toString();
    int k = m_kmerLength;
    int n = w.size();
    int nk = n - m_kmerLength + 1;
    int threshold = m_kmerThreshold;

    // Are all kmers in the read well-represented?
    bool allSolid = true;

    int i = 0;

    KmerWindow window;
    window.start = 0;
    window.end = 0;
    window.isInitialized = false;
    
    while(i < nk)
    {
        // Determine whether the kmer at position i
        // has been seen more than m times.
        if(window.isInitialized)
        {
            // Extend the current window by 1 base
            // to see if the spanning larger
            // kmer is seen more than m times.
            // If so, we don't need to check this kmer
            // individually
            int next = window.end + k;
            assert(next < (int)w.size());
            char b = w[next];
            char cb = complement(b);

            if(window.fwdIntervals.interval[0].isValid())
                BWTAlgorithms::updateBothR(window.fwdIntervals, b, m_pRBWT);
            if(window.rcIntervals.interval[1].isValid())
                BWTAlgorithms::updateBothR(window.rcIntervals, cb, m_pBWT);

            int64_t count = window.getCount();
            if(count <= threshold)
            {
                // The extended kmer didn't meet the threshold
                // recompute the interval for this kmer
                window.isInitialized = false;
            }
            else
            {
                window.end += 1;
            }
        }
        
        if(!window.isInitialized)
        {
            // initialize the window by computing the
            // BWTIntervals for the kmer starting at
            // i and its reverse complement
            char b = w[i];
            char cb = complement(b);
            BWTAlgorithms::initIntervalPair(window.fwdIntervals, b, m_pBWT, m_pRBWT);
            BWTAlgorithms::initIntervalPair(window.rcIntervals, cb, m_pRBWT, m_pBWT);

            for(int j = i + 1; j < i + k; ++j)
            {
                // Update intervals rightwards 
                b = w[j];
                cb = complement(b);

                if(window.fwdIntervals.interval[0].isValid())
                    BWTAlgorithms::updateBothR(window.fwdIntervals, b, m_pRBWT);
                if(window.rcIntervals.interval[1].isValid())
                    BWTAlgorithms::updateBothR(window.rcIntervals, cb, m_pBWT);
            }

            // record the start/end indices of the kmers spanned by this position
            window.start = i;
            window.end = i;
            window.isInitialized = true;
        }

        // Check the count of this interval
        int64_t count = window.getCount();
            
        if(count <= threshold)
        {
            allSolid = false;
            break;
        }
        else
        {
            i += 1;
        }
    }

    /*
    for(int i = 0; i < nk; ++i)
    {
        std::string kmer = w.substr(i, k);
        int count = BWTAlgorithms::countSequenceOccurrences(kmer, m_pBWT, m_pRBWT);
        if(count <= threshold)
        {
            allSolid = false;
            break;
        }
    }
    */
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
