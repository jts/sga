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
QCProcess::QCProcess(const BWT* pBWT, const BWT* pRBWT, BitVector* pSharedBV, bool checkDup, bool checkKmer, int kmerLength, int kmerThreshold) :
                     m_pBWT(pBWT),
                     m_pRBWT(pRBWT),
                     m_pSharedBV(pSharedBV),
                     m_checkDuplicate(checkDup),
                     m_checkKmer(checkKmer),
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
    QCResult result;

    if(m_checkDuplicate)
        result.dupPassed = performDuplicateCheck(workItem);
    else
        result.dupPassed = true;
    
    // Only perform the more expensive k-mer test if the dup check succeeded
    if(m_checkKmer && result.dupPassed)
        result.kmerPassed = performKmerCheck(workItem);
    else
        result.kmerPassed = true;


    return result;
}

bool QCProcess::performKmerCheck(const SequenceWorkItem& workItem)
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

    // Ensure the read is longer than the k-mer length
    if((int)w.size() < m_kmerLength)
        return false;

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

    return allSolid;
}

// Perform duplicate check
// Look up the interval of the read in the BWT. If the index of the read
bool QCProcess::performDuplicateCheck(const SequenceWorkItem& workItem)
{
    assert(m_pSharedBV != NULL);

    std::string w = workItem.read.seq.toString();
    std::string rc_w = reverseComplement(w);

    // Look up the interval of the sequence and its reverse complement
    BWTIntervalPair fwdIntervals = BWTAlgorithms::findIntervalPair(m_pBWT, m_pRBWT, w);
    BWTIntervalPair rcIntervals = BWTAlgorithms::findIntervalPair(m_pBWT, m_pRBWT, rc_w);

    // Check if this read is a substring of any other
    // This is indicated by the presence of a non-$ extension in the left or right direction
    AlphaCount64 fwdECL = BWTAlgorithms::getExtCount(fwdIntervals.interval[0], m_pBWT);
    AlphaCount64 fwdECR = BWTAlgorithms::getExtCount(fwdIntervals.interval[1], m_pRBWT);

    AlphaCount64 rcECL = BWTAlgorithms::getExtCount(rcIntervals.interval[0], m_pBWT);
    AlphaCount64 rcECR = BWTAlgorithms::getExtCount(rcIntervals.interval[1], m_pRBWT);

    if(fwdECL.hasDNAChar() || fwdECR.hasDNAChar() || rcECL.hasDNAChar() || rcECR.hasDNAChar())
    {
        // Substring reads are always removed so no need to update the bit vector
        return false;
    }

    // Calculate the lexicographic intervals for the fwd and reverse intervals
    BWTAlgorithms::updateBothL(fwdIntervals, '$', m_pBWT);
    BWTAlgorithms::updateBothL(rcIntervals, '$', m_pBWT);

    // Calculate the canonical index for this string - the lowest
    // value in the two lexicographic index
    int64_t fi = fwdIntervals.interval[0].isValid() ? fwdIntervals.interval[0].lower : std::numeric_limits<int64_t>::max();
    int64_t ri = rcIntervals.interval[0].isValid() ? rcIntervals.interval[0].lower : std::numeric_limits<int64_t>::max();
    int64_t canonicalIdx = std::min(fi, ri);

    // Check if the bit reprsenting the canonical index is set in the shared bit vector
    if(!m_pSharedBV->test(canonicalIdx))
    {
        // This read is not a duplicate
        // Attempt to atomically set the bit from false to true
        if(m_pSharedBV->updateCAS(canonicalIdx, false, true))
        {
            // Call succeed, return that this read is not a duplicate
            return true;
        }
        else
        {
            // Call failed, some other thread set the bit before
            // this thread. Return that the reead is a duplicate
            return false;
        }
    }
    else
    {
        // this read is duplicate
        return false;
    }
}

//
//
//
QCPostProcess::QCPostProcess(std::ostream* pCorrectedWriter,
                             std::ostream* pDiscardWriter) :
                                m_pCorrectedWriter(pCorrectedWriter),
                                m_pDiscardWriter(pDiscardWriter),
                                m_readsKept(0), m_readsDiscarded(0),
                                m_readsFailedKmer(0), m_readsFailedDup(0)
{

}

//
QCPostProcess::~QCPostProcess()
{
    std::cout << "Reads kept: " << m_readsKept << "\n";
    std::cout << "Reads discarded: " << m_readsDiscarded << "\n";
    std::cout << "Reads failed kmer check: " << m_readsFailedKmer << "\n";
    std::cout << "Reads failed duplicate check: " << m_readsFailedDup << "\n";
}

//
void QCPostProcess::process(const SequenceWorkItem& item, const QCResult& result)
{
    SeqRecord record = item.read;
    if(result.kmerPassed && result.dupPassed)
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

        if(!result.kmerPassed)
            m_readsFailedKmer += 1;
        else if(!result.dupPassed)
            m_readsFailedDup += 1;
    }
}
