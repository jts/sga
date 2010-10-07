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
#include "ErrorCorrect.h"

//#define KMER_TESTING 1

//
//
//
ErrorCorrectProcess::ErrorCorrectProcess(const OverlapAlgorithm* pOverlapper, 
                                         int minOverlap, int numRounds, 
                                         int conflictCutoff, ErrorCorrectAlgorithm algo,
                                         bool printMO) : 
                                            m_pOverlapper(pOverlapper), 
                                            m_minOverlap(minOverlap),
                                            m_numRounds(numRounds),
                                            m_conflictCutoff(conflictCutoff),
                                            m_algorithm(algo),
                                            m_printOverlaps(printMO),
                                            m_depthFilter(10000)
{
}

//
ErrorCorrectProcess::~ErrorCorrectProcess()
{

}

//
ErrorCorrectResult ErrorCorrectProcess::process(const SequenceWorkItem& workItem)
{

    if(m_algorithm == ECA_KMER)
        return kmerCorrection(workItem);
    else
        return overlapCorrection(workItem);
}

ErrorCorrectResult ErrorCorrectProcess::overlapCorrection(const SequenceWorkItem& workItem)
{
    // Overlap based correction
    static const double p_error = 0.01f;
    bool done = false;
    int rounds = 0;
    
    ErrorCorrectResult result;
    SeqRecord currRead = workItem.read;
    std::string originalRead = workItem.read.seq.toString();

    while(!done)
    {
        m_blockList.clear();
        OverlapResult overlap_result = m_pOverlapper->overlapRead(currRead, m_minOverlap, &m_blockList);
        int sumOverlaps = 0;
        // Sum the spans of the overlap blocks to calculate the total number of overlaps this read has
        for(OverlapBlockList::iterator iter = m_blockList.begin(); iter != m_blockList.end(); ++iter)
        {
            assert(iter->ranges.interval[0].size() == iter->ranges.interval[1].size());
            sumOverlaps += iter->ranges.interval[0].size();
        }

        if(m_depthFilter > 0 && sumOverlaps > m_depthFilter)
        {
            result.num_prefix_overlaps = sumOverlaps;
            result.num_suffix_overlaps = sumOverlaps;
            result.correctSequence = currRead.seq;
            break;
        }

        // Convert the overlap block list into a multi-overlap 
        MultiOverlap mo = blockListToMultiOverlap(currRead, m_blockList);

        result.num_prefix_overlaps = 0;
        result.num_suffix_overlaps = 0;
        mo.countOverlaps(result.num_prefix_overlaps, result.num_suffix_overlaps);

        if(m_printOverlaps)
        {
            std::cout << "Prefix overlap: " << result.num_prefix_overlaps 
                      << " suffix overlaps: " << result.num_suffix_overlaps << "\n";
            std::cout << "Coverage overlap: " << mo.calculateCoverageOverlap() << "\n";
            mo.print();
        }

        if(m_algorithm == ECA_TRIE)
        {
            if(mo.isConflicted(m_conflictCutoff))
            {
                // Perform simple correction
                SeqTrie leftTrie;
                SeqTrie rightTrie;
                mo.makeSeqTries(p_error, leftTrie, rightTrie);
                result.correctSequence = ErrorCorrect::trieCorrect(currRead.seq.toString(), p_error, leftTrie, rightTrie);
            }
            else
            {
                result.correctSequence = mo.consensusConflict(p_error, m_conflictCutoff);
            }
        }
        else if(m_algorithm == ECA_CC)
        {
            result.correctSequence = mo.consensusConflict(p_error, m_conflictCutoff);
        }
        else if(m_algorithm == ECA_SIMPLE)
        {
            result.correctSequence = mo.calculateConsensusFromPartition(p_error);
        }
        else
        {
            assert(false);
        }

        ++rounds;
        if(rounds == m_numRounds || result.correctSequence == currRead.seq)
            done = true;
        else
            currRead.seq = result.correctSequence;
    }

    if(m_printOverlaps)
    {
        std::string corrected_seq = result.correctSequence.toString();
        std::cout << "OS: " << originalRead << "\n";
        std::cout << "CS: " << corrected_seq << "\n";
        std::cout << "DS: " << getDiffString(originalRead, corrected_seq) << "\n";
        std::cout << "QS: " << currRead.qual << "\n";
    }
    
    // Quality checks
    if(result.num_prefix_overlaps > 0 && result.num_suffix_overlaps > 0)
        result.passedQC = true;
    else
        result.passedQC = false;

    return result;
}

//
MultiOverlap ErrorCorrectProcess::blockListToMultiOverlap(const SeqRecord& record, OverlapBlockList& blockList)
{
    std::string read_idx = record.id;
    std::string read_seq = record.seq.toString();
    MultiOverlap out(read_idx, read_seq);

    for(OverlapBlockList::iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
    {
        std::string overlap_string = iter->getOverlapString(read_seq);

        // Compute the endpoints of the overlap
        int s1 = read_seq.length() - iter->overlapLen;
        int e1 = s1 + iter->overlapLen - 1;
        SeqCoord sc1(s1, e1, read_seq.length());

        int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
        int e2 = s2 + iter->overlapLen - 1;
        SeqCoord sc2(s2, e2, overlap_string.length());

        // The coordinates are always with respect to the read, so flip them if
        // we aligned to/from the reverse of the read
        if(iter->flags.isQueryRev())
            sc1.flip();
        if(iter->flags.isTargetRev())
            sc2.flip();

        bool isRC = false; // since we transformed the original sequence, they are never RC
        if(sc1.isContained())
            continue; // skip containments

        // Add an overlap for each member of the block
        for(int64_t i = iter->ranges.interval[0].lower; i <= iter->ranges.interval[0].upper; ++i)
        {
            Overlap o(read_idx, sc1, makeIdxString(i), sc2, isRC, -1);
            out.add(overlap_string, o);
        }
    }
    return out;
}

// make an id string from a read index
std::string ErrorCorrectProcess::makeIdxString(int64_t idx)
{
    std::stringstream ss;
    ss << idx;
    return ss.str();
}

// Correct a read with a k-mer based corrector
ErrorCorrectResult ErrorCorrectProcess::kmerCorrection(const SequenceWorkItem& workItem)
{
    ErrorCorrectResult result;
    SeqRecord currRead = workItem.read;
    std::string readSequence = workItem.read.seq.toString();

#ifdef KMER_TESTING
    std::cout << "Kmer correcting read " << workItem.read.id << "\n";
#endif

    size_t count_threshold = 3;
    size_t k_size = 27;
    size_t n = readSequence.size();
    size_t nk = n - k_size + 1;
    
    // Are all kmers in the read well-represented?
    bool allSolid = false;

    bool done = false;
    int rounds = 0;
    int maxAttempts = 4;
    while(!done)
    {
        // Compute the kmer counts across the read
        // and determine the positions in the read that are not covered by any solid kmers
        // These are the candidate incorrect bases
        std::vector<size_t> countVector(nk, 0);
        std::vector<size_t> solidVector(n, 0);

        for(size_t i = 0; i < nk; ++i)
        {
            std::string kmer = readSequence.substr(i, k_size);
            size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, m_pOverlapper->getBWT(), m_pOverlapper->getRBWT());
            countVector[i] = count;

            if(count > count_threshold)
            {
                for(size_t j = i; j < i + k_size; ++j)
                {
                    solidVector[j] = 1;
                }
            }
        }

        allSolid = true;
        for(size_t i = 0; i < n; ++i)
        {
#ifdef KMER_TESTING
            std::cout << "Position[" << i << "] = " << solidVector[i] << "\n";
#endif
            if(solidVector[i] != 1)
                allSolid = false;
        }
        
#ifdef KMER_TESTING  
        std::cout << "Read " << workItem.read.id << (allSolid ? " is solid\n" : " has potential errors\n");
#endif

        // Stop if all kmers are well represented or we have exceeded the number of correction rounds
        if(allSolid || rounds++ > maxAttempts)
            break;

        // Attempt to correct the leftmost potentially incorrect base
        bool corrected = false;
        for(size_t i = 0; i < n; ++i)
        {
            if(solidVector[i] != 1)
            {
                // Attempt to correct the base using the leftmost covering kmer
                size_t left_k_idx = (i + 1 >= k_size ? i + 1 - k_size : 0);
                corrected = attemptKmerCorrection(i, left_k_idx, k_size, std::max(countVector[left_k_idx], count_threshold), readSequence);
                if(corrected)
                    break;

                // base was not corrected, try using the rightmost covering kmer
                size_t right_k_idx = std::min(i, n - k_size);
                corrected = attemptKmerCorrection(i, right_k_idx, k_size, std::max(countVector[right_k_idx], count_threshold), readSequence);
                if(corrected)
                    break;
            }
        }

        // If no base in the read was corrected, stop the correction process
        if(!corrected)
        {
            assert(!allSolid);
            done = true;
        }
    }

    if(allSolid)
    {
        result.correctSequence = readSequence;
        result.passedQC = true;
    }
    else
    {
        result.passedQC = false;
    }
    return result;
}

// Attempt to correct the base at position idx in readSequence. Returns true if a correction was made
// The correction is made only if the count of the corrected kmer is at least minCount
bool ErrorCorrectProcess::attemptKmerCorrection(size_t i, size_t k_idx, size_t k_size, size_t minCount, std::string& readSequence)
{
    assert(i >= k_idx && i < k_idx + k_size);
    size_t base_idx = i - k_idx;
    char originalBase = readSequence[i];
    std::string kmer = readSequence.substr(k_idx, k_size);
    size_t bestCount = 0;
    char bestBase = '$';

#if KMER_TESTING
    std::cout << "i: " << i << " k-idx: " << k_idx << " " << kmer << " " << reverseComplement(kmer) << "\n";
#endif

    for(int j = 0; j < DNA_ALPHABET::size; ++j)
    {
        char currBase = ALPHABET[j];
        if(currBase == originalBase)
            continue;
        kmer[base_idx] = currBase;
        size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, m_pOverlapper->getBWT(), m_pOverlapper->getRBWT());

#if KMER_TESTING
        printf("%c %zu\n", currBase, count);
#endif
        if(count > bestCount && count >= minCount)
        {
            // Multiple corrections exist, do not correct
            if(bestBase != '$')
                return false;

            bestCount = count;
            bestBase = currBase;
        }
    }

    if(bestCount >= minCount)
    {
        assert(bestBase != '$');
        readSequence[i] = bestBase;
        return true;
    }
    return false;
}


//
//
//
ErrorCorrectPostProcess::ErrorCorrectPostProcess(std::ostream* pCorrectedWriter,
                                                 std::ostream* pDiscardWriter,
                                                 bool bCollectMetrics) : 
                                                      m_pCorrectedWriter(pCorrectedWriter),
                                                      m_pDiscardWriter(pDiscardWriter),
                                                      m_bCollectMetrics(bCollectMetrics),
                                                      m_totalBases(0), m_totalErrors(0),
                                                      m_readsKept(0), m_readsDiscarded(0)

{

}

//
void ErrorCorrectPostProcess::writeMetrics(std::ostream* pWriter)
{
    m_positionMetrics.write(pWriter, "Bases corrected by position\n", "pos");
    m_originalBaseMetrics.write(pWriter, "\nOriginal base that was corrected\n", "base");
    m_precedingSeqMetrics.write(pWriter, "\nkmer preceding the corrected base\n", "kmer");
    m_qualityMetrics.write(pWriter, "\nBases corrected by quality value\n\n", "quality");
        
    std::cout << "ErrorCorrect -- Corrected " << m_totalErrors << " out of " << m_totalBases <<
                 " bases (" << (double)m_totalErrors / m_totalBases << ")\n";
    std::cout << "Kept " << m_readsKept << " reads. Discarded " << m_readsDiscarded <<
                 " reads (" << (double)m_readsDiscarded / (m_readsKept + m_readsDiscarded)<< ")\n";
}

//
void ErrorCorrectPostProcess::process(const SequenceWorkItem& item, const ErrorCorrectResult& result)
{
    
    // Determine if the read should be discarded
    bool discardRead = !result.passedQC;

    // Collect metrics for the reads that were actually corrected
    if(m_bCollectMetrics && !discardRead)
    {
        collectMetrics(item.read.seq.toString(), 
                       result.correctSequence.toString(), 
                       item.read.qual);
    }

    SeqRecord record = item.read;
    record.seq = result.correctSequence;
    std::stringstream ss;
    ss << "PO:" << result.num_prefix_overlaps;
    ss << " SO:" << result.num_suffix_overlaps;

    if(!discardRead || m_pDiscardWriter == NULL)
    {
        record.write(*m_pCorrectedWriter, ss.str());
        ++m_readsKept;
    }
    else
    {
        record.write(*m_pDiscardWriter, ss.str());
        ++m_readsDiscarded;
    }
}

void ErrorCorrectPostProcess::collectMetrics(const std::string& originalSeq,
                                             const std::string& correctedSeq,
                                             const std::string& qualityStr)
{
    size_t precedingLen = 2;
    for(size_t i = 0; i < originalSeq.length(); ++i)
    {
        char qc = !qualityStr.empty() ? qualityStr[i] : '\0';
        char ob = originalSeq[i];

        ++m_totalBases;
        
        m_positionMetrics.incrementSample(i);

        if(!qualityStr.empty())
            m_qualityMetrics.incrementSample(qc);

        m_originalBaseMetrics.incrementSample(ob);

        std::string precedingMer;
        if(i > precedingLen)
        {
            precedingMer = originalSeq.substr(i - precedingLen, precedingLen);
            m_precedingSeqMetrics.incrementSample(precedingMer);
        }

        if(originalSeq[i] != correctedSeq[i])
        {
            m_positionMetrics.incrementError(i);
            if(!qualityStr.empty())
                m_qualityMetrics.incrementError(qc);
            m_originalBaseMetrics.incrementError(ob);

            if(!precedingMer.empty())
            {
                m_precedingSeqMetrics.incrementError(precedingMer);
            }
            ++m_totalErrors;
        }
    }
}
