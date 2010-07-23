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
                                            m_printOverlaps(printMO)
{

}

//
ErrorCorrectProcess::~ErrorCorrectProcess()
{

}

//
ErrorCorrectResult ErrorCorrectProcess::process(const SequenceWorkItem& workItem)
{
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

        // Convert the overlap block list into a multi-overlap 
        MultiOverlap mo = blockListToMultiOverlap(currRead, m_blockList);

        result.num_prefix_overlaps = 0;
        result.num_suffix_overlaps = 0;
        mo.countOverlaps(result.num_prefix_overlaps, result.num_suffix_overlaps);

        if(m_printOverlaps)
        {
            std::cout << "Prefix overlap: " << result.num_prefix_overlaps << " suffix overlaps: " << result.num_suffix_overlaps << "\n";
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

//
//
//
ErrorCorrectPostProcess::ErrorCorrectPostProcess(std::ostream* pWriter) : m_pWriter(pWriter)
{

}

//
void ErrorCorrectPostProcess::process(const SequenceWorkItem& item, const ErrorCorrectResult& result)
{
    SeqRecord correctedRecord = item.read;
    correctedRecord.seq = result.correctSequence;
    std::stringstream ss;
    ss << "PO:" << result.num_prefix_overlaps;
    ss << " SO:" << result.num_suffix_overlaps;

    correctedRecord.write(*m_pWriter, ss.str());
    //m_pOverlapper->writeResultASQG(*m_pASQGWriter, item.read, result);
}
