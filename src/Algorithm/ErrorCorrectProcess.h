///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ErrorCorrectProcess - Wrapper to perform error correction
// for a sequence work item
//
#ifndef CORRECTROCESS_H
#define CORRECTPROCESS_H

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "MultiOverlap.h"
#include "Metrics.h"

enum ErrorCorrectAlgorithm
{
    ECA_TRIE, // aggressive trie-based correction of conflicted sequences
    ECA_CC, // conflict-aware consensus
    ECA_SIMPLE, // straightforward correct
    ECA_KMER // K-mer based correction algorithm
};

enum ECFlag
{
    ECF_NOTCORRECTED,
    ECF_CORRECTED,
    ECF_AMBIGIOUS,
    ECF_DUPLICATE
};

class ErrorCorrectResult
{
    public:
        DNAString correctSequence;
        ECFlag flag;

        // Metrics
        bool passedQC;
        size_t num_prefix_overlaps;
        size_t num_suffix_overlaps;
};

//
class ErrorCorrectProcess
{
    public:
        ErrorCorrectProcess(const OverlapAlgorithm* pOverlapper, 
                            int minOverlap, int numRounds, 
                            int conflictCutoff, ErrorCorrectAlgorithm algo,
                            bool printMO);

        ~ErrorCorrectProcess();

        ErrorCorrectResult process(const SequenceWorkItem& item);
    
    private:
        
        ErrorCorrectResult kmerCorrection(const SequenceWorkItem& item);
        ErrorCorrectResult overlapCorrection(const SequenceWorkItem& workItem);

        bool attemptKmerCorrection(size_t i, size_t k_idx, size_t k_size, size_t minCount, std::string& readSequence);

        MultiOverlap blockListToMultiOverlap(const SeqRecord& record, 
                                             OverlapBlockList& blockList);

        std::string makeIdxString(int64_t idx);

        OverlapBlockList m_blockList;
        const OverlapAlgorithm* m_pOverlapper;
        const int m_minOverlap;
        const int m_numRounds;
        const int m_conflictCutoff;
        const ErrorCorrectAlgorithm m_algorithm;
        const bool m_printOverlaps;
        const int m_depthFilter;
};

// Write the results from the overlap step to an ASQG file
class ErrorCorrectPostProcess
{
    public:
        ErrorCorrectPostProcess(std::ostream* pCorrectedWriter, 
                                std::ostream* pDiscardWriter, bool bCollectMetrics);

        void process(const SequenceWorkItem& item, const ErrorCorrectResult& result);
        void writeMetrics(std::ostream* pWriter);

    private:

        void collectMetrics(const std::string& originalSeq, 
                            const std::string& correctedSeq, const std::string& qualityStr);

        std::ostream* m_pCorrectedWriter;
        std::ostream* m_pDiscardWriter;
        bool m_bCollectMetrics;

        ErrorCountMap<char> m_qualityMetrics;
        ErrorCountMap<int> m_positionMetrics;
        ErrorCountMap<char> m_originalBaseMetrics;
        ErrorCountMap<std::string> m_precedingSeqMetrics;

        size_t m_totalBases;
        size_t m_totalErrors;
        size_t m_readsKept;
        size_t m_readsDiscarded;

};

#endif
