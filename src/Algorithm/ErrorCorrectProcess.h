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

enum ErrorCorrectAlgorithm
{
    ECA_TRIE, // aggressive trie-based correction of conflicted sequences
    ECA_CC, // conflict-aware consensus
    ECA_SIMPLE // straightforward correct
};

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "MultiOverlap.h"

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
};

// Write the results from the overlap step to an ASQG file
class ErrorCorrectPostProcess
{
    public:
        ErrorCorrectPostProcess(std::ostream* pCorrectedWriter, std::ostream* pDiscardWriter);
        void process(const SequenceWorkItem& item, const ErrorCorrectResult& result);

    private:
        std::ostream* m_pCorrectedWriter;
        std::ostream* m_pDiscardWriter;
};

#endif
