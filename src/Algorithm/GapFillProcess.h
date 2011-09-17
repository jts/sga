///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapFillProcess - Fill in gaps in a scaffold
//
#ifndef GAPFILL_PROCESS_H
#define GAPFILL_PROCESS_H

#include <list>
#include <stack>
#include <queue>
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "BitVector.h"
#include "HaplotypeBuilder.h"
#include "SequenceProcessFramework.h"
#include "BWTIntervalCache.h"
#include "SampledSuffixArray.h"

// Structures and typedefs

// Parameters structure
struct GapFillParameters
{
    // BWTS
    const BWT* pBWT; 
    const BWT* pRevBWT;

    const BWTIntervalCache* pBWTCache;
    const BWTIntervalCache* pRevBWTCache;
    
    size_t kmer;
    size_t kmerThreshold;

    int verbose;
};

enum GapFillReturnCode
{
    GFRC_OK,
    GFRC_NO_HAPLOTYPE,
    GFRC_NO_ANCHOR,
    GFRC_NUM_CODES
};

// Statistics tracking object
struct GapFillStats
{
    GapFillStats(); 
    size_t numGapsAttempted;
    size_t numGapsFilled;

    // Failure stats
    size_t numFails[GFRC_NUM_CODES];

    void print() const;
};

//
//
//
class GapFillProcess
{
    public:

        //
        // Functions
        //
        GapFillProcess(const GapFillParameters& params);
        ~GapFillProcess();
        
        // Generate haplotypes from chromosome refName, position [start, end]
        void processScaffold(const std::string& scaffold) const;

    private:
        
        //
        // Functions
        //

        GapFillReturnCode processGap(const std::string& scaffold, int gapStart, int gapEnd) const;
        AnchorSequence findAnchor(const std::string& scaffold, int64_t position, bool upstream) const;

        //
        // Data
        //
        GapFillParameters m_parameters;
        mutable GapFillStats m_stats;
};

#if 0
// 
class GapFillPostProcess
{

    public:
        GapFillPostProcess(const std::string& filename);
        ~GapFillPostProcess();

        //void process(const SequenceWorkItem& item, const GraphCompareResult& result);


    private:
        std::ostream* m_pWriter;
};
#endif

#endif
