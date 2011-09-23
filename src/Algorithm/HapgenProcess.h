///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HapgenProcess - Generate candidate haplotypes
// from an assembly graph for a stream of sites
//
#ifndef HAPGEN_PROCESS_H
#define HAPGEN_PROCESS_H

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
struct HapgenParameters
{
    // BWTS
    const BWT* pBWT; 
    const BWT* pRevBWT;
    const SampledSuffixArray* pSSA;

    const BWTIntervalCache* pBWTCache;
    const BWTIntervalCache* pRevBWTCache;
 
    const ReadTable* pRefTable;
    size_t kmer;
    size_t kmerThreshold;

    int verbose;
};

//
//
//
class HapgenProcess
{
    public:

        //
        // Functions
        //
        HapgenProcess(const HapgenParameters& params);
        ~HapgenProcess();
        
        // Process a read and all its kmers
        //HapgenResult process(const SequenceWorkItem& item);

        // Generate haplotypes from chromosome refName, position [start, end]
        void processSite(const std::string& refName, size_t start, size_t end, const std::string& comment);
        AnchorSequence findAnchorKmer(const std::string& refName, int64_t start, bool upstream);

    private:
        
        //
        // Functions
        //

        // Extract the reads from the FM-index that share a kmer with any given haplotype
        void extractHaplotypeReads(const StringVector& haplotypes, bool doReverseComp, 
                                   SeqItemVector* pOutReads, SeqItemVector* pOutMates = NULL) const;

        //
        // Data
        //
        HapgenParameters m_parameters;
};

#if 0
// 
class HapgenPostProcess
{

    public:
        HapgenPostProcess(const std::string& filename);
        ~HapgenPostProcess();

        //void process(const SequenceWorkItem& item, const GraphCompareResult& result);


    private:
        std::ostream* m_pWriter;
};
#endif

#endif
