///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MetAssembleProcess - Find contigs in an
// abstractly represented de Bruijn graph
// of a metagenome
//
#ifndef METASSEMBLE_PROCESS_H
#define METASSEMBLE_PROCESS_H

#include <list>
#include <stack>
#include <queue>
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "BitVector.h"
#include "VariationBubbleBuilder.h"
#include "SequenceProcessFramework.h"
#include "BWTIntervalCache.h"

// Parameters structure
struct MetAssembleParameters
{
    // BWTS
    const BWT* pBWT; 
    const BWT* pRevBWT;

    // FM-index
    const BWTIntervalCache* pBWTCache;
    const BWTIntervalCache* pRevBWTCache;
    
    size_t kmer;
    size_t kmerThreshold;
    size_t minLength;
    BitVector* pBitVector;
};

struct MetAssembleResult
{
    StringVector contigs;
};

//
//
//
class MetAssemble
{
    public:

        //
        // Functions
        //
        MetAssemble(const MetAssembleParameters& params);
        ~MetAssemble();
        
        // Process a read and all its kmers
        MetAssembleResult process(const SequenceWorkItem& item);

    private:
        
        //
        // Functions
        //

        // Perform an assemble starting at str
        std::string processKmer(const std::string& str, int count);

        // Mark all the kmers in str as being visited
        void markSequenceKmers(const std::string& str);

        // Break a contig into a lexicographically ordered set of kmers
        StringVector getLexicographicKmers(const std::string& contig) const;

        //
        // Data
        //
        MetAssembleParameters m_parameters;
};

// Shared result object that the threaded
// GraphCompare instances write to. Protected
// by a mutex
class MetAssembleAggregateResults
{

    public:
        MetAssembleAggregateResults(const std::string& filename);
        ~MetAssembleAggregateResults();

        void process(const SequenceWorkItem& item, const MetAssembleResult& result);


    private:

        std::ostream* m_pWriter;
        size_t m_numContigs;
        size_t m_basesWritten;
};

#endif
