///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ClusterProcess - Compute clusters of overlapping reads
//
#ifndef CLUSTERPROCESS_H
#define CLUSTERPROCESS_H

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"
#include "BitVector.h"
#include "Bigraph.h"
#include "SGUtil.h"

//
struct ClusterNode
{
    std::string sequence;
    int64_t fwdCanonicalID;
    BWTInterval interval;
    bool isReverseInterval;

    static inline bool compare(const ClusterNode& a, const ClusterNode& b)
    {
        return BWTInterval::compare(a.interval, b.interval);
    }

    static inline bool equal(const ClusterNode& a, const ClusterNode& b)
    {
        return BWTInterval::equal(a.interval, b.interval);
    }
};
typedef std::queue<ClusterNode> ClusterNodeQueue;

struct ClusterResult
{
    std::vector<ClusterNode> clusterNodes;
};


// Compute the overlap blocks for reads
class ClusterProcess
{
    public:
        ClusterProcess(const OverlapAlgorithm* pOverlapper, 
                       int minOverlap, BitVector* pMarkedReads);

        ~ClusterProcess();

        ClusterResult process(const SequenceWorkItem& item);
    
    private:

        const OverlapAlgorithm* m_pOverlapper;
        const int m_minOverlap;
        BitVector* m_pMarkedReads;
};

// Write the results from the overlap step to an ASQG file
class ClusterPostProcess
{
    public:
        ClusterPostProcess(std::ostream* pWriter, size_t minClusterSize, BitVector* pMarkedReads);
        ~ClusterPostProcess();
        
        void process(const SequenceWorkItem& item, const ClusterResult& result);

    private:
        size_t m_minClusterSize;
        size_t m_numClusters;
        size_t m_numTotalReads;
        size_t m_numTotalReadsClustered;
        
        std::ostream* m_pWriter;
        BitVector* m_pMarkedReads;
};

#endif
