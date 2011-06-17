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
#include "ReadCluster.h"

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

// Write the cluster results to a temporary output file
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
