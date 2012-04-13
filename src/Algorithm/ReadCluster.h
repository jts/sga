///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadCluster - Generate a cluster of overlapping
// reads using the FM-index
//
#ifndef READCLUSTER_H
#define READCLUSTER_H

#include "Util.h"
#include "OverlapAlgorithm.h"

// A node in the queue of reads to use to expand the cluster
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
typedef std::vector<ClusterNode> ClusterNodeVector;
typedef std::set<int64_t> ClusterIntervalSet;

//
class ReadCluster
{
    public:
        ReadCluster(const OverlapAlgorithm* pOverlapper, int minOverlap);
        ClusterNode addSeed(const std::string& sequence, bool bCheckInIndex);

        // Run the cluster process. If the number of total nodes
        // exceeds max, abort the search.
        void run(size_t max_size, int max_iterations);

        ClusterNodeVector getOutput() const;
    
    private:

        ClusterNodeQueue m_queue;
        const OverlapAlgorithm* m_pOverlapper;
        int m_minOverlap;
        ClusterNodeVector m_outCluster;
        ClusterIntervalSet m_usedIndex;
};

#endif
