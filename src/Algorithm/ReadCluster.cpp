///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadCluster - Generate a cluster of overlapping
// reads using the FM-index
//
#include "ReadCluster.h"

//
ReadCluster::ReadCluster(const OverlapAlgorithm* pOverlapper, int minOverlap) : m_pOverlapper(pOverlapper), 
                                                                                m_minOverlap(minOverlap)
{

}

// Add a seed read to the cluster. Overlaps will be found for 
// each seed read to grow the cluster
ClusterNode ReadCluster::addSeed(const std::string& sequence)
{
    // Check if this read is a substring
    SeqRecord tempRecord;
    tempRecord.id = "cluster-seed";
    tempRecord.seq = sequence;

    OverlapBlockList tempBlockList;
    OverlapResult overlapResult = m_pOverlapper->alignReadDuplicate(tempRecord, &tempBlockList);
    if(overlapResult.isSubstring)
    {
        std::cerr << "Error: substring reads found in sga-cluster. Please run sga filter (with the --no-kmer-check flag) before cluster\n";
        exit(1);
    }

    // Find the interval in the fm-index containing the read
    const BWT* pBWT = m_pOverlapper->getBWT();
    BWTInterval readInterval = BWTAlgorithms::findInterval(pBWT, sequence);
    BWTAlgorithms::updateInterval(readInterval, '$', pBWT);

    // The read must be present in the index
    assert(readInterval.isValid());

    ClusterNode node;
    node.sequence = sequence;
    node.interval = readInterval;
    node.isReverseInterval = false;
    m_usedIndex.insert(readInterval.lower);
    m_queue.push(node);
    return node;
}

// Run the cluster process. If the number of total nodes
// exceeds max, abort the search.
void ReadCluster::run(size_t max)
{
    while(!m_queue.empty())
    {
        if(m_queue.size() + m_outCluster.size() > max)
        {
            while(!m_queue.empty())
                m_queue.pop();
            m_outCluster.clear();
            return;
        }

        ClusterNode node = m_queue.front();
        m_queue.pop();

        // Add this node to the output
        m_outCluster.push_back(node);

        // Find overlaps for the current node
        SeqRecord tempRecord;
        tempRecord.id = "cluster";
        tempRecord.seq = node.sequence;
        OverlapBlockList blockList;
        m_pOverlapper->overlapRead(tempRecord, m_minOverlap, &blockList);
        
        // Parse each member of the block list and potentially expand the cluster
        for(OverlapBlockList::const_iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
        {
            // Check if the reads in this block are part of the cluster already
            BWTInterval canonicalInterval = iter->getCanonicalInterval();
            int64_t canonicalIndex = canonicalInterval.lower;
            if(m_usedIndex.count(canonicalIndex) == 0)
            {
                // This is a new node that isn't in the cluster. Add it.
                m_usedIndex.insert(canonicalIndex);

                ClusterNode newNode;
                newNode.sequence = iter->getFullString(node.sequence);
                newNode.interval = canonicalInterval;
                newNode.isReverseInterval = iter->flags.isTargetRev();
                m_queue.push(newNode);
            }
        }
    }
}

// 
ClusterNodeVector ReadCluster::getOutput() const
{
    // Sort the intervals into ascending order and remove any duplicate intervals (which can occur
    // if the subgraph has a simple cycle)
    ClusterNodeVector retVector = m_outCluster;
    std::sort(retVector.begin(), retVector.end(), ClusterNode::compare);

    std::vector<ClusterNode>::iterator newEnd = std::unique(retVector.begin(),
                                                            retVector.end(),
                                                            ClusterNode::equal);

    retVector.erase(newEnd, retVector.end());
    return retVector;
}
