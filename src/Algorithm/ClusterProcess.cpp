///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ClusterProcess - Merge unambiguously overlapping sequences
//
#include "ClusterProcess.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"

//
ClusterProcess::ClusterProcess(const OverlapAlgorithm* pOverlapper, int minOverlap, BitVector* pMarkedReads) : 
                                     m_pOverlapper(pOverlapper), 
                                     m_minOverlap(minOverlap), 
                                     m_pMarkedReads(pMarkedReads)
{

}

//
ClusterProcess::~ClusterProcess()
{


}

//
ClusterResult ClusterProcess::process(const SequenceWorkItem& item)
{
    // Calculate the intervals in the forward FM-index for this read
    const BWT* pBWT = m_pOverlapper->getBWT();

    // Find the interval in the fm-index containing the read
    std::string readString = item.read.seq.toString();
    BWTInterval readInterval = BWTAlgorithms::findInterval(pBWT, readString);
    BWTAlgorithms::updateInterval(readInterval, '$', pBWT);

    // The read must be present in the index
    assert(readInterval.isValid());

    // Check if this read has been used yet
    bool used = false;
    for(int64_t i = readInterval.lower; i <= readInterval.upper; ++i)
    {
        if(m_pMarkedReads->test(i))
        {
            used = true;
            break;
        }
    }

    ClusterResult result;
    if(used)
        return result; // already part of a cluster, return nothing

    // Compute a new cluster around this read
    std::set<int64_t> usedIndex;
    ClusterNodeQueue queue;
    ClusterNode node;
    node.sequence = item.read.seq.toString();
    node.interval = readInterval;
    node.isReverseInterval = false;
    usedIndex.insert(readInterval.lower);
    queue.push(node);

    while(!queue.empty())
    {
        ClusterNode node = queue.front();
        queue.pop();

        // Update the used index and the result structure with this node's data
        result.clusterNodes.push_back(node);

        SeqRecord tempRecord;
        tempRecord.id = "cluster";
        tempRecord.seq = node.sequence;

        OverlapBlockList blockList;
        OverlapResult result = m_pOverlapper->overlapRead(tempRecord, m_minOverlap, &blockList);

        // Parse each member of the block list and potentially expand the cluster
        for(OverlapBlockList::const_iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
        {
            // Check if the reads in this block are part of the cluster already
            BWTInterval canonicalInterval = iter->getCanonicalInterval();
            int64_t canonicalIndex = canonicalInterval.lower;
            if(usedIndex.count(canonicalIndex) == 0)
            {
                usedIndex.insert(canonicalIndex);
                ClusterNode newNode;
                newNode.sequence = iter->getFullString(node.sequence);
                newNode.interval = canonicalInterval;
                newNode.isReverseInterval = iter->flags.isTargetRev();
                queue.push(newNode);
            }
        }
    }

    // If some work was performed, update the bitvector so other threads do not try to merge the same set of reads.
    // This uses compare-and-swap instructions to ensure the uppdate is atomic. 
    // If some other thread has merged this set (and updated
    // the bitvector), we discard all the merged data.
    
    // As a given set of reads should all be merged together, we only need to make sure we atomically update
    // the bit for the read with the lowest index in the set.

    // Sort the intervals into ascending order and remove any duplicate intervals (which can occur
    // if the subgraph has a simple cycle)
    std::sort(result.clusterNodes.begin(), result.clusterNodes.end(), ClusterNode::compare);
    std::vector<ClusterNode>::iterator newEnd = std::unique(result.clusterNodes.begin(),
                                                            result.clusterNodes.end(),
                                                            ClusterNode::equal);

    size_t oldSize = result.clusterNodes.size();
    result.clusterNodes.erase(newEnd, result.clusterNodes.end());
    size_t newSize = result.clusterNodes.size();
    if(oldSize != newSize)
        std::cout << "Warning: duplicate cluster nodes were found\n";

    // Check if the bit in the vector has already been set for the lowest read index
    // If it has some other thread has already output this set so we do nothing
    int64_t lowestIndex = result.clusterNodes.front().interval.lower;
    bool currentValue = m_pMarkedReads->test(lowestIndex);
    bool updateSuccess = false;

    if(currentValue == false)
    {
        // Attempt to update the bit vector with an atomic CAS. If this returns false
        // the bit was set by some other thread
        updateSuccess = m_pMarkedReads->updateCAS(lowestIndex, currentValue, true);
    }

    if(updateSuccess)
    {
        // We successfully atomically set the bit for the first read in this set
        // to true. We can safely update the rest of the bits and keep the merged sequences
        // for output.
        std::vector<ClusterNode>::const_iterator iter = result.clusterNodes.begin();
        for(; iter != result.clusterNodes.end(); ++iter)
        {
            for(int64_t i = iter->interval.lower; i <= iter->interval.upper; ++i)
            {
                if(i == lowestIndex) //already set
                    continue;

                currentValue = m_pMarkedReads->test(i);
                if(currentValue)
                {
                    // This value should not be true, emit a warning
                    std::cout << "Warning: Bit " << i << " was set outside of critical section\n";
                }
                else
                {
                    m_pMarkedReads->updateCAS(i, currentValue, true);
                }
            }
        }
    }
    else
    {
        // Some other thread merged these reads already, discard the intermediate
        // data and set the result to false
        result.clusterNodes.clear();
    }
    return result;
}

//
ClusterPostProcess::ClusterPostProcess(std::ostream* pWriter, 
                                       size_t minClusterSize, 
                                       BitVector* pMarkedReads) : m_minClusterSize(minClusterSize),
                                                                  m_numClusters(0), 
                                                                  m_numTotalReads(0), 
                                                                  m_numTotalReadsClustered(0), 
                                                                  m_pWriter(pWriter), 
                                                                  m_pMarkedReads(pMarkedReads)
{

}

//
ClusterPostProcess::~ClusterPostProcess()
{
    printf("[sga cluster] Clustered %zu reads into %zu clusters (%zu total reads input)\n", m_numTotalReadsClustered, m_numClusters, m_numTotalReads);
}

//
void ClusterPostProcess::process(const SequenceWorkItem& item, const ClusterResult& result)
{
    (void)item;
    m_numTotalReads += 1;

    if(result.clusterNodes.size() > 0)
    {
        // Compute the cluster size
        size_t clusterSize = 0;
        for(size_t i = 0; i < result.clusterNodes.size(); ++i)
            clusterSize += result.clusterNodes[i].interval.size();
        
        if(clusterSize >= m_minClusterSize)
        {
            for(size_t i = 0; i < result.clusterNodes.size(); ++i)
                *m_pWriter << "cluster-" << m_numClusters << "\t" << clusterSize << "\t" << result.clusterNodes[i].sequence << "\t" << result.clusterNodes[i].interval << "\n";
            m_numClusters += 1;
            m_numTotalReadsClustered += clusterSize;
        }
    }
}
