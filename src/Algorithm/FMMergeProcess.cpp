///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// FMMergeProcess - Merge unambiguously overlapping sequences
//
#include "FMMergeProcess.h"
#include "SGAlgorithms.h"

//
FMMergeProcess::FMMergeProcess(const OverlapAlgorithm* pOverlapper, int minOverlap, const BitVector* pMarkedReads) : 
                                     m_pOverlapper(pOverlapper), 
                                     m_minOverlap(minOverlap), 
                                     m_pMarkedReads(pMarkedReads)
{

}

//
FMMergeProcess::~FMMergeProcess()
{


}

//
FMMergeResult FMMergeProcess::process(const SequenceWorkItem& item)
{
    (void)item;

    // Calculate the intervals in the forward FM-index for this read
    const BWT* pBWT = m_pOverlapper->getBWT();

    // Find the interval in the fm-index containing the read
    std::string readString = item.read.seq.toString();
    BWTInterval readInterval = BWTAlgorithms::findInterval(pBWT, readString);

    // Update the interval by looking for the $ to map the interval into indices in the bit array
    BWTAlgorithms::updateInterval(readInterval, '$', pBWT);

    // The read must be present in the index
    assert(readInterval.isValid());

    // Check if this read has been used yet
    bool used = false;
    for(int64_t i = readInterval.lower; i <= readInterval.upper; ++i)
    {
        used = m_pMarkedReads->test(i) || used;
    }

    FMMergeResult result;

    if(!used)
    {
        // Construct a new local graph around this read
        StringGraph* pGraph = new StringGraph;
        std::string rootID = "root";
        
        Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(rootID, readString);
        pGraph->addVertex(pVertex);

        // Add the root vertex to the result structure
        result.usedIntervals.push_back(readInterval);
        result.usedSequences.push_back(readString);

        // Enqueue the read for overlap detection in both directions
        FMMergeQueue queue;

        // Construct the overlap block list for this node
        SeqRecord record;
        record.id = pVertex->getID();
        record.seq = pVertex->getSeq().toString();
        OverlapBlockList blockList;
        m_pOverlapper->overlapRead(record, m_minOverlap, &blockList);
        removeContainmentBlocks(pVertex->getSeqLen(), &blockList);

        // Construct the initial list of candidates from the block list
        addCandidates(pGraph, pVertex, NULL, &blockList, queue);

        while(!queue.empty())
        {
            FMMergeCandidate currCandidate = queue.front();
            queue.pop();

            // Determine whether this is a valid vertex to merge or not.
            // It is valid if it has a single edge in the direction of the vertex
            // that added it to the candidate list
            SeqRecord record;
            record.id = currCandidate.pVertex->getID();
            record.seq = currCandidate.pVertex->getSeq().toString();

            OverlapBlockList candidateBlockList;
            m_pOverlapper->overlapRead(record, m_minOverlap, &candidateBlockList);
            removeContainmentBlocks(pVertex->getSeqLen(), &candidateBlockList);

            bool validMergeNode = checkCandidate(currCandidate, &candidateBlockList);
            if(validMergeNode)
            {
                addCandidates(pGraph, currCandidate.pVertex, currCandidate.pEdge, &candidateBlockList, queue);
                result.usedIntervals.push_back(currCandidate.interval);
                result.usedSequences.push_back(currCandidate.pVertex->getSeq().toString());
            }
            else
            {
                // Mark this vertex for later removal
                currCandidate.pVertex->setColor(GC_RED);
            }
        }
        
        pGraph->sweepVertices(GC_RED);
        delete pGraph;
    }

    return result;
}

// Check if the candidate node can be merged with the node it is linked to. Returns true if so
bool FMMergeProcess::checkCandidate(const FMMergeCandidate& candidate, const OverlapBlockList* pBlockList) const
{
    // Get the direction of the edge back to the node that generated this candidate
    // pEdge is the edge TO the candidate vertex so the twin direction is the direction away
    // from the candidate vertex
    EdgeDir mergeDir = candidate.pEdge->getTwinDir();

    size_t mergeDirCount = 0;
    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter)
    {
        if(iter->getEdgeDir() == mergeDir)
            mergeDirCount += 1;
    }

    assert(mergeDirCount > 0);
    return mergeDirCount == 1;
}

// Add the edges starting from pX as candidate vertices
// using the blockList.
// Precondition: pX is a valid vertex in the merge graph. In other words, there is a
// unique assembly that includes pX and the root vertex.
void FMMergeProcess::addCandidates(StringGraph* pGraph, const Vertex* pX, const Edge* pEdgeToX, const OverlapBlockList* pBlockList, FMMergeQueue& candidateQueue)
{
    // Count the number of edges in each direction. 
    size_t numAnti = 0;
    size_t numSense = 0;

    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter)
    {
        if(iter->getEdgeDir() == ED_SENSE)
            ++numSense;
        if(iter->getEdgeDir() == ED_ANTISENSE)
            ++numAnti;
    }

    std::cout << "Num AS: " << numAnti << "\n";
    std::cout << "Num S: " << numSense << "\n";

    // For each edge block, if it is unique for the direction add the vertex it describes as a candidate
    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter)
    {
        EdgeDir currDir = iter->getEdgeDir();
        if((currDir == ED_SENSE && numSense == 1) || 
           (currDir == ED_ANTISENSE && numAnti == 1))
        {
            // Skip edges in the direction to X
            if(pEdgeToX != NULL && pEdgeToX->getTwinDir() == currDir)
                continue;

            // Construct new candidate vertices and add them to the graph
            std::string vertexID = iter->toCanonicalID();
            assert(vertexID != pX->getID());
            std::string vertexSeq = iter->getFullString(pX->getSeq().toString());
            Overlap o = iter->toOverlap(pX->getID(), vertexID, pX->getSeqLen(), vertexSeq.length());

            // Ensure this vertex does not exist in the graph
            Vertex* pY = pGraph->getVertex(vertexID);
            assert(pY == NULL);

            // Generate the new vertex
            pY = new(pGraph->getVertexAllocator()) Vertex(vertexID, vertexSeq);
            pGraph->addVertex(pY);

            // Construct the found edge and add it to the graph
            Edge* pXY = SGAlgorithms::createEdgesFromOverlap(pGraph, o, false);

            // Add the new vertex as a candidate
            FMMergeCandidate candidate;
            candidate.pVertex = pY;
            candidate.pEdge = pXY;
            candidate.interval = iter->getCanonicalInterval();
            candidateQueue.push(candidate);
        }
    }
}

//
FMMergePostProcess::FMMergePostProcess(std::ostream* pWriter, BitVector* pMarkedReads) : m_pWriter(pWriter), m_pMarkedReads(pMarkedReads)
{

}

//
void FMMergePostProcess::process(const SequenceWorkItem& item, const FMMergeResult& result)
{
    (void)item;
    (void)result;

    std::cout << "Assembled strings:\n";
    std::copy(result.usedSequences.begin(), result.usedSequences.end(), std::ostream_iterator<std::string>(std::cout, "\n"));

    // Set a bit mask for the indicated values
    for(std::vector<BWTInterval>::const_iterator iter = result.usedIntervals.begin();
           iter != result.usedIntervals.end(); ++iter)
    {
        for(int64_t i = iter->lower; i <= iter->upper; ++i)
            m_pMarkedReads->set(i, true);
    }

}
