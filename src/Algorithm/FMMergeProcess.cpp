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
#include "SGVisitors.h"

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
        if(m_pMarkedReads->test(i))
        {
            used = true;
            break;
        }
    }

    FMMergeResult result;
    //std::cout << "Processing read " << item.read.id << "\n";

    if(!used)
    {
        // Construct a new local graph around this read
        StringGraph* pGraph = new StringGraph;
        std::string rootID = "root";
        
        Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(rootID, readString);
        pGraph->addVertex(pVertex);

        // Add the root vertex to the result structure
        result.usedIntervals.push_back(readInterval);

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
            removeContainmentBlocks(currCandidate.pVertex->getSeqLen(), &candidateBlockList);

            bool validMergeNode = checkCandidate(currCandidate, &candidateBlockList);
            if(validMergeNode)
            {
                addCandidates(pGraph, currCandidate.pVertex, currCandidate.pEdge, &candidateBlockList, queue);
                result.usedIntervals.push_back(currCandidate.interval);
            }
            else
            {
                // Mark this vertex for later removal
                currCandidate.pVertex->setColor(GC_RED);
            }
        }
        
        // The graph has now been constructed. Remove all the nodes that are marked invalid for merging
        pGraph->sweepVertices(GC_RED);

        WARN_ONCE("Removing duplicate edge - consider doing this in the overlap step");
        SGDuplicateVisitor dupVisit(true);
        pGraph->visit(dupVisit);
        
        // Merge nodes
        pGraph->simplify();

        // If there was a cycle in the graph, it is possible that more than 1 vertex 
        // remains in the graph. Copy the vertex sequences into the result object.
        pGraph->getVertexSequences(result.mergedSequences);
        result.isMerged = true;
        delete pGraph;
    }
    else
    {
        result.isMerged = false;
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
            Overlap ovrXY = iter->toOverlap(pX->getID(), vertexID, pX->getSeqLen(), vertexSeq.length());

            // The vertex may already exist in the graph if the graph contains a loop
            Vertex* pY = pGraph->getVertex(vertexID);

            // Generate the new vertex
            if(pY == NULL)
            {
                pY = new(pGraph->getVertexAllocator()) Vertex(vertexID, vertexSeq);
                pGraph->addVertex(pY);
            }

            // Construct a description of the edge based on the overlap
            EdgeDesc ed = SGAlgorithms::overlapToEdgeDesc(pY, ovrXY);

            // If an edge with the same description as XY exists for X do not add a new edge or candidate
            if(!pX->hasEdge(ed))
            {
                // Construct the found edge and add it to the graph
                Edge* pXY = SGAlgorithms::createEdgesFromOverlap(pGraph, ovrXY, false);

                // Add the new vertex as a candidate
                FMMergeCandidate candidate;
                candidate.pVertex = pY;
                candidate.pEdge = pXY;
                candidate.interval = iter->getCanonicalInterval();
                candidateQueue.push(candidate);
            }
        }
    }
}

//
FMMergePostProcess::FMMergePostProcess(std::ostream* pWriter, BitVector* pMarkedReads) : m_numMerged(0), 
                                                                                          m_numTotal(0), 
                                                                                          m_totalLength(0), 
                                                                                          m_pWriter(pWriter), 
                                                                                          m_pMarkedReads(pMarkedReads)
{

}

//
FMMergePostProcess::~FMMergePostProcess()
{
    printf("[sga fm-merge] Merged %zu reads into %zu sequences\n", m_numTotal, m_numMerged);
    printf("[sga fm-merge] Reduction factor: %lf\n", (double)m_numTotal / m_numMerged);
    printf("[sga fm-merge] Mean merged size: %lf\n", (double)m_totalLength / m_numMerged);
}

//
void FMMergePostProcess::process(const SequenceWorkItem& item, const FMMergeResult& result)
{
    (void)item;
    (void)result;

    m_numTotal += 1;

    WARN_ONCE("check for collisions");
    if(result.isMerged)
    {
        // Write out the merged sequences
        for(std::vector<std::string>::const_iterator iter = result.mergedSequences.begin();
                iter != result.mergedSequences.end(); ++iter)
        {
            std::stringstream nameSS;
            nameSS << "merged-" << m_numMerged++;
            SeqRecord record;
            record.id = nameSS.str();
            record.seq = *iter;
            record.write(*m_pWriter);

            m_totalLength += iter->size();
        }

        // Set a bit mask for the indicated values
        for(std::vector<BWTInterval>::const_iterator iter = result.usedIntervals.begin();
               iter != result.usedIntervals.end(); ++iter)
        {
            for(int64_t i = iter->lower; i <= iter->upper; ++i)
                m_pMarkedReads->set(i, true);
        }
    }
}
