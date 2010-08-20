///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringGraphGenerator - Iteratively construct 
// a local String Graph using the FM-index
//
#include "StringGraphGenerator.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"

StringGraphGenerator::StringGraphGenerator(const OverlapAlgorithm* pOverlapper, 
                                           const SeqRecord& startRead, 
                                           const SeqRecord& endRead, 
                                           int minOverlap,
                                           EdgeDir startDir) : m_pOverlapper(pOverlapper), m_minOverlap(minOverlap), m_pGraph(NULL)
{
    m_pGraph = new StringGraph;

    // Add the start and end vertices to the graph
    m_pStartVertex = addTerminalVertex(startRead);
    m_pEndVertex = addTerminalVertex(endRead);

    // Set the color of the starting node to be UNEXPLORED
    // and the color of the end node to be EXPLORED.
    // This indicates to the subsequent search which vertices
    // should be expanded
    m_pStartVertex->setColor(UNEXPLORED_COLOR);
    m_pEndVertex->setColor(EXPLORED_COLOR);

    // Set up the expansion frontier
    FrontierQueue queue;
    GraphFrontier node;
    node.pVertex = m_pStartVertex;
    node.dir = startDir;
    node.distance = m_pStartVertex->getSeqLen();
    queue.push(node);

    buildGraph(queue, 300);

    m_pGraph->writeDot("local.dot");

    SGDuplicateVisitor dupVisit;
    m_pGraph->visit(dupVisit);

    SGContainRemoveVisitor containVisit;
    m_pGraph->visit(containVisit);
    m_pGraph->writeDot("local-final.dot");

}

//
StringGraphGenerator::~StringGraphGenerator()
{
    delete m_pGraph;
    m_pGraph = NULL;

    // m_pStartVertex and m_pEndVertex are deleted by the graph destructor
    // so they do not need to be explicitly freed here.
}

// Build the graph by expanding nodes on the frontier
void StringGraphGenerator::buildGraph(FrontierQueue& queue, int maxDistance)
{
    while(!queue.empty())
    {
        GraphFrontier node = queue.front();
        queue.pop();
        if(node.pVertex->getColor() == EXPLORED_COLOR)
            continue; // node has been visited already
        
        std::cout << "expanding " << node.pVertex->getID() << "\n";

        // Search the FM-index for the current vertex
        SeqRecord record;
        record.id = node.pVertex->getID();
        record.seq = node.pVertex->getSeq().toString();
        
        OverlapBlockList blockList;
        assert(blockList.empty());
        m_pOverlapper->overlapRead(record, m_minOverlap, &blockList);

        // Update the graph and the frontier queue with newly found vertices
        updateGraphAndQueue(node, queue, blockList, maxDistance);
        node.pVertex->setColor(EXPLORED_COLOR);
    }

    m_pGraph->setColors(GC_WHITE);
}

// 
void StringGraphGenerator::updateGraphAndQueue(GraphFrontier& currNode, FrontierQueue& queue, OverlapBlockList& blockList, int maxDistance)
{
    // Partition the block list into containment blocks and extension (valid) blocks
    // We do not add containment edges to the graph so the containments are discarded
    OverlapBlockList containList;
    OverlapBlockList overlapList;

    Vertex* pX = currNode.pVertex;

    //partitionBlockList(pX->getSeqLen(), &blockList, &overlapList, &containList);

    // Process the overlap blocks, adding new vertices and edges where necessary
    for(OverlapBlockList::iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
    {
        if(iter->getEdgeDir() != currNode.dir)
            continue;

        std::string vertexID = overlapBlockToCanonicalID(*iter);
        if(vertexID == pX->getID())
            continue; // skip self-edges


        std::string vertexSeq = iter->getFullString(pX->getSeq().toString());
        Overlap o = iter->toOverlap(pX->getID(), vertexID, pX->getSeqLen(), vertexSeq.length());
        


        std::cout << "has overlap to: " << vertexID << " len: " << iter->overlapLen << " flags: " << iter->flags << "\n";
        std::cout << "Overlap string: " << iter->getOverlapString(pX->getSeq().toString()) << "\n";

        // Check if a vertex with endVertexID exists in the graph
        Vertex* pVertex = m_pGraph->getVertex(vertexID);
        if(pVertex == NULL)
        {
            std::cout << "Vertex with ID: " << vertexID << " does not exist\n";

            // Generate the new vertex
            vertexSeq = iter->getFullString(pX->getSeq().toString());
            std::cout << "CREATING VERTEX: " << vertexID << "\n";
            pVertex = new Vertex(vertexID, vertexSeq);
            pVertex->setColor(UNEXPLORED_COLOR);
            m_pGraph->addVertex(pVertex);
        }

        // Construct the found edge
        std::cout << "Vertex sequence: " << vertexSeq << "\n";
        Edge* pXY = SGAlgorithms::createEdgesFromOverlap(m_pGraph, o, true);
        std::cout << "CREATED EDGE DIR: " << pXY->getDir() << "\n";
        // If the endpoint vertex is unexplored, queue it
        if(pVertex->getColor() == UNEXPLORED_COLOR)
        {
            GraphFrontier node;
            node.pVertex = pVertex;
            node.dir = !pXY->getTwin()->getDir(); // continuation direction
            node.distance = currNode.distance + pXY->getSeqLen();
            if(node.distance <= maxDistance)
                queue.push(node);
        }
    }
}

Vertex* StringGraphGenerator::addTerminalVertex(const SeqRecord& record)
{
    assert(m_pGraph != NULL);

    // Build the vertex by performing a full-length search for the
    // sequence in the FM-index. We set the ID of the vertex to be the 
    // lowest index in the returned block list.
    OverlapBlockList endBlockList;
    m_pOverlapper->alignReadDuplicate(record, &endBlockList);

    // Search the block list for the exact match to the end read. This must exist
    OverlapBlockList::iterator matchIter = endBlockList.begin();
    while(matchIter != endBlockList.end())
    {
        if(matchIter->numDiff == 0 && !matchIter->flags.isQueryRev())
            break; // this block corresponds to the actual sequence of endRead
    }
    assert(matchIter != endBlockList.end());
    
    // Construct the canonical ID from the matching interval
    std::string endID = overlapBlockToCanonicalID(*matchIter);

    Vertex* pVertex = new Vertex(endID, record.seq.toString());
    m_pGraph->addVertex(pVertex);
    return pVertex;
}

//
std::string StringGraphGenerator::overlapBlockToCanonicalID(OverlapBlock& block)
{
    std::stringstream ss;
    int ci = block.getCanonicalIntervalIndex();
    std::cout << "CI: " << ci << "\n";
    std::cout << "ID[0]: " << block.ranges.interval[0].lower << "\n";
    std::cout << "ID[1]: " << block.ranges.interval[1].lower << "\n";
    ss << "IDX-" << block.ranges.interval[ci].lower;
    return ss.str();
}

