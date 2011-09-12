//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGVisitors - Algorithms that visit
// each vertex in the graph and perform some
// operation
//
#include "SGVisitors.h"
#include "ErrorCorrect.h"
#include "CompleteOverlapSet.h"
#include "SGSearch.h"
#include "stdaln.h"

//
// SGFastaVisitor - output the vertices in the graph in 
// fasta format
//
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    m_fileHandle << ">" << pVertex->getID() << " " <<  pVertex->getSeq().length() 
                 << " " << 0 << "\n";
    m_fileHandle << pVertex->getSeq() << "\n";
    return false;
}


//
// SGOverlapWriterVisitor - write all the overlaps in the graph to a file 
//
bool SGOverlapWriterVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    EdgePtrVec edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Overlap ovr = edges[i]->getOverlap();
        if(ovr.id[0] < ovr.id[1])
            m_fileHandle << ovr << "\n";
    }
    return false;
}

//
// SGTransRedVisitor - Perform a transitive reduction about this vertex
// This uses Myers' algorithm (2005, The fragment assembly string graph)
// Precondition: the edge list is sorted by length (ascending)
void SGTransitiveReductionVisitor::previsit(StringGraph* pGraph)
{
    // The graph must not have containments
    assert(!pGraph->hasContainment());

    // Set all the vertices in the graph to "vacant"
    pGraph->setColors(GC_WHITE);
    pGraph->sortVertexAdjListsByLen();

    marked_verts = 0;
    marked_edges = 0;
}

bool SGTransitiveReductionVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    size_t trans_count = 0;
    static const size_t FUZZ = 10; // see myers

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pVertex->getEdges(dir); // These edges are already sorted
        if(edges.size() == 0)
            continue;

        for(size_t i = 0; i < edges.size(); ++i)
            (edges[i])->getEnd()->setColor(GC_GRAY);

        Edge* pLongestEdge = edges.back();
        size_t longestLen = pLongestEdge->getSeqLen() + FUZZ;
        
        // Stage 1
        for(size_t i = 0; i < edges.size(); ++i)
        {
            Edge* pVWEdge = edges[i];
            Vertex* pWVert = pVWEdge->getEnd();

            //std::cout << "Examining edges from " << pWVert->getID() << " longest: " << longestLen << "\n";
            //std::cout << pWVert->getID() << " w_edges: \n";
            EdgeDir transDir = !pVWEdge->getTwinDir();
            if(pWVert->getColor() == GC_GRAY)
            {
                EdgePtrVec w_edges = pWVert->getEdges(transDir);
                for(size_t j = 0; j < w_edges.size(); ++j)
                {
                    Edge* pWXEdge = w_edges[j];
                    size_t trans_len = pVWEdge->getSeqLen() + pWXEdge->getSeqLen();
                    if(trans_len <= longestLen)
                    {
                        if(pWXEdge->getEnd()->getColor() == GC_GRAY)
                        {
                            // X is the endpoint of an edge of V, therefore it is transitive
                            pWXEdge->getEnd()->setColor(GC_BLACK);
                            //std::cout << "Marking " << pWXEdge->getEndID() << " as transitive to " << pVertex->getID() << "\n";
                        }
                    }
                    else
                        break;
                }
            }
        }
        
        // Stage 2
        for(size_t i = 0; i < edges.size(); ++i)
        {
            Edge* pVWEdge = edges[i];
            Vertex* pWVert = pVWEdge->getEnd();

            //std::cout << "Examining edges from " << pWVert->getID() << " longest: " << longestLen << "\n";
            //std::cout << pWVert->getID() << " w_edges: \n";
            EdgeDir transDir = !pVWEdge->getTwinDir();
            EdgePtrVec w_edges = pWVert->getEdges(transDir);
            for(size_t j = 0; j < w_edges.size(); ++j)
            {
                //std::cout << "    edge: " << *w_edges[j] << "\n";
                Edge* pWXEdge = w_edges[j];
                size_t len = pWXEdge->getSeqLen();

                if(len < FUZZ || j == 0)
                {
                    if(pWXEdge->getEnd()->getColor() == GC_GRAY)
                    {
                        // X is the endpoint of an edge of V, therefore it is transitive
                        pWXEdge->getEnd()->setColor(GC_BLACK);
                        //std::cout << "Marking " << pWXEdge->getEndID() << " as transitive to " << pVertex->getID() << " in stage 2";
                        //std::cout << " via " << pWVert->getID() << "\n";
                    }
                }
                else
                {
                    break;
                }
            }
        }

        for(size_t i = 0; i < edges.size(); ++i)
        {
            if(edges[i]->getEnd()->getColor() == GC_BLACK)
            {
                // Mark the edge and its twin for removal
                if(edges[i]->getColor() != GC_BLACK || edges[i]->getTwin()->getColor() != GC_BLACK)
                {
                    edges[i]->setColor(GC_BLACK);
                    edges[i]->getTwin()->setColor(GC_BLACK);
                    marked_edges += 2;
                    trans_count++;
                }
            }
            edges[i]->getEnd()->setColor(GC_WHITE);
        }
    }

    if(trans_count > 0)
        ++marked_verts;

    return false;
}

// Remove all the marked edges
void SGTransitiveReductionVisitor::postvisit(StringGraph* pGraph)
{
    printf("TR marked %d verts and %d edges\n", marked_verts, marked_edges);
    pGraph->sweepEdges(GC_BLACK);
    pGraph->setTransitiveFlag(false);
    assert(pGraph->checkColors(GC_WHITE));
}

//
// SGIdenticalRemoveVisitor - Removes identical vertices
// from the graph. This algorithm is less complex 
// than SGContainRemoveVisitor because we do not have to
// remodel the graph in this case because no irreducible
// edges need to be moved.
//
void SGIdenticalRemoveVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    count = 0;
}

// 
bool SGIdenticalRemoveVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    if(!pVertex->isContained())
        return false;

    // Check if this vertex is identical to any other vertex
    EdgePtrVec neighborEdges = pVertex->getEdges();
    for(size_t i = 0; i < neighborEdges.size(); ++i)
    {
        Edge* pEdge = neighborEdges[i];
        Vertex* pOther = pEdge->getEnd();
        if(pVertex->getSeqLen() != pOther->getSeqLen())
            continue;
        
        Overlap ovr = pEdge->getOverlap();
        if(!ovr.isContainment() || ovr.getContainedIdx() != 0)
            continue;

        if(pVertex->getSeq() == pOther->getSeq())
        {
            pVertex->setColor(GC_BLACK);
            ++count;
            break;
        }
    }
            
    return false;
}

void SGIdenticalRemoveVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_BLACK);
}

//
// SGContainRemoveVisitor - Removes contained
// vertices from the graph
//
void SGContainRemoveVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);

    // Clear the containment flag, if any containments are added
    // during this algorithm the flag will be reset and another
    // round must be re-run
    pGraph->setContainmentFlag(false);    
}

//
bool SGContainRemoveVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    if(!pVertex->isContained())
        return false;
    // Add any new irreducible edges that exist when pToRemove is deleted
    // from the graph
    EdgePtrVec neighborEdges = pVertex->getEdges();
    
    // If the graph has been transitively reduced, we have to check all
    // the neighbors to see if any new edges need to be added. If the graph is a
    // complete overlap graph we can just remove the edges to the deletion vertex
    if(!pGraph->hasTransitive() && !pGraph->isExactMode())
    {
        // This must be done in order of edge length or some transitive edges
        // may be created
        EdgeLenComp comp;
        std::sort(neighborEdges.begin(), neighborEdges.end(), comp);

        for(size_t j = 0; j < neighborEdges.size(); ++j)
        {
            Vertex* pRemodelVert = neighborEdges[j]->getEnd();
            Edge* pRemodelEdge = neighborEdges[j]->getTwin();
            SGAlgorithms::remodelVertexForExcision2(pGraph, 
                                                    pRemodelVert, 
                                                    pRemodelEdge);
        }
    }
            
    // Delete the edges from the graph
    for(size_t j = 0; j < neighborEdges.size(); ++j)
    {
        Vertex* pRemodelVert = neighborEdges[j]->getEnd();
        Edge* pRemodelEdge = neighborEdges[j]->getTwin();
        pRemodelVert->deleteEdge(pRemodelEdge);
        pVertex->deleteEdge(neighborEdges[j]);
    }
    pVertex->setColor(GC_BLACK);
    return false;
}

void SGContainRemoveVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_BLACK);
}

//
// Validate the structure of the graph by detecting missing
// or erroneous edges
//

//
bool SGValidateStructureVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    SGAlgorithms::EdgeDescOverlapMap irreducibleMap;
    SGAlgorithms::EdgeDescOverlapMap transitiveMap;
    
    // Construct the set of overlaps reachable within the current parameters
    CompleteOverlapSet vertexOverlapSet(pVertex, pGraph->getErrorRate(), pGraph->getMinOverlap());
    vertexOverlapSet.computeIrreducible(NULL, NULL);

    SGAlgorithms::EdgeDescOverlapMap missingMap;
    SGAlgorithms::EdgeDescOverlapMap extraMap;
    
    vertexOverlapSet.getDiffMap(missingMap, extraMap);  
    
    if(!missingMap.empty())
    {
        std::cout << "Missing irreducible for " << pVertex->getID() << ":\n";
        SGAlgorithms::printOverlapMap(missingMap);
    }

    if(!extraMap.empty())
    {
        std::cout << "Extra irreducible for " << pVertex->getID() << ":\n";
        SGAlgorithms::printOverlapMap(extraMap);
    }
    return false;
}

//
// SGRemodelVisitor - Remodel the graph to infer missing edges or remove erroneous edges
//
void SGRemodelVisitor::previsit(StringGraph* pGraph)
{
    m_remodelER = 0.02;
    pGraph->setColors(GC_WHITE);
}

bool SGRemodelVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    bool graph_changed = false;

    // Construct the set of overlaps reachable within the current parameters
    CompleteOverlapSet vertexOverlapSet(pVertex, m_remodelER, pGraph->getMinOverlap());
    SGAlgorithms::EdgeDescOverlapMap containMap;
    vertexOverlapSet.computeIrreducible(NULL, &containMap);
    SGAlgorithms::EdgeDescOverlapMap irreducibleMap = vertexOverlapSet.getOverlapMap();

    // Construct the set of edges that should be added
    EdgePtrVec edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        SGAlgorithms::EdgeDescOverlapMap::iterator iter = irreducibleMap.find(edges[i]->getDesc());
        if(iter != irreducibleMap.end())
        {
            // Edge exists already
            irreducibleMap.erase(iter);
        }
        else
        {
            edges[i]->setColor(GC_BLACK);
            edges[i]->getTwin()->setColor(GC_BLACK);
            //std::cout << "Marking edge for deletion: " << edges[i]->getOverlap() << "\n";
        }
    }

    // Add remaining edges in the irreducible map
    SGAlgorithms::EdgeDescOverlapMap::iterator iter;
    for(iter = irreducibleMap.begin(); iter != irreducibleMap.end(); ++iter)
    {
        Overlap& ovr = iter->second;
        //std::cout << "Adding overlap: " << ovr << "\n";
        SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, false);
        graph_changed = true;
    }

    // Update the containment flags in the graph to ensure that we can subsequently remove containment verts
    SGAlgorithms::updateContainFlags(pGraph, pVertex, containMap);

    return graph_changed;
}

//
void SGRemodelVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepEdges(GC_BLACK);
    pGraph->setErrorRate(m_remodelER);
}

//
// SGErrorCorrectVisitor - Run error correction on the reads
//
bool SGErrorCorrectVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    static size_t numCorrected = 0;

    if(numCorrected > 0 && numCorrected % 50000 == 0)
        std::cerr << "Corrected " << numCorrected << " reads\n";

    std::string corrected = ErrorCorrect::correctVertex(pGraph, pVertex, 5, 0.01);
    pVertex->setSeq(corrected);
    ++numCorrected;
    return false;
}

//
// SGEdgeStatsVisitor - Compute and display summary statistics of 
// the overlaps in the graph, including edges that were potentially missed
//
void SGEdgeStatsVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    maxDiff = 0;
    minOverlap = pGraph->getMinOverlap();
    maxOverlap = 0;

}

bool SGEdgeStatsVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    const int MIN_OVERLAP = pGraph->getMinOverlap();
    const double MAX_ERROR = pGraph->getErrorRate();

    static int visited = 0;
    ++visited;
    if(visited % 50000 == 0)
        std::cout << "visited: " << visited << "\n";

    // Add stats for the found overlaps
    EdgePtrVec edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Overlap ovr = edges[i]->getOverlap();
        int numDiff = ovr.match.countDifferences(pVertex->getStr(), edges[i]->getEnd()->getStr());
        int overlapLen = ovr.match.getMinOverlapLength();
        addOverlapToCount(overlapLen, numDiff, foundCounts);
    }
        
    // Explore the neighborhood around this graph for potentially missing overlaps
    CandidateVector candidates = getMissingCandidates(pGraph, pVertex, MIN_OVERLAP);
    MultiOverlap addedMO(pVertex->getID(), pVertex->getStr());
    for(size_t i = 0; i < candidates.size(); ++i)
    {
        Candidate& c = candidates[i];
        int numDiff = c.ovr.match.countDifferences(pVertex->getStr(), c.pEndpoint->getStr());
        double error_rate = double(numDiff) / double(c.ovr.match.getMinOverlapLength());

        if(error_rate < MAX_ERROR)
        {
            int overlapLen = c.ovr.match.getMinOverlapLength();
            addOverlapToCount(overlapLen, numDiff, missingCounts);
        }
    }
    
    return false;
}

//
void SGEdgeStatsVisitor::postvisit(StringGraph* /*pGraph*/)
{    
    printf("FoundOverlaps\n");
    printCounts(foundCounts);

    printf("\nPotentially Missing Overlaps\n\n");
    printCounts(missingCounts);
}

//
void SGEdgeStatsVisitor::printCounts(CountMatrix& matrix)
{
    // Header row
    printf("OL\t");
    for(int j = 0; j <= maxDiff; ++j)
    {
        printf("%d\t", j);
    }

    printf("sum\n");
    IntIntMap columnTotal;
    for(int i = minOverlap; i <= maxOverlap; ++i)
    {
        printf("%d\t", i);
        int sum = 0;
        for(int j = 0; j <= maxDiff; ++j)
        {
            int v = matrix[i][j];
            printf("%d\t", v);
            sum += v;
            columnTotal[j] += v;
        }
        printf("%d\n", sum);
    }

    printf("total\t");
    int total = 0;
    for(int j = 0; j <= maxDiff; ++j)
    {
        int v = columnTotal[j];
        printf("%d\t", v);
        total += v;
    }
    printf("%d\n", total);
}

//
void SGEdgeStatsVisitor::addOverlapToCount(int ol, int nd, CountMatrix& matrix)
{
    matrix[ol][nd]++;

    if(nd > maxDiff)
        maxDiff = nd;

    if(ol > maxOverlap)
        maxOverlap = ol;
}

// Explore the neighborhood around a vertex looking for missing overlaps
SGEdgeStatsVisitor::CandidateVector SGEdgeStatsVisitor::getMissingCandidates(StringGraph* /*pGraph*/, 
                                                                             Vertex* pVertex, 
                                                                             int minOverlap) const
{
    CandidateVector out;

    // Mark the vertices that are reached from this vertex as black to indicate
    // they already are overlapping
    EdgePtrVec edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        edges[i]->getEnd()->setColor(GC_BLACK);
    }
    pVertex->setColor(GC_BLACK);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pXY = edges[i];
        EdgePtrVec neighborEdges = pXY->getEnd()->getEdges();
        for(size_t j = 0; j < neighborEdges.size(); ++j)
        {
            Edge* pYZ = neighborEdges[j];
            if(pYZ->getEnd()->getColor() != GC_BLACK)
            {
                // Infer the overlap object from the edges
                Overlap ovrXY = pXY->getOverlap();
                Overlap ovrYZ = pYZ->getOverlap();

                if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
                {
                    Overlap ovr_xz = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
                    if(ovr_xz.match.getMinOverlapLength() >= minOverlap)
                    {
                        out.push_back(Candidate(pYZ->getEnd(), ovr_xz));
                        pYZ->getEnd()->setColor(GC_BLACK);
                    }
                }
            }
        }
    }

    // Reset colors
    for(size_t i = 0; i < edges.size(); ++i)
        edges[i]->getEnd()->setColor(GC_WHITE);
    pVertex->setColor(GC_WHITE);
    for(size_t i = 0; i < out.size(); ++i)
        out[i].pEndpoint->setColor(GC_WHITE);
    return out;
}

//
// SGTrimVisitor - Remove "dead-end" vertices from the graph
//
void SGTrimVisitor::previsit(StringGraph* pGraph)
{
    num_island = 0;
    num_terminal = 0;
    pGraph->setColors(GC_WHITE);
}

// Mark any nodes that either dont have edges or edges in only one direction for removal
bool SGTrimVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    if(pVertex->countEdges() == 0)
    {
        // Is an island, remove if the sequence length is less than the threshold
        if(pVertex->getSeqLen() < m_minLength)
        {
            pVertex->setColor(GC_BLACK);
            ++num_island;
        }
    }
    else
    {
        // Check if this node is a dead-end
        for(size_t idx = 0; idx < ED_COUNT; idx++)
        {
            EdgeDir dir = EDGE_DIRECTIONS[idx];
            if(pVertex->countEdges(dir) == 0 && pVertex->getSeqLen() < m_minLength)
            {
                pVertex->setColor(GC_BLACK);
                ++num_terminal;
            }
        }
    }

    /*
    if(isTrimmed)
    {
        std::stringstream ss;
        ss << pVertex->getID() << "-trimmed";
        writeFastaRecord(&m_tmpFile, ss.str(), pVertex->getSeq().toString());
    }
    */
    return false;
}

// Remove all the marked edges
void SGTrimVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_BLACK);
    printf("StringGraphTrim: Removed %d island and %d dead-end short vertices\n", num_island, num_terminal);
}

//
// SGDuplicateVisitor - Detect and remove duplicate edges
//
void SGDuplicateVisitor::previsit(StringGraph* pGraph)
{
    assert(pGraph->checkColors(GC_WHITE));
    (void)pGraph;
    m_hasDuplicate = false;
}

bool SGDuplicateVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    m_hasDuplicate = pVertex->markDuplicateEdges(GC_RED) || m_hasDuplicate;
    return false;
}


void SGDuplicateVisitor::postvisit(StringGraph* pGraph)
{
    assert(pGraph->checkColors(GC_WHITE));
    if(m_hasDuplicate)
    {
        int numRemoved = pGraph->sweepEdges(GC_RED);
        if(!m_bSilent)
            std::cerr << "Warning: removed " << numRemoved << " duplicate edges\n";
    }
}

//
// SGIslandVisitor - Remove island (unconnected) vertices
//
void SGIslandVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
}

// Mark any nodes that dont have edges
bool SGIslandVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    if(pVertex->countEdges() == 0)
    {
        pVertex->setColor(GC_BLACK);
        return true;
    }
    return false;
}

// Remove all the marked vertices
void SGIslandVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_BLACK);
}

//
// Small repeat resolver - Remove edges induced from small (sub-read length)
// repeats
// 
void SGSmallRepeatResolveVisitor::previsit(StringGraph*)
{

}

//
bool SGSmallRepeatResolveVisitor::visit(StringGraph* /*pGraph*/, Vertex* pX)
{
    bool changed = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec x_edges = pX->getEdges(dir); // These edges are already sorted
        if(x_edges.size() < 2)
            continue;

        // Try to eliminate the shortest edge from this vertex (let this be X->Y)
        // If Y has a longer edge than Y->X in the same direction, we remove X->Y

        // Edges are sorted by length so the last edge is the shortest
        Edge* pXY = x_edges.back();
        size_t xy_len =  pXY->getOverlap().getOverlapLength(0);
        size_t x_longest_len = x_edges.front()->getOverlap().getOverlapLength(0);
        if(xy_len == x_longest_len)
            continue;

        Edge* pYX = pXY->getTwin();
        Vertex* pY = pXY->getEnd();

        EdgePtrVec y_edges = pY->getEdges(pYX->getDir());
        size_t yx_len = pYX->getOverlap().getOverlapLength(0);

        size_t y_longest_len = 0;
        for(size_t i = 0; i < y_edges.size(); ++i)
        {
            Edge* pYZ = y_edges[i];
            if(pYZ == pYX)
                continue; // skip Y->X

            size_t yz_len = pYZ->getOverlap().getOverlapLength(0);
            if(yz_len > y_longest_len)
                y_longest_len = yz_len;
        }


        if(y_longest_len > yx_len)
        {
            // Delete the edge if the difference between the shortest and longest is greater than minDiff
            int x_diff = x_longest_len - xy_len;
            int y_diff = y_longest_len - yx_len;

            if(x_diff > m_minDiff && y_diff > m_minDiff)
            {
                /*
                printf("Edge %s -> %s is likely a repeat\n", pX->getID().c_str(), pY->getID().c_str());
                printf("Actual overlap lengths: %zu and %zu\n", xy_len, yx_len);
                printf("Spanned by longer edges of size: %zu and %zu\n", x_longest_len, y_longest_len);
                printf("Differences: %d and %d\n", x_diff, y_diff);
                */
                pX->deleteEdge(pXY);
                pY->deleteEdge(pYX);
                changed = true;
            }
        }
    }

    return changed;
}

//
void SGSmallRepeatResolveVisitor::postvisit(StringGraph*)
{

}

//
// SGBubbleVisitor - Find and collapse variant
// "bubbles" in the graph
//
void SGBubbleVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    num_bubbles = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGBubbleVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    bool bubble_found = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pVertex->getEdges(dir);
        if(edges.size() > 1)
        {
            Vertex* pStart = pVertex;
            Vertex* pEnd = NULL;
            // Check the vertices
            for(size_t i = 0; i < edges.size(); ++i)
            {
                Edge* pVWEdge = edges[i];
                Vertex* pWVert = pVWEdge->getEnd();

                // Get the edges from w in the same direction
                EdgeDir transDir = !pVWEdge->getTwinDir();
                EdgePtrVec wEdges = pWVert->getEdges(transDir);

                if(pWVert->getColor() == GC_RED)
                    return false;

                // If the bubble has collapsed, there should only be one edge
                if(wEdges.size() == 1)
                {
                    Vertex* pBubbleEnd = wEdges.front()->getEnd();
                    if(pBubbleEnd->getColor() == GC_RED)
                        return false;
                }
            }

            // Mark the vertices
            for(size_t i = 0; i < edges.size(); ++i)
            {
                Edge* pVWEdge = edges[i];
                Vertex* pWVert = pVWEdge->getEnd();

                // Get the edges from w in the same direction
                EdgeDir transDir = !pVWEdge->getTwinDir();
                EdgePtrVec wEdges = pWVert->getEdges(transDir);

                // If the bubble has collapsed, there should only be one edge
                if(wEdges.size() == 1)
                {
                    Vertex* pBubbleEnd = wEdges.front()->getEnd();
                    if(pBubbleEnd->getColor() == GC_BLACK)
                    {
                        // The endpoint has been visited, set this vertex as needing removal
                        // and set the endpoint as unvisited
                        pWVert->setColor(GC_RED);
                        bubble_found = true;
                        pEnd = pBubbleEnd;
                    }
                    else
                    {
                        pBubbleEnd->setColor(GC_BLACK);
                        pWVert->setColor(GC_BLUE);
                    }
                }
            }
            
            // Unmark vertices
            for(size_t i = 0; i < edges.size(); ++i)
            {
                Edge* pVWEdge = edges[i];
                Vertex* pWVert = pVWEdge->getEnd();

                // Get the edges from w in the same direction
                EdgeDir transDir = !pVWEdge->getTwinDir();
                EdgePtrVec wEdges = pWVert->getEdges(transDir);

                // If the bubble has collapsed, there should only be one edge
                if(wEdges.size() == 1)
                {
                    Vertex* pBubbleEnd = wEdges.front()->getEnd();
                    pBubbleEnd->setColor(GC_WHITE);
                }
                if(pWVert->getColor() == GC_BLUE)
                    pWVert->setColor(GC_WHITE);
            }

            (void)pStart;
            (void)pEnd;
            if(bubble_found)
            {
                /*
                SGWalkVector walkVector;
                SGSearch::findWalks(pStart, pEnd, dir, 1000, 20, walkVector);
                
                if(walkVector.size() == 2)
                {
                    SGWalk& walk1 = walkVector[0];
                    SGWalk& walk2 = walkVector[1];

                    int len1 = walk1.getStartToEndDistance();
                    int len2 = walk2.getStartToEndDistance();
                    int diff = len1 - len2;
                    std::string type = "SNP";
                    if(diff != 0)
                    {
                        type = "INDEL";
                    }
                    std::cout << "Bubble " << pStart->getID() << " to " << pEnd->getID() << " is a "
                              << type << "(d: " << diff << ")\n";
                }
                */
                ++num_bubbles;
            }
        }
    }
    return bubble_found;
}

// Remove all the marked edges
void SGBubbleVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_RED);
    printf("bubbles: %d\n", num_bubbles);
    assert(pGraph->checkColors(GC_WHITE));
}

//
// SGBubbleEdgeVisitor 
//
void SGBubbleEdgeVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    num_bubbles = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGBubbleEdgeVisitor::visit(StringGraph* /*pGraph*/, Vertex* pX)
{    
    bool bubble_found = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pX->getEdges(dir);
        if(edges.size() == 2) // di-bubbles only for now
        {
            // Determine which edge has a shorter overlap to pX
            // Call the longer overlap pY, the shorter pZ
            Edge* pXY;
            Edge* pXZ;

            if(edges[0]->getOverlap().getOverlapLength(0) > edges[1]->getOverlap().getOverlapLength(0))
            {
                pXY = edges[0];
                pXZ = edges[1];
            }
            else if(edges[1]->getOverlap().getOverlapLength(0) > edges[0]->getOverlap().getOverlapLength(0))

            {
                pXY = edges[1];
                pXZ = edges[0];
            }
            else
            {
                break; // equal length overlaps, cannot be a bubble or else the vertices would be contained
            }
            
            // Mark the neighbors of pZ as the "target" vertices
            // if they can be reached by pY we mark pY as being unreliable and remove it
            typedef std::list<Vertex*> VertexPtrList;
            VertexPtrList targetList;

            EdgeDir targetDir = pXZ->getTransitiveDir();
            EdgePtrVec targetEdges = pXZ->getEnd()->getEdges(targetDir);
            for(size_t i = 0; i < targetEdges.size(); ++i)
                targetList.push_back(targetEdges[i]->getEnd());

            // Start exploring from pY
            ExploreQueue queue;
            Overlap ovrXY = pXY->getOverlap();
            EdgeDesc edXY = pXY->getDesc();
            queue.push(ExploreElement(edXY, ovrXY));

            int numSteps = 100;
            WARN_ONCE("USING FIXED NUMBER OF STEPS IN BUBBLE EDGE");
            while(!queue.empty() && numSteps-- > 0)
            {
                ExploreElement ee = queue.front();
                EdgeDesc& edXY = ee.ed;
                Vertex* pY = edXY.pVertex;
                Overlap& ovrXY = ee.ovr;

                queue.pop();

                // Check if Y is on the target list
                VertexPtrList::iterator iter = targetList.begin();
                while(iter != targetList.end())
                {
                    if(*iter == edXY.pVertex)
                        targetList.erase(iter++);
                    else
                        ++iter;
                }
                
                if(targetList.empty())
                    break;

                // Enqueue the neighbors of pY
                EdgeDir dirY = edXY.getTransitiveDir();
                EdgePtrVec edges = pY->getEdges(dirY);
                for(size_t i = 0; i < edges.size(); ++i)
                {
                    Edge* pEdge = edges[i];
                    Vertex* pZ = pEdge->getEnd();

                    // Compute the edgeDesc and overlap on pX for this edge
                    Overlap ovrYZ = pEdge->getOverlap();

                    if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
                    {
                        Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
                        EdgeDesc edXZ = SGAlgorithms::overlapToEdgeDesc(pZ, ovrXZ);
                        queue.push(ExploreElement(edXZ, ovrXZ));
                    }
                }
            }

            if(targetList.empty())
            {
                // bubble found
                pXZ->getEnd()->deleteEdges();
                pXZ->getEnd()->setColor(GC_RED);
                bubble_found = true;
                ++num_bubbles;
            }
        }
    }
    return bubble_found;
}

// Remove all the marked vertices
void SGBubbleEdgeVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_RED);
    printf("bubbles: %d\n", num_bubbles);
    assert(pGraph->checkColors(GC_WHITE));
}

//
// SGSmoothingVisitor - Find branches in the graph
// which arise from variation and remove them
//
void SGSmoothingVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    m_simpleBubblesRemoved = 0;
    m_complexBubblesRemoved = 0;
}

//
bool SGSmoothingVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    (void)pGraph;
    if(pVertex->getColor() == GC_RED)
        return false;

    bool found = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pVertex->getEdges(dir);
        if(edges.size() <= 1)
            continue;

        for(size_t i = 0; i < edges.size(); ++i)
        {
            if(edges[i]->getEnd()->getColor() == GC_RED)
                return false;
        }

        //std::cout << "Smoothing " << pVertex->getID() << "\n";

        const int MAX_WALKS = 10;
        const int MAX_DISTANCE = 5000;
        bool bIsDegenerate = false;
        bool bFailGapCheck = false;
        bool bFailDivergenceCheck = false;
        bool bFailIndelSizeCheck = false;

        SGWalkVector variantWalks;
        SGSearch::findVariantWalks(pVertex, dir, MAX_DISTANCE, MAX_WALKS, variantWalks);

        if(variantWalks.size() > 0)
        {
            found = true;
            size_t selectedIdx = -1;
            size_t selectedCoverage = 0;

            // Calculate the minimum amount overlapped on the start/end vertex.
            // This is used to properly extract the sequences from walks that represent the variation.
            int minOverlapX = std::numeric_limits<int>::max();
            int minOverlapY = std::numeric_limits<int>::max();

            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                if(variantWalks[i].getNumEdges() <= 1)
                    bIsDegenerate = true;

                // Calculate the walk coverage using the internal vertices of the walk. 
                // The walk with the highest coverage will be retained
                size_t walkCoverage = 0;
                for(size_t j = 1; j < variantWalks[i].getNumVertices() - 1; ++j)
                    walkCoverage += variantWalks[i].getVertex(j)->getCoverage();

                if(walkCoverage > selectedCoverage)
                {
                    selectedIdx = i;
                    selectedCoverage = walkCoverage;
                }
                
                Edge* pFirstEdge = variantWalks[i].getFirstEdge();
                Edge* pLastEdge = variantWalks[i].getLastEdge();

                if((int)pFirstEdge->getMatchLength() < minOverlapX)
                    minOverlapX = pFirstEdge->getMatchLength();

                if((int)pLastEdge->getTwin()->getMatchLength() < minOverlapY)
                    minOverlapY = pLastEdge->getTwin()->getMatchLength();
            }

            // Calculate the strings for each walk that represent the region of variation
            StringVector walkStrings;
            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                Vertex* pStartVertex = variantWalks[i].getStartVertex();
                Vertex* pLastVertex = variantWalks[i].getLastVertex();
                assert(pStartVertex != NULL && pLastVertex != NULL);
                
                std::string full = variantWalks[i].getString(SGWT_START_TO_END);
                int posStart = 0;
                int posEnd = 0;

                if(dir == ED_ANTISENSE)
                {
                    // pLast   -----------
                    // pStart          ------------
                    // full    --------------------
                    // out             ----
                    posStart = pLastVertex->getSeqLen() - minOverlapY;
                    posEnd = full.size() - (pStartVertex->getSeqLen() - minOverlapX);
                }
                else
                {
                    // pStart         --------------
                    // pLast   -----------
                    // full    ---------------------
                    // out            ----
                    posStart = pStartVertex->getSeqLen() - minOverlapX; // match start position
                    posEnd = full.size() - (pLastVertex->getSeqLen() - minOverlapY); // match end position
                }
                
                std::string out;
                if(posEnd > posStart)
                    out = full.substr(posStart, posEnd - posStart);
                walkStrings.push_back(out);
            }

            assert(selectedIdx != (size_t)-1);
            SGWalk& selectedWalk = variantWalks[selectedIdx];
            assert(selectedWalk.isIndexed());

            // Check the divergence of the other walks to this walk
            StringVector cigarStrings;
            std::vector<int> maxIndel;
            std::vector<double> gapPercent; // percentage of matching that is gaps
            std::vector<double> totalPercent; // percent of total alignment that is mismatch or gap

            cigarStrings.resize(variantWalks.size());
            gapPercent.resize(variantWalks.size());
            totalPercent.resize(variantWalks.size());
            maxIndel.resize(variantWalks.size());

            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                if(i == selectedIdx)
                    continue;

                // We want to compute the total gap length, total mismatches and percent
                // divergence between the two paths.
                int matchLen = 0;
                int totalDiff = 0;
                int gapLength = 0;
                int maxGapLength = 0;
                // We have to handle the degenerate case where one internal string has zero length
                // this can happen when there is an isolated insertion/deletion and the walks are like:
                // x -> y -> z
                // x -> z
                if(walkStrings[selectedIdx].empty() || walkStrings[i].empty())
                {
                    matchLen = std::max(walkStrings[selectedIdx].size(), walkStrings[i].size());
                    totalDiff = matchLen;
                    gapLength = matchLen;
                }
                else
                {
                    AlnAln *aln_global;
                    aln_global = aln_stdaln(walkStrings[selectedIdx].c_str(), walkStrings[i].c_str(), &aln_param_blast, 1, 1);

                    // Calculate the alignment parameters
                    while(aln_global->outm[matchLen] != '\0')
                    {
                        if(aln_global->outm[matchLen] == ' ')
                            totalDiff += 1;
                        matchLen += 1;
                    }

                    std::stringstream cigarSS;
                    for (int j = 0; j != aln_global->n_cigar; ++j)
                    {
                        char cigarOp = "MID"[aln_global->cigar32[j]&0xf];
                        int cigarLen = aln_global->cigar32[j]>>4;
                        if(cigarOp == 'I' || cigarOp == 'D')
                        {
                            gapLength += cigarLen;
                            if(gapLength > maxGapLength)
                                maxGapLength = gapLength;
                        }

                        cigarSS << cigarLen;
                        cigarSS << cigarOp;
                    }
                    cigarStrings[i] = cigarSS.str();

                    /*
                    printf("1: %s\n", aln_global->out1);
                    printf("M: %s\n", aln_global->outm);
                    printf("2: %s\n", aln_global->out2);
                    printf("CIGAR: %s\n", cigarStrings[i].c_str());
                    */

                    aln_free_AlnAln(aln_global);
                }

                double percentDiff = (double)totalDiff / matchLen;
                double percentGap = (double)gapLength / matchLen;

                if(percentDiff > m_maxTotalDivergence)
                    bFailDivergenceCheck = true;
                
                if(percentGap > m_maxGapDivergence)
                    bFailGapCheck = true;

                if(maxGapLength > m_maxIndelLength)
                    bFailIndelSizeCheck = true;

                gapPercent[i] = percentGap;
                totalPercent[i] = percentDiff;
                maxIndel[i] = maxGapLength;

                //printf("ml: %d tmm: %d pd: %lf pg: %lf\n", matchLen, totalDiff, percentDiff, percentGap);
            }

            if(bIsDegenerate || bFailGapCheck || bFailDivergenceCheck || bFailIndelSizeCheck)
                continue;

            // Write the selected path to the variants file as variant 0
            int variantIdx = 0;
            std::string selectedSequence = selectedWalk.getString(SGWT_START_TO_END);
            std::stringstream ss;
            ss << "variant-" << m_numRemovedTotal << "/" << variantIdx++;
            writeFastaRecord(&m_outFile, ss.str(), selectedSequence);


            // The vertex set for each walk is not necessarily disjoint,
            // the selected walk may contain vertices that are part
            // of other paths. We handle this be initially marking all
            // vertices of the 
            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                if(i == selectedIdx)
                    continue;

                SGWalk& currWalk = variantWalks[i];
                for(size_t j = 0; j < currWalk.getNumEdges() - 1; ++j)
                {
                    Edge* currEdge = currWalk.getEdge(j);
                    
                    // If the vertex is also on the selected path, do not mark it
                    Vertex* currVertex = currEdge->getEnd();
                    if(!selectedWalk.containsVertex(currVertex->getID()))
                    {
                        currEdge->getEnd()->setColor(GC_RED);
                    }
                }

                // Write the variant to a file
                std::string variantSequence = currWalk.getString(SGWT_START_TO_END);
                std::stringstream variantID;
                std::stringstream ss;
                ss << "variant-" << m_numRemovedTotal << "/" << variantIdx++;
                ss << " IGD:" << (double)gapPercent[i] << " ITD:" << totalPercent[i] << " MID: " << maxIndel[i] << " InternalCigar:" << cigarStrings[i];
                writeFastaRecord(&m_outFile, ss.str(), variantSequence);
            }

            if(variantWalks.size() == 2)
                m_simpleBubblesRemoved += 1;
            else
                m_complexBubblesRemoved += 1;
            ++m_numRemovedTotal;
        }
    }
    return found;
}

// Remove all the marked edges
void SGSmoothingVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_RED);
    assert(pGraph->checkColors(GC_WHITE));

    printf("VariationSmoother: Removed %d simple and %d complex bubbles\n", m_simpleBubblesRemoved, m_complexBubblesRemoved);
}


//
// Coverage analysis
//
void SGCoverageVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    m_numRemoved = 0;
}

//
bool SGCoverageVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pVertex->getEdges(dir);
        if(edges.size() <= 1)
            continue;

        int bestCoverage = 0;
        int worstCoverage = std::numeric_limits<int>::max();
        int popIndex = -1;

        std::cout << "Bubble coverage:\n";
        for(size_t i = 0; i < edges.size(); ++i)
        {
            int coverage = SGSearch::countSpanningCoverage(edges[i], 20);

            std::cout << "\t" << i << " cov: " << coverage << "\n";
            
            if(coverage > bestCoverage)
                bestCoverage = coverage;

            if(coverage < worstCoverage)
            {
                worstCoverage = coverage;
                popIndex = i;
            }            
        }

        if(worstCoverage > 0 && worstCoverage < m_cutoff && bestCoverage > m_cutoff)
        {
            assert(popIndex != -1);
            if(edges[popIndex]->getEnd()->getColor() != GC_RED)
            {
                edges[popIndex]->getEnd()->setColor(GC_RED);
                ++m_numRemoved;
            }
        }
    }
    return false;
}

//
void SGCoverageVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_RED);
    printf("Removed %d low-coverage nodes\n", m_numRemoved);
}

//
// SGGraphStatsVisitor - Collect summary stasitics
// about the graph
//
void SGGraphStatsVisitor::previsit(StringGraph* /*pGraph*/)
{
    num_terminal = 0;
    num_island = 0;
    num_monobranch = 0;
    num_dibranch = 0;
    num_transitive = 0;
    num_edges = 0;
    num_vertex = 0;
    sum_edgeLen = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGGraphStatsVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    int s_count = pVertex->countEdges(ED_SENSE);
    int as_count = pVertex->countEdges(ED_ANTISENSE);
    if(s_count == 0 && as_count == 0)
    {
        ++num_island;
    }
    else if(s_count == 0 || as_count == 0)
    {
        ++num_terminal;
    }

    if(s_count > 1 && as_count > 1)
        ++num_dibranch;
    else if(s_count > 1 || as_count > 1)
        ++num_monobranch;

    if(s_count == 1 || as_count == 1)
        ++num_transitive;

    num_edges += (s_count + as_count);
    ++num_vertex;

    EdgePtrVec edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
        sum_edgeLen += edges[i]->getSeqLen();

    return false;
}

//
void SGGraphStatsVisitor::postvisit(StringGraph* /*pGraph*/)
{
    printf("island: %d terminal: %d monobranch: %d dibranch: %d transitive: %d\n", num_island, num_terminal,
                                                                                   num_monobranch, num_dibranch, num_transitive);
    printf("Total Vertices: %d Total Edges: %d Sum edge length: %zu\n", num_vertex, num_edges, sum_edgeLen);
}

//
bool SGBreakWriteVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    int s_count = pVertex->countEdges(ED_SENSE);
    int as_count = pVertex->countEdges(ED_ANTISENSE);

    if(s_count == 0 && as_count == 0)
    {
        writeBreak("ISLAND", pVertex);
    }
    else if(s_count == 0)
    {
        writeBreak("STIP", pVertex);
    }
    else if(as_count == 0)
    {
        writeBreak("ASTIP", pVertex);
    }

    if(s_count > 1)
    {
        std::stringstream text;
        text << "SBRANCHED," << calculateOverlapLengthDifference(pVertex, ED_SENSE);
        writeBreak(text.str(), pVertex);
    }
    
    if(as_count > 1)
    {
        std::stringstream text;
        text << "ASBRANCHED," << calculateOverlapLengthDifference(pVertex, ED_ANTISENSE);
        writeBreak(text.str(), pVertex);
    }
    return false;
}

int SGBreakWriteVisitor::calculateOverlapLengthDifference(const Vertex* pVertex, EdgeDir dir)
{
    EdgePtrVec edges = pVertex->getEdges(dir);
    if(edges.size() < 2)
        return 0;
    int shortestLen = edges[edges.size() - 1]->getOverlap().getOverlapLength(0);
    int secondLen = edges[edges.size() - 2]->getOverlap().getOverlapLength(0);
    return secondLen - shortestLen;
}

void SGBreakWriteVisitor::writeBreak(const std::string& type, Vertex* pVertex)
{
    *m_pWriter << "BREAK\t" << type << "\t" << pVertex->getID() << "\t" << pVertex->getSeq() << "\n";
}
