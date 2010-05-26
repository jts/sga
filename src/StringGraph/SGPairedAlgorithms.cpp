//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGPairedAlgorithms - Collection of algorithms
// using paired end data in string graphs
//
#include "SGPairedAlgorithms.h"
#include <queue>

//
// Generic paired end algorithms
//

// Perform a depth-first search for paths between pX and pY
// Terminate the search when the total path length (the length of the
// string implied by the path) is greater than maxDistance
struct SearchNode
{
    Path path;
    int distance;
};

typedef std::queue<SearchNode> SearchQueue;
void SGPairedAlgorithms::searchPaths(const Vertex* pX, const Vertex* pY, 
                                     int maxDistance, PathVector& outPaths)
{
    //std::cout << "Resolving path from " << pX->getID() << " to " << pY->getID() << "\n";
    // Get the initial search direction
    EdgeDir dir = getDirectionToPair(pX->getID());

    // Create the initial path nodes
    SearchQueue queue;
    EdgePtrVec edges = pX->getEdges(dir);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        assert(!pEdge->getOverlap().isContainment());

        SearchNode node;
        node.path.push_back(pEdge);
        node.distance = pX->getSeqLen() + pEdge->getSeqLen();
        queue.push(node);
    }

    while(!queue.empty())
    {
        SearchNode& currNode = queue.front();
        
        // Check if we have found pY or exceeded the distance
        Edge* pWZ = currNode.path.back(); 
        Vertex* pZ = pWZ->getEnd();
        if(pZ == pY)
        {
            outPaths.push_back(currNode.path);
            queue.pop();
        }
        else if(currNode.distance > maxDistance)
        {
            queue.pop();
        }
        else
        {
            // Path is still valid, continue
            EdgeDir continueDir = pWZ->getTransitiveDir();
            EdgePtrVec zEdges = pZ->getEdges(continueDir);

            if(!zEdges.empty())
            {
                // Update curr node with the first edge
                Edge* pNextEdge = zEdges[0];
                currNode.distance += pNextEdge->getSeqLen();
                currNode.path.push_back(pNextEdge);

                // 
                for(size_t i = 1; i < zEdges.size(); ++i)
                {
                    SearchNode branch = currNode;
                    Edge* pBranchEdge = zEdges[i];
                    branch.distance += pBranchEdge->getSeqLen();
                    branch.path.push_back(pBranchEdge);
                    queue.push(branch);
                }
            }
            else
            {
                queue.pop();
            }
        }
    }
}

// Return the direction to the pair of this sequence based on the ID
EdgeDir SGPairedAlgorithms::getDirectionToPair(const std::string& /*id*/)
{
    // Assumes illumina PE format for now, both reads
    // point in the sense direction towards search other
    return ED_SENSE;
}

std::string SGPairedAlgorithms::pathToString(const Vertex* pX, const Path& path)
{
    std::string out = pX->getSeq().toString();
    EdgeComp currComp = EC_SAME;

    for(size_t i = 0; i < path.size(); ++i)
    {
        Edge* pYZ = path[i];
        EdgeComp ecYZ = pYZ->getComp();

        // Calculate the next comp, between X and Z
        EdgeComp ecXZ;
        if(ecYZ == EC_SAME)
            ecXZ = currComp;
        else
            ecXZ = !currComp;
        
        // The string return is in Y's frame
        // If the lastComp (which is
        std::string edge_str = pYZ->getLabel();
        assert(edge_str.size() != 0);
        if(currComp == EC_REVERSE)
            edge_str = reverseComplement(edge_str);

        out.append(edge_str);
        currComp = ecXZ;
    }
    return out;
}

//
// SGPairedPathResolveVisitor
//
SGPairedPathResolveVisitor::SGPairedPathResolveVisitor()
{
    m_pWriter = createWriter("fragments.fa");
}

//
SGPairedPathResolveVisitor::~SGPairedPathResolveVisitor()
{
    delete m_pWriter;
    m_pWriter = NULL;
}

//
void SGPairedPathResolveVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
}


bool SGPairedPathResolveVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    if(pVertex->getColor() == GC_BLACK)
        return false; // has been resolved already

    // Get the vertex of the pair
    std::string pairID = getPairID(pVertex->getID());
    Vertex* pPair = pGraph->getVertex(pairID);
    if(pPair != NULL)
    {
        PathVector paths;
        SGPairedAlgorithms::searchPaths(pVertex, pPair, 300, paths);   
        pVertex->setColor(GC_BLACK);
        pPair->setColor(GC_BLACK);

        std::cout << "Found " << paths.size() << " paths from " << pVertex->getID()
                  << " to " << pPair->getID() << "\n";

        
        if(paths.size() == 1)
        {
            std::string fragment = SGPairedAlgorithms::pathToString(pVertex, paths[0]);
            SeqRecord record;
            record.id = pVertex->getID();
            record.seq = fragment;
            record.write(*m_pWriter);
        }
        else
        {
            SeqRecord recordX;
            recordX.id = pVertex->getID();
            recordX.seq = pVertex->getSeq().toString();
            recordX.write(*m_pWriter);

            SeqRecord recordY;
            recordY.id = pVertex->getID();
            recordY.seq = pVertex->getSeq().toString();
            recordY.write(*m_pWriter);
        }
    }

    return false;
}

void SGPairedPathResolveVisitor::postvisit(StringGraph*)
{
}

//
// SGVertexPairingVisitor - links paired vertices in the SG
//
void SGVertexPairingVisitor::previsit(StringGraph* pGraph)
{
    num_paired = 0;
    num_unpaired = 0;
    pGraph->setColors(GC_WHITE);
}

// Visit each vertex in the graph, find its pair and link them
bool SGVertexPairingVisitor::visit(StringGraph* /*pGraph*/, Vertex* /*pVertex*/)
{
#if 0
    // do nothing if the pairing was already set
    if(pVertex->getPairVertex() == NULL)
    {
        std::string id = pVertex->getID();
        std::string pid = getPairID(id);
        Vertex* pPairV = pGraph->getVertex(pid);
        if(pPairV != NULL)
        {
            pVertex->setPairVertex(pPairV);
            pPairV->setPairVertex(pVertex);
            num_paired++;
        }
        else
        {
            //pVertex->setColor(GC_BLACK);
            num_unpaired++;
        }
    }
#endif
    return false;
}

void SGVertexPairingVisitor::postvisit(StringGraph* /*pGraph*/)
{
    printf("Graph has %d paired vertices, %d verts with no pairs\n", 
           num_paired, num_unpaired);
}

//
// SGPETrustVisitor - determines which edges in the 
// string graph are "trusted" - the reads overlapping
// in the edge have pairs that also overlap
//
bool SGPETrustVisitor::visit(StringGraph* /*pGraph*/, Vertex* /*pVertex*/)
{
#if 0
    Vertex* pPairVertex = pVertex->getPairVertex();
    if(pPairVertex == NULL)
        return false;

    // First, mark all pair vertices that overlap the pair of this node
    // The set of marked vertices that overlap pVertex are the trusted vertices
    EdgePtrVec pairEdgeVec = pPairVertex->getEdges();
    for(size_t i = 0; i < pairEdgeVec.size(); ++i)
    {
        // Get the pair of the endpoint of this edge
        Vertex* pBackVertex = pairEdgeVec[i]->getEnd()->getPairVertex();
        if(pBackVertex != NULL)
            pBackVertex->setColor(GC_RED);
    }

    EdgePtrVec vertEdgeVec = pVertex->getEdges();
    
    bool changed = true;
    while(changed)
    {
        changed = false;
        // Propogate trust
        for(size_t i = 0; i < vertEdgeVec.size(); ++i)
        {
            Vertex* pCurr = vertEdgeVec[i]->getEnd();
            if(pCurr->getColor() != GC_RED)
            {
                // If any vertex that pCurr overlaps with is red, mark it red too
                EdgePtrVec currEdgeVec = pCurr->getEdges();
                for(size_t j = 0; j < currEdgeVec.size(); ++j)
                {
                    if(currEdgeVec[j]->getEnd()->getColor() == GC_RED)
                    {
                        pCurr->setColor(GC_RED);
                        changed = true;
                        break;
                    }
                }
            }
        }
    }

    // 
    int trusted = 0;
    int nottrusted = 0;
    int diffstrand = 0;
    for(size_t i = 0; i < vertEdgeVec.size(); ++i)
    {
        if(vertEdgeVec[i]->getEnd()->getColor() == GC_RED)
        {
            trusted++;
            vertEdgeVec[i]->isTrusted = true;
        }
        else
        {
            nottrusted++;
        }
    }

    (void)diffstrand;
    //printf("TOKEN\t%d\t%d\t%d\t%zu\n", trusted, nottrusted, diffstrand, vertEdgeVec.size());

    // Reset all the vertex colors
    for(size_t i = 0; i < pairEdgeVec.size(); ++i)
    {
        // Get the pair of the endpoint of this edge
        Vertex* pBackVertex = pairEdgeVec[i]->getEnd()->getPairVertex();
        if(pBackVertex)
            pBackVertex->setColor(GC_WHITE);
    }

    for(size_t i = 0; i < vertEdgeVec.size(); ++i)
        vertEdgeVec[i]->getEnd()->setColor(GC_WHITE);
#endif
    return false;
}

//
// SGPEConflictRemover - experimental algorithm
// to remove conflicted edges from the graph based
// on paired end data
//
void SGPEConflictRemover::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    num_same = 0;
    num_diff = 0;
}

bool SGPEConflictRemover::visit(StringGraph* pGraph, Vertex* pVertex)
{
    (void)pGraph;
    (void)pVertex;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pVertex->getEdges(dir);

        if(edges.size() > 1)
        {
            bool hasTrusted = false;
            for(size_t j = 0; j < edges.size(); ++j)
            {
                if(edges[j]->isTrusted)
                {
                    hasTrusted = true;
                }
            }
            
            if(hasTrusted)
            {
                for(size_t j = 0; j < edges.size(); ++j)
                {
                    if(!edges[j]->isTrusted)
                    {
                        edges[j]->setColor(GC_BLACK);
                        edges[j]->getTwin()->setColor(GC_BLACK);
                    }
                    if(edges[j]->getComp() == EC_SAME)
                        num_same++;
                    else
                        num_diff++;
                }
            }
        }    
    }
    return 0;
}

void SGPEConflictRemover::postvisit(StringGraph* pGraph)
{
    pGraph->sweepEdges(GC_BLACK);
    printf("Removed %d diff %d same\n", num_diff, num_same);
}

//
// SGPairedOverlapVisitor - print a formatted report to stdout
// detailing how much overlap there is between both end of a paired
// read.
//
bool SGPairedOverlapVisitor::visit(StringGraph* /*pGraph*/, Vertex* /*pVertex*/)
{
#if 0
    Vertex* pPairSV = pVertex->getPairVertex();
    if(pPairSV == NULL)
        return false;

    EdgePtrVec edges = pVertex->getEdges();

    // Determine which vertices that are paired to pVertex
    // have a pair that overlaps with pPairVertex
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pVWEdge = edges[i];
        Vertex* pW = pVWEdge->getEnd();
        Vertex* pPairW = pW->getPairVertex();
        if(pPairW == NULL)
            continue;

        EdgePtrVec ppw_edges = pPairW->findEdgesTo(pPairSV->getID());
        size_t overlap_len = pVWEdge->getMatchLength();

        if(pVWEdge->getComp() == EC_SAME)
        {
            if(ppw_edges.size() == 1)
            {
                Edge* pPPEdge = ppw_edges.front();
                size_t pair_overlap_len = pPPEdge->getMatchLength();
                printf("pairoverlap\t%s\t%s\t%zu\t%zu\n", pVertex->getID().c_str(), pW->getID().c_str(), overlap_len, pair_overlap_len);
            }
            else
            {
                printf("pairoverlap\t%s\t%s\t%zu\t%d\n", pVertex->getID().c_str(), pW->getID().c_str(), overlap_len, 0);
            }
        }
    }
#endif
    return false;
}

