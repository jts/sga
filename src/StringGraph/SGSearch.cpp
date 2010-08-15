//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGSearch - Algorithms and data structures
// for searching a string graph
//
#include "SGSearch.h"
#include <queue>

//
SGWalk::SGWalk(const Vertex* pStartVertex) : m_pStartVertex(pStartVertex)
{

}

//
void SGWalk::addEdge(Edge* pEdge)
{
    m_edges.push_back(pEdge);
}

//
Edge* SGWalk::getLastEdge() const
{
    assert(!m_edges.empty());
    return m_edges.back();
}

// 
void SGWalk::print() const
{
    std::cout << "Walk start: " << m_pStartVertex->getID() << "\nWalk: ";
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        std::cout << *(*iter) << " ";
    std::cout << "\n";
}   

//
struct SGSearchNode
{
    SGSearchNode(const Vertex* pStart) : walk(pStart), distance(0) {}
    SGWalk walk;
    int distance;
};

void SGSearch::findWalks(const Vertex* pX, const Vertex* pY, EdgeDir initialDir,
                         int maxDistance, SGWalkVector& outWalks)
{
    std::cout << "Resolving path from " << pX->getID() << " to " << pY->getID() << "\n";

    // Create the initial path nodes
    typedef std::queue<SGSearchNode> SearchQueue;
    EdgeDir dir = initialDir;
    SearchQueue queue;
    EdgePtrVec edges = pX->getEdges(dir);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        assert(!pEdge->getOverlap().isContainment());

        SGSearchNode node(pX);
        node.walk.addEdge(pEdge);
        node.distance = pEdge->getSeqLen();
        queue.push(node);
    }

    while(!queue.empty())
    {
        SGSearchNode& currNode = queue.front();
        // Check if we have found pY or exceeded the distance
        Edge* pWZ = currNode.walk.getLastEdge(); 
        Vertex* pZ = pWZ->getEnd();
        
        if(pZ == pY)
        {
            outWalks.push_back(currNode.walk);
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
                currNode.walk.addEdge(pNextEdge);

                // 
                for(size_t i = 1; i < zEdges.size(); ++i)
                {
                    SGSearchNode branch = currNode;
                    Edge* pBranchEdge = zEdges[i];
                    branch.distance += pBranchEdge->getSeqLen();
                    branch.walk.addEdge(pBranchEdge);
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

